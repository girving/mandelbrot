#!/bin/bash
# Build mandelbrot with CUDA on an H200 pod and benchmark escape_tree (run by cluster/bench.yaml).
#
# The conda-forge build environment lives on the PVC mounted at /data and is created on first use
# (a few minutes); the code is cloned from GitHub at $COMMIT onto node-local disk.  Results are also
# written to /data/results, but the log (kubectl logs) is the record.
set -euo pipefail
: "${COMMIT:=hybrid-area}"
: "${MANDELBROT_THREADS:=22}"
export MANDELBROT_THREADS
T0=$(date +%s)
step() { echo; echo "=== $* [$(( $(date +%s) - T0 )) s]"; }

# Build environment
export MAMBA_ROOT_PREFIX=/data/mamba
ENV=/data/env-v2
if [ ! -x $ENV/bin/clang++ ]; then
  step "Creating conda environment $ENV"
  timeout 1800 /data/bin/micromamba create -y -q -p $ENV -c conda-forge --override-channels \
    clangxx=22 libcxx=22 libcxx-devel=22 llvm-openmp=22 lld=22 meson ninja pkg-config \
    "libflint>=3" gmp libpng zlib "openssl>=3" git curl
fi
export PATH=$ENV/bin:/usr/local/cuda/bin:$PATH
export PKG_CONFIG_PATH=$ENV/lib/pkgconfig
export CPATH=$ENV/include
export LIBRARY_PATH=$ENV/lib:/usr/local/cuda/lib64
export LD_LIBRARY_PATH=$ENV/lib:/usr/local/cuda/lib64
export CXX=clang++ CUDA_PATH=/usr/local/cuda

step "Machine"
nvidia-smi --query-gpu=name,memory.total,driver_version,clocks.max.sm --format=csv
echo "CPUs: $(nproc) visible, using $MANDELBROT_THREADS threads"
clang++ --version | head -1

step "Cloning $COMMIT"
cd /tmp
timeout 300 git clone -q --depth 1 --branch "$COMMIT" https://github.com/girving/mandelbrot
cd mandelbrot
echo "commit $(git rev-parse HEAD)"

step "Configuring and building"
CXXFLAGS="-O3 -march=native" timeout 600 meson setup build/release --buildtype=release > /tmp/setup.log \
  || { cat /tmp/setup.log; exit 1; }
grep -i -E "cuda|openmp" /tmp/setup.log || true
timeout 1800 meson compile -C build/release escape_batch_test escape_tree

mkdir -p /data/results
OUT=/data/results/bench-$(date +%Y%m%d-%H%M%S).txt
run() { step "$*"; timeout 1800 "$@" 2>&1 | tee -a "$OUT"; }

run ./build/release/escape_batch_test
if [ -n "${RUNS:-}" ]; then
  # Custom escape_tree runs: RUNS holds one argument list per line
  while IFS= read -r args; do
    [ -n "$args" ] && run ./build/release/escape_tree $args
  done <<< "$RUNS"
else
  run ./build/release/escape_tree --cuda 16384 1048568                 # Laptop config: 3.4e-7 in 44 s on M5 Pro
  run ./build/release/escape_tree 16384 1048568                        # Same on this pod's CPUs
  run ./build/release/escape_tree --cuda --depth 6 16384 1048568       # 4x the leaves
fi
step "Done; results in $OUT"
