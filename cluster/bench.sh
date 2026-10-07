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
ENV=/data/env-v3
if [ ! -x $ENV/bin/clang++ ]; then
  step "Creating conda environment $ENV"
  timeout 1800 /data/bin/micromamba create -y -q -p $ENV -c conda-forge --override-channels \
    clangxx=22 libcxx=22 libcxx-devel=22 llvm-openmp=22 lld=22 meson ninja pkg-config \
    "libflint>=3" gmp libpng zlib "openssl>=3" git curl ca-certificates
fi
export PATH=$ENV/bin:/usr/local/cuda/bin:$PATH
export PKG_CONFIG_PATH=$ENV/lib/pkgconfig
export CPATH=$ENV/include
export LIBRARY_PATH=$ENV/lib:/usr/local/cuda/lib64
export LD_LIBRARY_PATH=$ENV/lib:/usr/local/cuda/lib64
export CXX=clang++ CUDA_PATH=/usr/local/cuda

step "Machine"
nvidia-smi --query-gpu=name,memory.total,driver_version,clocks.max.sm --format=csv 2>/dev/null || echo "no GPU"
echo "CPUs: $(nproc) visible, using $MANDELBROT_THREADS threads"
clang++ --version | head -1

step "Cloning $COMMIT"
cd /tmp
# GitHub sometimes refuses anonymous clones from the cluster for a while: retry, then fall back to a tarball
cloned=
for attempt in 1 2 3 4; do
  if timeout 300 git clone -q --depth 1 --branch "$COMMIT" https://github.com/girving/mandelbrot; then cloned=1; break; fi
  rm -rf mandelbrot; echo "clone attempt $attempt failed"; sleep $((30 * attempt))
done
if [ -z "$cloned" ]; then
  mkdir mandelbrot
  timeout 300 curl -fsSL "https://codeload.github.com/girving/mandelbrot/tar.gz/$COMMIT" \
    | tar -xz -C mandelbrot --strip-components=1 || { echo "tarball download failed"; exit 1; }
  echo "fetched $COMMIT as a tarball (no git history)"
fi
cd mandelbrot
echo "commit $(git rev-parse HEAD 2>/dev/null || echo "$COMMIT (tarball)")"

step "Configuring and building"
CXXFLAGS="-O3 -march=native" timeout 600 meson setup build/release --buildtype=release > /tmp/setup.log \
  || { cat /tmp/setup.log; exit 1; }
grep -i -E "cuda|openmp" /tmp/setup.log || true
timeout 1800 meson compile -C build/release tree_test escape_tree orbit_bench orbit_census orbit_regimes ${EXTRA_TARGETS:-}
if [ -n "${BUILD64:-}" ]; then  # 64-bit orbit step counts, for max_iter beyond 2^30
  CXXFLAGS="-O3 -march=native -DMANDELBROT_ORBIT64" timeout 600 meson setup build/release64 --buildtype=release \
    > /tmp/setup64.log || { cat /tmp/setup64.log; exit 1; }
  timeout 1800 meson compile -C build/release64 escape_tree
fi

mkdir -p /data/results
OUT=/data/results/bench-$(date +%Y%m%d-%H%M%S).txt
run() { step "$*"; timeout "${RUN_TIMEOUT:-1800}" "$@" 2>&1 | tee -a "$OUT"; }

run ./build/release/tree_test
if [ -n "${SHARDED:-}" ]; then
  # One escape_tree run split over the node's GPUs: SHARDED is the command (options before the thresholds),
  # SHARDS the number of GPUs, and RUN_NAME a directory under /data/results for shard results, kept across
  # jobs so that a rerun skips shards already done.  Each shard gets its share of the CPUs for its CPU tail.
  # With SHARD set, run only that shard (on this pod's GPU, with all its CPUs), leaving the merge to a
  # later job without SHARD once every shard is saved (cluster/idle_shards.py).
  : "${SHARDS:=8}" "${RUN_NAME:?RUN_NAME names the sharded run}"
  DIR=/data/results/$RUN_NAME
  mkdir -p "$DIR"
  set -- $SHARDED
  prog=$1; shift
  step "Sharded run $RUN_NAME: $SHARDS shards of $SHARDED"
  pids=()
  first=0 last=$((SHARDS - 1)) gpus=$SHARDS
  if [ -n "${SHARD:-}" ]; then first=$SHARD last=$SHARD gpus=1; fi
  for ((s = first; s <= last; s++)); do
    if [ -s "$DIR/shard-$s.txt" ]; then echo "shard $s already done"; continue; fi
    CUDA_VISIBLE_DEVICES=$((s - first)) MANDELBROT_THREADS=$(( MANDELBROT_THREADS / gpus )) timeout "${RUN_TIMEOUT:-1800}" \
      $prog --shard $s/$SHARDS --save "$DIR/shard-$s.tmp" "$@" > "$DIR/shard-$s.log" 2>&1 \
      && mv "$DIR/shard-$s.tmp" "$DIR/shard-$s.txt" &
    pids+=($!)
  done
  failed=0
  for pid in ${pids[@]+"${pids[@]}"}; do wait "$pid" || failed=1; done
  for ((s = first; s <= last; s++)); do
    echo "--- shard $s"
    # A single shard's whole log (progress lines included); otherwise the end of each
    if [ -n "${SHARD:-}" ]; then cat "$DIR/shard-$s.log" || true; else tail -n 12 "$DIR/shard-$s.log" || true; fi
  done
  [ $failed = 0 ] || { echo "some shards failed; rerun to finish them"; exit 1; }
  if [ -n "${SHARD:-}" ]; then
    # The saved result, for runs whose /data does not outlive the job (idle_shards.py "storage": "local")
    echo "=== RESULT BEGIN $RUN_NAME shard $SHARD"
    cat "$DIR/shard-$SHARD.txt"
    echo "=== RESULT END"
    step "Done with shard $SHARD of $RUN_NAME"
    exit 0
  fi
  files=$(IFS=,; f=(); for ((s = 0; s < SHARDS; s++)); do f+=("$DIR/shard-$s.txt"); done; echo "${f[*]}")
  run $prog --merge "$files" "$@"
elif [ -n "${RUNS:-}" ]; then
  # Custom runs: RUNS holds one command per line, optionally prefixed by VAR=value settings
  while IFS= read -r cmd; do
    [ -n "$cmd" ] && run env $cmd
  done <<< "$RUNS"
else
  run ./build/release/escape_tree --cuda 16384 1048568                 # Laptop config: 3.4e-7 in 44 s on M5 Pro
  run ./build/release/escape_tree 16384 1048568                        # Same on this pod's CPUs
  run ./build/release/escape_tree --cuda --depth 6 16384 1048568       # 4x the leaves
fi
step "Done; results in $OUT"
