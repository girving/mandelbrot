#!/bin/bash
# The general-root strip census on one machine (a GPU pod via cluster/bench.sh RUNS, or a laptop):
#   bash analysis/lavaurs/gl_pipeline.sh p q gate theta levels outdir [check]
# single transits (gl_single.py), then per level r = 2..levels: children (gl_census.py; on the GPU when there is one),
# tunings (gl_tune.py), limbs (gl_limbsum.py; CLASSIFY_MIN: classify only components above it), with timings.  With
# "check": first a small census on CPU and GPU, compared.  GL_START: resume at that level from $out's files.  GL_SOURCES, GL_TARGETS: census sizes (default 300, 100;
# not TARGETS, which bench.sh uses for build targets).  GL_NO_CENSUS=1: redo only tunings and limbs from $out's census files.
set -euo pipefail
p=$1 q=$2 gate=$3 theta=$4 levels=$5 out=$6 check=${7:-}
T=${MANDELBROT_THREADS:-8}
D=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$out"
t0=$(date +%s); stamp() { echo "[$(( $(date +%s) - t0 )) s] $*"; }
start=${GL_START:-2}   # resume at this level from $out (its earlier levels' files must be there)
if [ $start = 2 ]; then
  stamp "single transits"
  python3 "$D/gl_single.py" $p $q $gate --cmin 1e-13 --exact 1e-11 --threads $T --out "$out/single.txt" | head -4
fi
if [ -n "$check" ]; then
  stamp "GPU vs CPU check (20 sources x 20 targets, level 2)"
  for g in 0 1; do
    GL_GPU=$g python3 "$D/gl_census.py" $p $q $gate "$out/single.txt" --levels 2 --sources 20 --targets 20 --threads $T \
      --starts 26 --big-src 1e-9 --out "$out/check$g" > "$out/check$g.log" 2>&1
    stamp "  GL_GPU=$g done: $(grep 'primitive mass' "$out/check$g.log" | head -1)"
  done
  python3 "$D/gl_gpucheck.py" "$out/check0_raw_r2.txt" "$out/check1_raw_r2.txt"
fi
for ((r = start; r <= levels; r++)); do
  stamp "level $r census"
  if [ -n "${GL_NO_CENSUS:-}" ]; then
    :
  elif [ $r = 2 ]; then
    python3 "$D/gl_census.py" $p $q $gate "$out/single.txt" --levels 2 --threads $T --starts 26 --big-src 1e-9 --out "$out/c" \
      --sources ${GL_SOURCES:-300} --targets ${GL_TARGETS:-100}
  else
    python3 "$D/gl_census.py" $p $q $gate "$out/single.txt" --levels $r --resume $r --threads $T --starts 26 --big-src 1e-9 \
      --src-min ${GL_SRC_MIN:-1e-10} --out "$out/c" --targets ${GL_TARGETS:-100}
  fi
  stamp "level $r tunings"
  python3 "$D/gl_tune.py" $p $q $gate "$out/single.txt" "$out/c" --level $r --threads $T | head -12
  stamp "level $r limbs (tunings removed)"
  grep -v -F -f <(awk '{print $1" "}' "$out/c_tuned_r$r.txt") "$out/c_r$r.txt" > "$out/c_untuned_r$r.txt" || true
  python3 "$D/gl_limbsum.py" $p $q $gate $theta "$out/c_untuned_r$r.txt" --threads $T --cmin ${CLASSIFY_MIN:-0}
done
stamp "done"
