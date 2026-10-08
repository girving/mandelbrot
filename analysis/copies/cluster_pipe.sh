#!/bin/bash
# Family pipeline for the cluster (bench.sh RUNS line: bash analysis/copies/cluster_pipe.sh), from the repository root
# after the build.  Waits for $IN/lines.txt and then $IN/ready (copy them in with kubectl cp), runs the size engine on
# the lines (limb_families size), converts its centers to `bulb_batch --local` jobs (c = -3/4 + δ by exact two-sums,
# as size_to_local.py), computes the areas, and writes $OUT/size.txt.gz and $OUT/area.txt.gz; then holds HOLD seconds
# for kubectl cp.
set -euo pipefail
: "${IN:=/tmp/in}" "${OUT:=/tmp/out}" "${HOLD:=0}" "${MANDELBROT_THREADS:=2}"
mkdir -p "$IN" "$OUT"
t0=$(date +%s)
stamp() { echo "[$(( $(date +%s) - t0 )) s] $*"; }
stamp "waiting for $IN/ready"
until [ -f "$IN/ready" ]; do sleep 5; done
stamp "$(wc -l < "$IN/lines.txt") families"
./build/release/limb_families size "$MANDELBROT_THREADS" < "$IN/lines.txt" > "$OUT/size.txt" 2> "$OUT/size_err.txt" || true
tail -2 "$OUT/size_err.txt"
stamp "$(wc -l < "$OUT/size.txt") centers"
awk 'function two_sum(a, b) { S = a + b; BB = S - a; E = (a - (S - BB)) + (b - BB) }
  { two_sum(-0.75, $2); h = S; l = E + $3; two_sum(h, l); printf "%s %s %.17g %.17g %.17g %.17g\n", $1, $6, S, E, $4, $5 }' \
  "$OUT/size.txt" > "$OUT/jobs.txt"
./build/release/bulb_batch --local --out "$OUT/area.txt" < "$OUT/jobs.txt"
stamp "areas done"
gzip -f "$OUT/size.txt" "$OUT/area.txt"
md5sum "$OUT"/*.gz
if [ "$HOLD" -gt 0 ]; then echo "holding $HOLD s for kubectl cp"; sleep "$HOLD"; fi
