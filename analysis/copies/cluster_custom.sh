#!/bin/bash
# Explicit-key family runs for the cluster (bench.sh RUNS line: bash analysis/copies/cluster_custom.sh), from the
# repository root after the build.  Waits for $IN/lines.txt ("name lo hi k" lines for limb_families custom: keys
# already specific to k, e.g. the two-transit sector of diag_census.py) and $IN/ready, finds centers from both rays,
# computes areas with bulb_batch, and writes $OUT/area.txt.gz; then holds HOLD seconds for kubectl cp.
set -euo pipefail
: "${IN:=/tmp/in}" "${OUT:=/tmp/out}" "${HOLD:=0}" "${MANDELBROT_THREADS:=2}"
mkdir -p "$IN" "$OUT"
t0=$(date +%s)
stamp() { echo "[$(( $(date +%s) - t0 )) s] $*"; }
stamp "waiting for $IN/ready"
until [ -f "$IN/ready" ]; do sleep 5; done
stamp "$(wc -l < "$IN/lines.txt") jobs"
./build/release/limb_families custom "$MANDELBROT_THREADS" < "$IN/lines.txt" > "$OUT/jobs.txt" 2> "$OUT/err.txt" || true
tail -1 "$OUT/err.txt"
./build/release/bulb_batch --out "$OUT/area.txt" < "$OUT/jobs.txt"
stamp "areas done"
gzip -f "$OUT/area.txt"
md5sum "$OUT"/*.gz
if [ "$HOLD" -gt 0 ]; then echo "holding $HOLD s for kubectl cp"; sleep "$HOLD"; fi
