#!/bin/bash
# Lavaurs-model component areas on the cluster (bench.sh RUNS line: bash analysis/lavaurs/cluster_model.sh), from the
# repository root after the build.  Waits for $IN/jobs.txt and $IN/ready (copy them in with kubectl cp), runs
# lavaurs_area, writes $OUT/model.txt.gz, then holds HOLD seconds for kubectl cp.  RUN_SECONDS (if set) stops
# lavaurs_area early; --island output streams per pair, so a stopped run keeps its finished pairs.
set -euo pipefail
: "${IN:=/tmp/in}" "${OUT:=/tmp/out}" "${HOLD:=0}" "${MANDELBROT_THREADS:=2}" "${RUN_SECONDS:=0}"
mkdir -p "$IN" "$OUT"
t0=$(date +%s)
echo "[0 s] waiting for $IN/ready"
until [ -f "$IN/ready" ]; do sleep 5; done
echo "[$(( $(date +%s) - t0 )) s] $(wc -l < "$IN/jobs.txt") jobs"
timeout "$RUN_SECONDS" ./build/release/lavaurs_area ${MODE:-} "$MANDELBROT_THREADS" < "$IN/jobs.txt" > "$OUT/model.txt" \
  || echo "lavaurs_area stopped (status $?)"
echo "[$(( $(date +%s) - t0 )) s] done"
gzip -f "$OUT/model.txt"
md5sum "$OUT"/*.gz
if [ "$HOLD" -gt 0 ]; then echo "holding $HOLD s for kubectl cp"; sleep "$HOLD"; fi
