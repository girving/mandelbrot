#!/bin/bash
# Deep two-pass size estimates for the cluster (bench.sh RUNS line: bash analysis/copies/cluster_deep.sh), from the
# repository root after the build: the sequences of two_pass.py at m in MS, k dense from m+1 to m+16 then geometric
# (×1.05) to 16m (two_pass_deep.py's ks up to rounding of halves), through `limb_families size`.  Writes $OUT/deep.txt.gz, then holds HOLD s
# for kubectl cp.
set -euo pipefail
: "${OUT:=/tmp/deep}" "${HOLD:=0}" "${MANDELBROT_THREADS:=2}"
: "${MS:=1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 20 24 28 32 40 48 56 64 80 96 112 128 160 192 224 256 320 384 448 512 640 768 1024}"
mkdir -p "$OUT"
t0=$(date +%s)
# Largest m first, so the slowest families start first
for m in $MS; do echo "$m"; done | sort -rn | awk 'function rep(s, n,   r, i) { r = ""; for (i = 0; i < n; i++) r = r s; return r }
  function round(x) { return int(x + 0.5) }
  {
    m = $1; ks = ""
    for (k = m + 1; k <= m + 16; k++) ks = ks (ks == "" ? "" : ",") k
    k = m + 16
    while (k < 16 * m) { k2 = round(k * 1.05); k = (k2 > k + 1 ? k2 : k + 1); ks = ks "," k }
    print "E0m" m, "00" rep("10", m) "1", "00" rep("10", m - 1) "110", ks
    print "E2m" m, "0011" rep("01", m - 1) "001", "0011" rep("01", m - 1) "010", ks
    print "O1m" m, "0011" rep("01", m), "0011" rep("01", m - 1) "10", ks
    print "OLm" m, "010" rep("01", m - 1) "001", "010" rep("01", m - 1) "010", ks
  }' > "$OUT/lines.txt"
echo "[$(( $(date +%s) - t0 )) s] $(wc -l < "$OUT/lines.txt") families"
./build/release/limb_families size "$MANDELBROT_THREADS" < "$OUT/lines.txt" > "$OUT/deep.txt" 2> "$OUT/err.txt" || true
tail -3 "$OUT/err.txt"
echo "[$(( $(date +%s) - t0 )) s] $(wc -l < "$OUT/deep.txt") lines"
gzip -f "$OUT/deep.txt"
md5sum "$OUT/deep.txt.gz"
if [ "$HOLD" -gt 0 ]; then echo "holding $HOLD s for kubectl cp"; sleep "$HOLD"; fi
