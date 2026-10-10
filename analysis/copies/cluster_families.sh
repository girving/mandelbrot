#!/bin/bash
# Seahorse-family runs for the cluster (cluster/bench.sh RUNS line: bash analysis/copies/cluster_families.sh), from the
# repository root after the build.  Writes $OUT/census.txt.gz (limb_families jobs J from limb K0 at k in KS, through
# bulb_batch) and $OUT/grid.txt.gz (the two-pass sequences of two_pass.py at every m ≤ M and m < k ≤ KMAX), then
# holds for HOLD seconds so the results can be copied out with kubectl cp.
set -euo pipefail
: "${OUT:=/tmp/families}" "${J:=14}" "${K0:=7}" "${KS:=16,20,24,28,32,40,48,56,64}" "${M:=64}" "${KMAX:=128}"
: "${HOLD:=0}" "${MANDELBROT_THREADS:=2}"
mkdir -p "$OUT"
T=$MANDELBROT_THREADS
t0=$(date +%s)
stamp() { echo "[$(( $(date +%s) - t0 )) s] $*"; }

stamp "census: J=$J from limb $K0 at k = $KS"
K0=$K0 ./build/release/limb_families jobs "$J" "$KS" "$T" > "$OUT/census_jobs.txt" 2> "$OUT/census_err.txt"
tail -1 "$OUT/census_err.txt"
./build/release/bulb_batch --out "$OUT/census.txt" < "$OUT/census_jobs.txt"
stamp "census done"

# Two-pass sequences (two_pass.py's SEQ): name, lo and hi key words for m, and all k with m < k ≤ KMAX
awk -v M="$M" -v KMAX="$KMAX" 'function rep(s, n,   r, i) { r = ""; for (i = 0; i < n; i++) r = r s; return r }
  BEGIN {
    for (m = 1; m <= M; m++) {
      ks = ""; for (k = m + 1; k <= KMAX; k++) ks = ks (ks == "" ? "" : ",") k
      print "E0m" m, "00" rep("10", m) "1", "00" rep("10", m - 1) "110", ks
      print "E2m" m, "0011" rep("01", m - 1) "001", "0011" rep("01", m - 1) "010", ks
      print "O1m" m, "0011" rep("01", m), "0011" rep("01", m - 1) "10", ks
      print "OLm" m, "010" rep("01", m - 1) "001", "010" rep("01", m - 1) "010", ks
    }
  }' > "$OUT/grid_lines.txt"
stamp "grid: $(wc -l < "$OUT/grid_lines.txt") sequences x k"
./build/release/limb_families custom "$T" < "$OUT/grid_lines.txt" > "$OUT/grid_jobs.txt" 2> "$OUT/grid_err.txt"
tail -1 "$OUT/grid_err.txt"
grep -c disagree "$OUT/grid_err.txt" || true
./build/release/bulb_batch --out "$OUT/grid.txt" < "$OUT/grid_jobs.txt"
stamp "grid done"

grep -c failed "$OUT/census.txt" "$OUT/grid.txt" || true
gzip -f "$OUT/census.txt" "$OUT/grid.txt" "$OUT/census_err.txt" "$OUT/grid_err.txt"
ls -la "$OUT"
md5sum "$OUT"/*.gz
if [ "$HOLD" -gt 0 ]; then stamp "holding $HOLD s for kubectl cp"; sleep "$HOLD"; fi
