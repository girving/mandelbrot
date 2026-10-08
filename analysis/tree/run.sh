#!/bin/bash
# Run bulb_areas for every parent in jobs_$1.txt (double-double), one parent at a time on all threads
D=$BULB_DATA/tree/ext
while read name P cr ci; do
  [ -s $D/$name.out ] && continue
  BULB_EXP=1 BULB_TOL=1e-9 MANDELBROT_THREADS=10 timeout 3600 $(dirname $0)/../../build/release/bulb_areas $P $cr $ci 64 < $D/$name.in > $D/$name.out
done < $D/jobs_$1.txt
