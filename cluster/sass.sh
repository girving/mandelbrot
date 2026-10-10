#!/bin/bash
# SASS of the kernels in a binary whose names match a pattern: bash cluster/sass.sh <binary> <pattern>
set -euo pipefail
/usr/local/cuda/bin/cuobjdump -sass "$1" | awk -v pat="$2" '
  /Function : / { keep = index($0, pat) > 0 }
  keep { gsub(/ *\/\* 0x[0-9a-f]+ \*\/ *$/, ""); gsub(/^ +/, ""); if ($0 != "") print }'
