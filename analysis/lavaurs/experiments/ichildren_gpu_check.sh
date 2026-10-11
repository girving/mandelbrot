#!/bin/bash
# GPU vs CPU island children (glavaurs ichildren): identical jobs on the GPU and on CPU threads, timed, compared.
#   bash ichildren_gpu_check.sh single.txt [sources] [targets] [threads]
set -e
SINGLE=$1; NS=${2:-60}; NT=${3:-40}; TH=${4:-6}
B=build/release/glavaurs
T=$(mktemp -d)
python3 - "$SINGLE" "$NS" "$NT" > $T/jobs.txt <<'PY'
import sys
single, ns, nt = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
comps = []
for l in open(single):
    f = l.split(); sat = f[0].startswith('sat'); n = int(f[0][3:]) if sat else int(f[0])
    comps.append(('B%d' % len(comps) if sat else 'F%d' % len(comps), n, complex(float(f[1]), float(f[2])), float(f[3]), sat))
cut = (-0.16135210336424927 - 4.1583377953217475) / 2
own = [c for c in comps if c[2].imag > cut]
targets = [c for c in own if c[4]] + sorted([c for c in own if not c[4]], key=lambda c: -c[3])[:nt]
for nm, n, s, C, _ in sorted([c for c in own if not c[4]], key=lambda c: -c[3])[:ns]:
    for tn, ntt, ct, _, _ in targets:
        print('%s~%s 2 %.17g %.17g %d %.17g %.17g %d %d' % (nm, tn, s.real, s.imag, ntt, ct.real, ct.imag, -(n + 1) - 20, -(n + 1) + 4))
PY
echo "jobs: $(wc -l < $T/jobs.txt) (paths: 2 x 25 shifts each)"
S0=$(date +%s.%N); GL_SIDE=-1 GL_GPU=1 ICH_DEBUG=1 $B 1 2 ichildren $TH < $T/jobs.txt > $T/gpu.txt 2> $T/gpu.err; S1=$(date +%s.%N)
GL_SIDE=-1 GL_GPU=0 $B 1 2 ichildren $TH < $T/jobs.txt > $T/cpu.txt 2> $T/cpu.err; S2=$(date +%s.%N)
grep "ichildren:" $T/gpu.err || true
echo "GPU $(echo "$S1 - $S0" | bc) s, CPU ($TH threads) $(echo "$S2 - $S1" | bc) s"
python3 - $T/gpu.txt $T/cpu.txt <<'PY'
import sys
g = {l.split()[0]: l.split()[1:] for l in open(sys.argv[1])}; c = {l.split()[0]: l.split()[1:] for l in open(sys.argv[2])}
both = set(g) & set(c)
ident = sum(g[k] == c[k] for k in both)
ds = max((abs(complex(float(g[k][0]), float(g[k][1])) - complex(float(c[k][0]), float(c[k][1]))) for k in both), default=0)
dw = max((abs(float(g[k][2]) / float(c[k][2]) - 1) for k in both if float(c[k][2]) < 1), default=0)
print('GPU %d children, CPU %d, common %d, bit-identical %d; only GPU %d, only CPU %d; max |Δσ| %.2e, max |Δw/w| %.2e' % (
    len(g), len(c), len(both), ident, len(set(g) - set(c)), len(set(c) - set(g)), ds, dw))
PY
