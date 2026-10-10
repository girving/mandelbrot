"""Compare two `glavaurs children` outputs (CPU and GPU): the same children (name, position) and weights.

  python3 gl_gpucheck.py cpu_raw.txt gpu_raw.txt"""
import sys
from collections import defaultdict

def load(f):
    d = defaultdict(list)
    for l in open(f):
        x = l.split(); d[x[0]].append((complex(float(x[1]), float(x[2])), float(x[3]), x[4]))
    return d

a, b = load(sys.argv[1]), load(sys.argv[2])
na, nb = sum(map(len, a.values())), sum(map(len, b.values()))
matched, dpos, dw, only_a, only_b = 0, 0.0, 0.0, 0, 0
for k in set(a) | set(b):
    A, B = a.get(k, []), b.get(k, [])
    used = set()
    for s, w, sat in A:
        best = min(range(len(B)), key=lambda i: abs(B[i][0] - s)) if B else None
        if best is None or abs(B[best][0] - s) > 1e-8 or best in used: only_a += 1; continue
        used.add(best); matched += 1
        dpos = max(dpos, abs(B[best][0] - s)); dw = max(dw, abs(B[best][1] / w - 1))
    only_b += len(B) - len(used)
print('GPU vs CPU children: %d CPU, %d GPU, %d matched (max |Δσ| %.1e, max |Δw/w| %.1e), only CPU %d, only GPU %d' % (
    na, nb, matched, dpos, dw, only_a, only_b))
