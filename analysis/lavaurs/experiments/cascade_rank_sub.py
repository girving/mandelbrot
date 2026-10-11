"""Mass-weighted spectrum of the island-children matrix (cascade_rank.py --out file) over parent subsets, to see
whether the rank at a given precision saturates as parents are added.   python3 cascade_rank_sub.py single.txt out"""
import sys, math
import numpy as np
comps = []
for l in open(sys.argv[1]):
    f = l.split(); sat = f[0].startswith('sat'); n = int(f[0][3:]) if sat else int(f[0])
    comps.append(('B%d' % len(comps) if sat else 'F%d' % len(comps), n, float(f[3]), float(f[5]), sat))
info = {c[0]: c for c in comps}
M = {}; cols = {}; rows = {}
for l in open(sys.argv[2]):
    f = l.split(); nm, br, j = f[0].rsplit('|', 2); src, t = nm.rsplit('~', 1); w = float(f[3])
    if w >= 1: continue
    key = (t, int(j) + info[src][1] + 1); cols.setdefault(key, len(cols)); rows.setdefault(src, len(rows))
    M[(rows[src], cols[key])] = M.get((rows[src], cols[key]), 0) + w * info[t][3]
A = np.zeros((len(rows), len(cols)))
for (i, j), v in M.items(): A[i, j] = v
order = np.argsort([-info[s][2] for s in rows])      # parents by their own mass
A = A[order]
for n in (75, 150, len(A)):
    B = A[:n] / A[:n].sum()
    s = np.linalg.svd(B, compute_uv=False); s /= s[0]
    rk = lambda e: int((s > e).sum())
    print('%3d parents: rank at 1e-4/1e-6/1e-8/1e-10: %d %d %d %d   tail %s' % (n, rk(1e-4), rk(1e-6), rk(1e-8), rk(1e-10), ' '.join('%.1e' % x for x in s[20:60:4])))
# also: mass captured by the top-k parents, to show the matrix covers the mass
m = A.sum(1); print('parents covering 1 - 1e-4 / 1e-6 of this block\'s mass: %d %d of %d' % (np.searchsorted(np.cumsum(m) / m.sum(), 1 - 1e-4), np.searchsorted(np.cumsum(m) / m.sum(), 1 - 1e-6), len(m)))
