"""Approach C on the one-level cascade: split the island-children matrix (cascade_rank.py --out, with its .err island
data) into a near field, |y - κ_c| ≈ |κ''/2| |σ_child - σ_c|^2 < d0 (the parent's critical value close to the target),
and the far field; report the far part's mass-weighted rank and the near part's size and sparsity.

  python3 cascade_split.py single.txt ich.txt"""
import re, sys
import numpy as np
comps = []
for l in open(sys.argv[1]):
    f = l.split(); sat = f[0].startswith('sat'); n = int(f[0][3:]) if sat else int(f[0])
    comps.append(('B%d' % len(comps) if sat else 'F%d' % len(comps), n, float(f[3]), float(f[5]), sat))
info = {c[0]: c for c in comps}
isl = {}
for l in open(sys.argv[2] + '.err'):
    m = re.match(r"(\S+)~\S+: σ_c (\S+?)([+-][0-9.e+-]+)i .* κ'' (\S+)", l)
    if m: isl[m.group(1)] = (complex(float(m.group(2)), float(m.group(3))), float(m.group(4)))
ent = []   # (src, col, mass, q)
for l in open(sys.argv[2]):
    f = l.split(); nm, br, j = f[0].rsplit('|', 2); src, t = nm.rsplit('~', 1); w = float(f[3])
    if w >= 1 or src not in isl: continue
    sc, k2 = isl[src]; sg = complex(float(f[1]), float(f[2]))
    ent.append((src, (t, int(j) + info[src][1] + 1), w * info[t][3], 0.5 * k2 * abs(sg - sc) ** 2))
rows = {}; cols = {}
for s, c, m, q in ent: rows.setdefault(s, len(rows)); cols.setdefault(c, len(cols))
tot = sum(e[2] for e in ent)
for d0 in (1e9, 1.0, 0.3, 0.1, 0.03, 0.01):
    A = np.zeros((len(rows), len(cols))); nn = 0; nm_ = 0.0; ntar = set()
    for s, c, m, q in ent:
        if q >= d0: A[rows[s], cols[c]] += m
        else: nn += 1; nm_ += m; ntar.add(c[0])
    sv = np.linalg.svd(A / tot, compute_uv=False); sv /= sv[0]
    rk = lambda e: int((sv > e).sum())
    print('d0 %-6g near: %6d entries (%5.2f%%), mass share %.2e, %d targets;  far rank at 1e-4/1e-5/1e-6/1e-8: %d %d %d %d' % (
        d0, nn, 100 * nn / len(ent), nm_ / tot, len(ntar), rk(1e-4), rk(1e-5), rk(1e-6), rk(1e-8)))
