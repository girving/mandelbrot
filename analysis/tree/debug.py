import os, sys
import numpy as np
from fractions import Fraction
exec(open('learn.py').read().split("\nif __name__ == '__main__':")[0])
rows, cols, H, Fcard, ext = build_hankel()
good_r = ~np.isnan(H).any(1); rows = [r for r, g in zip(rows, good_r) if g]; H = H[good_r]
keep_c = ~np.isnan(H).any(0); cols = [c for c, g in zip(cols, keep_c) if g]; H = H[:, keep_c]
cidx = {c: j for j, c in enumerate(cols)}; ridx = {r: i for i, r in enumerate(rows)}
Uu, Ss, Vt = np.linalg.svd(H, full_matrices=False)
K = 30; P = Uu[:, :K] * Ss[:K]; Q = Vt[:K]; Pp = np.linalg.pinv(P)
print('rank-K reconstruction error of H: %.1e' % np.abs(P @ Q - H).max())
# Digit closure for b = 1: residual of H[:, 1·s] = P N Q[:, s]
for b in (1, 2, 3):
    pairs = [(cidx[('dig', s)], cidx[('dig', (b,) + s)]) for s in S40 if ('dig', s) in cidx and ('dig', (b,) + s) in cidx]
    N = Pp @ H[:, [p[1] for p in pairs]] @ np.linalg.pinv(Q[:, [p[0] for p in pairs]])
    r = P @ N @ Q[:, [p[0] for p in pairs]] - H[:, [p[1] for p in pairs]]
    print('b=%d: %d column pairs, closure residual max %.1e; cond(Q[:, s]) %.1e; rows: root %.1e ctx %.1e' % (
        b, len(pairs), np.abs(r).max(), np.linalg.cond(Q[:, [p[0] for p in pairs]]),
        np.abs(r[[i for i, rr in enumerate(rows) if rr[0] == 'root']]).max(), np.abs(r[[i for i, rr in enumerate(rows) if rr[0] == 'ctx']]).max()))
pairs = [(cidx[('dig', cf(s2))], cidx[('sep', s2)]) for s2 in S20 if ('sep', s2) in cidx and ('dig', cf(s2)) in cidx]
Ns = Pp @ H[:, [p[1] for p in pairs]] @ np.linalg.pinv(Q[:, [p[0] for p in pairs]])
r = P @ Ns @ Q[:, [p[0] for p in pairs]] - H[:, [p[1] for p in pairs]]
print('#: %d pairs, closure residual max %.1e, cond %.1e' % (len(pairs), np.abs(r).max(), np.linalg.cond(Q[:, [p[0] for p in pairs]])))
# Row consistency: P[row u·b] vs P[row u] N(b) for root rows (row closure check)
pairs = [(ridx[('root', u)], ridx[('root', u + (1,))]) for u in U20 if ('root', u) in ridx and ('root', u + (1,)) in ridx]
N1 = Pp @ H[:, [cidx[('dig', (1,) + s)] for s in S40 if ('dig', (1,) + s) in cidx and ('dig', s) in cidx]] @ np.linalg.pinv(Q[:, [cidx[('dig', s)] for s in S40 if ('dig', (1,) + s) in cidx and ('dig', s) in cidx]])
d = P[[p[1] for p in pairs]] - P[[p[0] for p in pairs]] @ N1
print('row consistency P[u·1] - P[u] N(1): max %.1e (|P| max %.1e)' % (np.abs(d).max(), np.abs(P).max()))
# Separator row consistency: P[ctx r1, (u)] vs P[root u'] for u' = digits of r1 ... : P[('ctx', r1, u)] = P[('root', cf(r1))] N(#) N(u)?
for r1 in R1[:3]:
    if ('root', cf(r1)) not in ridx: continue
    for u in [(2,), (3,), (1, 2)]:
        if ('ctx', r1, u) not in ridx: continue
        v = P[ridx[('root', cf(r1))]] @ Ns
        Nu = {}
        print('  r1=%s u=%s: |P[ctx] - P[root r1] N(#) N(u)| = %s' % (r1, u, 'n/a'))
