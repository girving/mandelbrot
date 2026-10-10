"""Fit k^4 area(k) = C + Σ_{j≤m} b_j k^-j per family (bulb_batch output of family_asym.py jobs, keys NAME_k), in Decimal:
the limit C against the model, the coefficient growth |b_j|^(1/j) (the expansion's radius in 1/k) and the stability
of the b_j across fit windows.   python3 family_fit.py areas.txt comps.txt"""
import sys, math
from collections import defaultdict
from decimal import Decimal as D, getcontext
getcontext().prec = 80
fam = defaultdict(list)
for l in open(sys.argv[1]):
    r = l.split()
    if 'failed' in r: continue
    nm, k = r[0].rsplit('_', 1); k = int(k)
    fam[nm].append((k, (D(r[5]) + D(r[10])) * D(k) ** 4, abs(float(r[8]))))
model = {l.split()[0]: (int(l.split()[1]), float(l.split()[4])) for l in open(sys.argv[2])}
def lsq(M, y):
    n = len(M[0])
    N = [[sum(M[r][i] * M[r][j] for r in range(len(M))) for j in range(n)] + [sum(M[r][i] * y[r] for r in range(len(M)))] for i in range(n)]
    for i in range(n):
        p = max(range(i, n), key=lambda r: abs(N[r][i])); N[i], N[p] = N[p], N[i]
        for r in range(n):
            if r != i and N[r][i]:
                f = N[r][i] / N[i][i]; N[r] = [a - f * b for a, b in zip(N[r], N[i])]
    return [N[i][n] / N[i][i] for i in range(n)]
def fit(pts, m):
    M = [[D(1)] + [D(k) ** (-j) for j in range(1, m + 1)] for k, _, _ in pts]
    c = lsq(M, [g for _, g, _ in pts])
    res = max(abs(sum(a * b for a, b in zip(row, c)) - g) / g for row, (k, g, _) in zip(M, pts))
    return c, float(res)
for nm, pts in fam.items():
    pts.sort(); r, Cm = model[nm]
    print('%s (level %d, model C %.7e): k %d..%d, %d members, max conv %.1e' % (nm, r, Cm, pts[0][0], pts[-1][0], len(pts), max(p[2] for p in pts)))
    best = None
    for K0 in (4, 6, 9, 13, 19):
        for m in (8, 12, 16, 20):
            sel = [p for p in pts if p[0] >= K0]
            if len(sel) < m + 5: continue
            c, res = fit(sel, m)
            print('   K0 %2d m %2d  C %.15e (vs model %+.1e)  resid %.1e  b1/C %+.4f  b2/C %+.3f  b3/C %+.2f' % (
                K0, m, c[0], float(c[0]) / Cm - 1, res, float(c[1] / c[0]), float(c[2] / c[0]), float(c[3] / c[0])))
            if best is None or res < best[1]: best = (c, res, K0, m)
    c = best[0]
    print('   best (K0 %d m %d) |b_j/C|^(1/j):' % (best[2], best[3]), ' '.join('%.2f' % (abs(float(c[j] / c[0])) ** (1 / j)) for j in range(1, len(c))))
