"""In-model m-asymptotics: the Farey bulbs B_m of the q = 2 strip (60-digit model constants from glavaurs hp, lines
"m center_re center_im C ..."): fit m^4 C(m) = a_0 + Σ_{j≤J} a_j m^-j (optionally + Σ (log m) m^-j), in Decimal, for
m ≥ M0; the coefficient growth and the residual decide whether the in-model families are as analytic as the k ones.
  python3 farey_m_fit.py farey_hp.txt"""
import sys
from decimal import Decimal as D, getcontext
getcontext().prec = 80
pts = []
for l in open(sys.argv[1]):
    f = l.split()
    if len(f) >= 4: pts.append((int(f[0]), D(f[3]) * D(int(f[0])) ** 4))
def lsq(M, y):
    n = len(M[0])
    N = [[sum(M[r][i] * M[r][j] for r in range(len(M))) for j in range(n)] + [sum(M[r][i] * y[r] for r in range(len(M)))] for i in range(n)]
    for i in range(n):
        p = max(range(i, n), key=lambda r: abs(N[r][i])); N[i], N[p] = N[p], N[i]
        for r in range(n):
            if r != i and N[r][i]:
                f = N[r][i] / N[i][i]; N[r] = [a - f * b for a, b in zip(N[r], N[i])]
    return [N[i][n] / N[i][i] for i in range(n)]
for logs in (0, 1, 2):
    for M0 in (4, 8, 12):
        for J in (6, 10, 14):
            sel = [(m, g) for m, g in pts if m >= M0]
            def row(m):
                mm = D(m); r = [D(1)] + [mm ** (-j) for j in range(1, J + 1)]
                if logs: r += [mm.ln() ** e * mm ** (-j) for e in range(1, logs + 1) for j in range(1, J // 2 + 1)]
                return r
            M = [row(m) for m, _ in sel]
            if len(sel) < len(M[0]) + 3: continue
            c = lsq(M, [g for _, g in sel])
            res = max(abs(sum(a * b for a, b in zip(r, c)) - g) / g for r, (m, g) in zip(M, sel))
            print('logs %d M0 %2d J %2d: a0 %.15f  resid %.1e  a1 %+.5f a2 %+.4f  |a_j/a0|^(1/j) %s' % (
                logs, M0, J, c[0], res, c[1], c[2], ' '.join('%.2f' % (abs(float(c[j] / c[0])) ** (1 / j)) for j in range(1, J + 1))))
