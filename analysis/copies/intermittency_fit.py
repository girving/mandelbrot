"""Fits for the intermittency families (intermittency.py): parents, and their satellite children.

  python3 intermittency.py 120 > jobs; bulb_batch --out parents < jobs
  (children: one line "key P c_re c_im p q" per parent line of parents and p/q with q ≤ 6)
  bulb_batch --out kids < kid_jobs
  python3 intermittency_fit.py parents [kids]

Parents: the phase σ_w = lim π/(7√d) − k (d = c + 1.75 ≈ π²/(49 (k + σ)²)), the expansion of a(k)·(k+σ)^6 in
1/(k+σ) (fit through k = 22..100), and Σ_k a(k) with the tail from that expansion.  Children: limits of the
area ratio child/parent in 1/(k+σ) (ratios, not F: F's normalization w is computed in double, which is noisy at
1e-9 for these ~1e-10-sized parents)."""
import sys
from decimal import Decimal as D, getcontext
getcontext().prec = 50
PI = D('3.14159265358979323846264338327950288419716939937510')

def neville(xs, ys):
    p = list(ys)
    for m in range(1, len(xs)):
        p = [(xs[i+m] * p[i] - xs[i] * p[i+1]) / (xs[i+m] - xs[i]) for i in range(len(p) - 1)]
    return p[0]

def fit(xs, ys):  # polynomial coefficients through the points
    n = len(xs); A = [[x**j for j in range(n)] + [y] for x, y in zip(xs, ys)]
    for i in range(n):
        piv = max(range(i, n), key=lambda r: abs(A[r][i])); A[i], A[piv] = A[piv], A[i]
        for r in range(n):
            if r != i:
                f = A[r][i] / A[i][i]; A[r] = [a - f * b for a, b in zip(A[r], A[i])]
    return [A[i][n] / A[i][i] for i in range(n)]

area = lambda r: D(float(r[5])) + D(float(r[10]))
fam, par = {}, {}
for l in open(sys.argv[1]):
    r = l.split(); par[r[0]] = area(r)
    fam.setdefault(r[0][0], {})[int(r[0][1:])] = (area(r), D(float(r[3])) + D('1.75'))
sig = {}
for f, F in sorted(fam.items()):
    s = {k: PI / (7 * F[k][1].sqrt()) - k for k in F}
    ks = [100 - 6 * i for i in range(8)]
    sig[f] = neville([1 / D(k) for k in ks], [s[k] for k in ks])
    ks = [100 - 6 * i for i in range(14)]
    for shift, nm in ((D(0), 'k'), (sig[f], 'k+σ')):
        c = fit([1 / (D(k) + shift) for k in ks], [F[k][0] * (D(k) + shift)**6 for k in ks])
        print('%s: a(%s)^6 = %s' % (f, nm, ' '.join('%+.6e' % x for x in c[:6])))
    g = lambda k: sum(cj / (D(k) + sig[f])**j for j, cj in enumerate(c)) / (D(k) + sig[f])**6
    head = sum(F[k][0] for k in F); kmax = max(F)
    tail = sum(g(k) for k in range(kmax + 1, 20001)) + c[0] / (5 * (D(20000.5) + sig[f])**5)
    print('%s: σ = %.8f  C = %.12e  Σ a(k) = %.15e (k ≤ %d: %.15e, tail %.3e)' %
          (f, sig[f], c[0], head + tail, kmax, head, tail))
if len(sys.argv) > 2:
    R = {}
    for l in open(sys.argv[2]):
        r = l.split()
        if r[3] != 'failed': R[(r[0][0], int(r[0][1:]), int(r[1]), int(r[2]))] = area(r) / par[r[0]]
    for f in sorted(fam):
        for pq in sorted({k[2:] for k in R if k[0] == f}, key=lambda x: (x[1], x[0])):
            ks = sorted(k for k in range(10, 200) if (f, k) + pq in R)
            est = []
            for J in (6, 9, 12):
                sub = [ks[int(i * (len(ks) - 1) / (J - 1))] for i in range(J)]
                est.append(neville([1 / (D(k) + sig[f]) for k in sub], [R[(f, k) + pq] for k in sub]))
            print('%s %d/%d (k ≤ %d): area/parent → %s' % (f, pq[0], pq[1], ks[-1], '  '.join('%.13f' % e for e in est)))
