"""Fits for the seahorse-valley families (seahorse.py): parents, and their satellite children.

  python3 seahorse_fit.py parents [kids]      (bulb_batch outputs; kids: "key P c_re c_im p q" jobs per parent)

Parents: phase σ_w = lim iπ/δ − 2k (δ = c + 3/4), C_w = lim a(k) k^4, the expansion of a(k)|2k+σ|^4/16 in 1/k
(fit through k = 30..120), and Σ_k a(k) with the tail from that expansion.  Children: limits of child/parent area
ratios by Neville in 1/k."""
import sys, math
import numpy as np

def neville(xs, ys):
    p = list(ys)
    for m in range(1, len(xs)):
        p = [(xs[i+m] * p[i] - xs[i] * p[i+1]) / (xs[i+m] - xs[i]) for i in range(len(p) - 1)]
    return p[0]

def lim(F, g, kmax, J=9, step=10):
    ks = [kmax - step * i for i in range(J)]
    return neville([1 / k for k in ks], [g(k) for k in ks])

fam, par = {}, {}
for l in open(sys.argv[1]):
    r = l.split(); n, k = r[0].split('_')
    a = float(r[5]) + float(r[10]); par[r[0]] = a
    fam.setdefault(n, {})[int(k)] = (complex(float(r[3]), float(r[4])) + 0.75, a)
for n, F in sorted(fam.items()):
    kmax = max(F)
    sig = lim(F, lambda k: 1j * math.pi / F[k][0] - 2 * k, kmax)
    C = lim(F, lambda k: F[k][1] * k**4, kmax)
    ks = np.arange(30, kmax + 1)
    h = np.array([F[k][1] * abs(2 * k + sig)**4 / 16 for k in ks])
    c = np.polynomial.polynomial.polyfit(1 / ks, h, 10)
    model = lambda k: np.polynomial.polynomial.polyval(1 / k, c) * 16 / abs(2 * k + sig)**4
    kk = np.arange(kmax + 1, 2000001, dtype=float)
    tail = model(kk).sum() + c[0] / (3 * (kk[-1] + 0.5)**3)
    head = sum(F[k][1] for k in F)
    print('%s: σ = %.9f%+.9fi  C = %.11e  a|2k+σ|^4/16 = %s' % (n, sig.real, sig.imag, C,
          ' '.join('%+.4e' % x for x in c[:4])))
    print('    Σ a(k), k ≥ %d: %.12e (k ≤ %d: %.12e, tail %.6e)' % (min(F), head + tail, kmax, head, tail))
if len(sys.argv) > 2:
    R = {}
    for l in open(sys.argv[2]):
        r = l.split()
        if r[3] == 'failed': continue
        n, k = r[0].split('_')
        R[(n, int(k), int(r[1]), int(r[2]))] = (float(r[5]) + float(r[10])) / par[r[0]]
    for n in sorted(fam):
        for pq in sorted({x[2:] for x in R if x[0] == n}, key=lambda x: (x[1], x[0])):
            kmax = max(x[1] for x in R if x[0] == n and x[2:] == pq)
            est = [lim(None, lambda k: R[(n, k) + pq], kmax, J) for J in (5, 7, 9)]
            print('%s %d/%d: area/parent → %s   (k = 20: %.10f)' % (n, pq[0], pq[1], '  '.join('%.12f' % e for e in est),
                  R[(n, 20) + pq]))
