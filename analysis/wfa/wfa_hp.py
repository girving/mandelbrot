"""High-precision automaton: double-double bulbs (bulb_areas BULB_EXP=1, F = col6 + col11), suffix-closure learning,
held-out error against K"""
import os, sys, pickle
import numpy as np
from fractions import Fraction
D = os.environ['BULB_DATA'] + '/'
def load(*names):
    F = {}
    for n in names:
        for line in open(D + n):
            t = line.split()
            if t[2] == 'failed': continue
            p, q = int(t[0]), int(t[1]); f = float(t[5]) + (float(t[10]) if len(t) > 10 else 0.0)
            F[Fraction(p, q)] = F[Fraction(q - p, q)] = f
    return F
def val(w):
    x = Fraction(0)
    for a in reversed(w): x = 1 / (a + x)
    return x
def cf(x):
    a = []; p, q = x.numerator, x.denominator
    while q: a.append(p // q); p, q = q, p % q
    return a[1:]
if __name__ == '__main__':
    F = load(sys.argv[1] if len(sys.argv) > 1 else 'all_e.out')
    U, S = pickle.load(open(D + 'hankel20_sets.pkl', 'rb'))
    H = np.array([[F[val(list(u) + list(s))] for s in S] for u in U])
    Uidx = {u: i for i, u in enumerate(U)}; Sidx = {s: i for i, s in enumerate(S)}
    Uu, Ss, Vt = np.linalg.svd(H, full_matrices=False)
    print('Hankel %dx%d singular values:' % H.shape, ' '.join('%.1e' % x for x in Ss[:70:2]), '(every other)')
    vw = []
    for x in F:
        w = cf(x)
        if len(w) >= 4 and all(a <= 4 for a in w[:-1]) and w[-1] <= 60 and x.denominator > 2000: vw.append((w, F[x]))
    print(len(vw), 'held-out long words (q > 2000, interior digits ≤ 4)')
    for K in list(range(10, 41, 5)) + list(range(44, 81, 4)):
        P = Uu[:, :K] * Ss[:K]; Q = Vt[:K]; Pp = np.linalg.pinv(P); N = {}
        for b in range(1, 5):
            cols = [(Sidx[s], Sidx[(b,) + s]) for s in S if (b,) + s in Sidx]
            N[b] = Pp @ H[:, [c[1] for c in cols]] @ np.linalg.pinv(Q[:, [c[0] for c in cols]])
        al = P[Uidx[()]]; be = {s[0]: Q[:, Sidx[s]] for s in S if len(s) == 1}
        e = []
        for w, f in vw:
            v = al.copy()
            for a in w[:-1]: v = v @ N[a]
            e.append(abs(v @ be[w[-1]] - f))
        print('K=%2d: median %.1e  90%% %.1e  max %.1e' % (K, np.median(e), np.quantile(e, 0.9), max(e)))
