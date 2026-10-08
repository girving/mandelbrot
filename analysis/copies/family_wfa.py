"""Seahorse family constants as a weighted automaton over kneading tails.

  limb_families keys J > keys  (K0 = the census limb);  python3 family_wfa.py keys census.txt[.gz] K0 [a b]

Every family's kneading sequence is 1^{2k} 0 t with a free binary tail t of length j - 1 (all 2^{j-1} words occur, for
|t| < 2k), so G(t) = C_w is a function on the free binary monoid.  C_w = lim a k^4 by Neville in 1/k over the census k.
Spectral learning: Hankel block H[u, s] = G(us) over |u| ≤ a, |s| ≤ b, rank-r SVD H ≈ P Q, A_x = P⁺ H_x Q⁺ with
H_x[u, s] = G(u x s), α = Q⁺-row of the empty prefix, β = P⁺ G(U).  Tested on all tails longer than a + b + 1 (never
in the block): relative errors, and the error of their total, which is what the area needs."""
import sys, gzip, itertools, math
from collections import defaultdict
from fractions import Fraction as Fr
import numpy as np

def kneading(w):
    """ν_i = [2^{i-1}θ ∈ (θ/2, (θ+1)/2)], i = 1..p-1, for θ = .(w)^∞"""
    p = len(w); t = Fr(int(w, 2), 2**p - 1); lo, hi = t / 2, (t + 1) / 2
    out, x = '', t
    for _ in range(1, p):
        out += '1' if lo < x < hi else '0'
        x = (2 * x) % 1
    return out

def neville(xs, ys):
    p = list(ys)
    for m in range(1, len(xs)):
        p = [(xs[i+m] * p[i] - xs[i] * p[i+1]) / (xs[i+m] - xs[i]) for i in range(len(p) - 1)]
    return p[0]

def constants(keys_path, census_path, K0):
    F = defaultdict(dict)
    for l in (gzip.open(census_path, 'rt') if census_path.endswith('.gz') else open(census_path)):
        r = l.split()
        if r[3] == 'failed': continue
        j, i, k = r[0].split('_')
        F[(int(j[1:]), int(i))][int(k)] = float(r[5]) + float(r[10])
    head, pre = '1' * (2 * K0) + '0', '01' * (K0 - 1)
    G, err = {}, {}
    for l in open(keys_path):
        j, i, lo, hi = l.split(); f = F[(int(j), int(i))]; ks = sorted(f)
        kn = kneading(pre + lo); assert kn.startswith(head), kn
        C = neville([1 / k for k in ks], [f[k] * k**4 for k in ks])
        C2 = neville([1 / k for k in ks[1:]], [f[k] * k**4 for k in ks[1:]])
        G[kn[len(head):]] = C; err[kn[len(head):]] = abs(C2 / C - 1)
    return G, err

words = lambda n: [''.join(w) for m in range(n + 1) for w in itertools.product('01', repeat=m)]

if __name__ == '__main__':
    K0 = int(sys.argv[3]); a, b = (int(sys.argv[4]), int(sys.argv[5])) if len(sys.argv) > 5 else (5, 6)
    G, err = constants(sys.argv[1], sys.argv[2], K0)
    L = 0  # the longest length with every tail present (|t| = 2k is the bulb's period doubling, not a family)
    while all(w in G for w in itertools.product('01', repeat=L + 1) for w in [''.join(w)]): L += 1
    G = {w: v for w, v in G.items() if len(w) <= L}
    e = sorted(err.values())
    print('%d tails of length ≤ %d; C_w extrapolation error median %.1e, 99%% %.1e, max %.1e' %
          (len(G), L, e[len(e) // 2], e[int(.99 * len(e))], e[-1]))
    U, V = words(a), words(b)
    H = np.array([[G[u + s] for s in V] for u in U])
    Hx = {x: np.array([[G[u + x + s] for s in V] for u in U]) for x in '01'}
    hU = np.array([G[u] for u in U])
    Uu, S, Vt = np.linalg.svd(H, full_matrices=False)
    print('Hankel %dx%d: σ/σ0 = %s' % (*H.shape, ' '.join('%.0e' % x for x in S[:40:2] / S[0])))
    tests = [w for w in G if len(w) > a + b + 1]
    tot = sum(G[w] for w in tests)
    print('tests: %d tails of length %d..%d, total %.6e' % (len(tests), a + b + 2, L, tot))
    for r in (4, 8, 12, 16, 20, 24, 28, 32, 40, 48):
        if r > len(S): break
        P = Uu[:, :r] * S[:r]; Q = Vt[:r]
        Pp, Qp = np.linalg.pinv(P), np.linalg.pinv(Q)
        A = {x: Pp @ Hx[x] @ Qp for x in '01'}
        al = P[0]; be = Pp @ hU   # row 0 of U is the empty prefix
        def pred(w):
            v = al
            for x in w: v = v @ A[x]
            return v @ be
        rel = [abs(pred(w) / G[w] - 1) for w in tests]
        ptot = sum(pred(w) for w in tests)
        rho = max(abs(np.linalg.eigvals(A['0'] + A['1'])))
        print('  r=%2d: relative error median %.1e, 90%% %.1e; total error %.1e; spectral radius of A0+A1 %.3f' %
              (r, np.median(rel), np.quantile(rel, .9), abs(ptot / tot - 1), rho))
