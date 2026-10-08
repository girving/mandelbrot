"""Is the area of the component landing at angle w̄ a rational series in the binary word w?

g(w) = area of the hyperbolic component whose root receives the ray of angle w̄ = (word w repeated), for every
nonzero word w of length ≤ 16 (non-primitive words land where their primitive root lands; 0^p, 1^p at the cusp).
Hankel H[u, s] = g̃(u·s) with g̃(w) = g(w) β^{|w|} (a geometric weight leaves rationality unchanged), rows |u| ≤ a,
columns |s| ≤ b."""
import sys, math, collections
import numpy as np
from load import load, A_CARD
roots = load(sys.argv[1])
area = {}
for r in roots:
    area[(r.p, r.lo)] = r.area; area[(r.p, r.hi)] = r.area
def prim(bits, n):
    """minimal period of the word (bits, n) and its value there"""
    for d in range(1, n + 1):
        if n % d: continue
        b = bits >> (n - d)
        if all(((bits >> (n - d * (j + 1))) & ((1 << d) - 1)) == b for j in range(n // d)): return d, b
def g(bits, n):
    d, b = prim(bits, n)
    if b == 0 or b == (1 << d) - 1: return A_CARD
    return area.get((d, b))
a, b = int(sys.argv[2]), int(sys.argv[3])
words = lambda m: [(bits, n) for n in range(0, m + 1) for bits in range(1 << n)]
U, S = words(a), [w for w in words(b) if w[1] > 0]
miss = 0
for beta in (1.0, 2.0, 4.0):
    H = np.zeros((len(U), len(S)))
    for i, (ub, un) in enumerate(U):
        for j, (sb, sn) in enumerate(S):
            v = g((ub << sn) | sb, un + sn)
            if v is None: miss += 1; v = 0
            H[i, j] = v * beta ** (un + sn)
    sv = np.linalg.svd(H, compute_uv=False); sv = sv / sv[0]
    r = lambda t: int((sv > t).sum())
    print('β=%.0f: %d x %d, missing %d; rank at 1e-2/1e-4/1e-6/1e-8: %d %d %d %d; σ: %s' % (beta, len(U), len(S), miss, r(1e-2), r(1e-4), r(1e-6),
          r(1e-8), ' '.join('%.0e' % x for x in sv[:60:4])))
