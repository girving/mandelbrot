import os, sys, math
import numpy as np
from fractions import Fraction
from math import gcd
exec(open('wfa_hp_sum.py').read().split("\nif __name__ == '__main__'")[0])
K = int(sys.argv[1]) if len(sys.argv) > 1 else 40
al, Nb, Bt = build(K)
groups = {}
for q in range(401, 501):
    for p in range(1, q):
        if gcd(p, q) != 1: continue
        x = Fraction(p, q); w = cf(x); v = al.copy()
        for a in w[:-1]: v = v @ Nb(a)
        W = np.pi * math.sin(math.pi * p / q) ** 2 / q ** 4
        err = W * ((v @ Bt(w[-1])) - F[min(x, 1 - x)])
        mi = max(w[:-1]) if len(w) > 1 else 0
        key = ('interior ≤ 3' if mi <= 3 else 'interior 4..64' if mi <= 64 else 'interior > 64') + (', last > 60' if w[-1] > 60 else '')
        g = groups.setdefault(key, [0, 0.0, 0.0]); g[0] += 1; g[1] += err; g[2] = max(g[2], abs(err / W))
for k, (n, e, m) in sorted(groups.items()):
    print('%-28s n=%5d  Σ W·(F_model - F) %+.2e   max |F_model - F| %.1e' % (k, n, e, m))
