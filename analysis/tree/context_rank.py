"""Rank of the context map, one level up: the child factor L(u·s; r2) = A(u·s, r2) q_{us}^4 q_2^4 / (3π/8)... over rows u
(prefixes of a cardioid-level parent word, q_u ≤ 12) and columns (s, r2) (s q ≤ 12, r2 q ≤ 8), against the parent's
own normalized area F_card(u·s) over rows u and columns s."""
import os, sys
import numpy as np
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from store import Store, val
from gen import words, rationals
st = Store()
U = [()] + [u for u in words(12, False) if val(list(u)) != 1]
S12 = words(12, True); R2 = rationals(8)
v = lambda u, s: val(list(u) + list(s))
st.request([(v(u, s), r2) for u in U for s in S12 for r2 in R2])
n = sum(len(x) for x in st.wanted.values()); print('computing', n, 'bulbs under', len(st.wanted), 'parents'); sys.stdout.flush()
st.compute()
nz = lambda x: np.nan if x is None else x
q4 = lambda x: x.denominator ** 4.0
Fc = np.array([[nz(st.F((v(u, s),))) for s in S12] for u in U])
L = np.array([[nz(st.A((v(u, s), r2))) * q4(r2) / nz(st.A((v(u, s),))) for s in S12 for r2 in R2] for u in U])
ok = ~(np.isnan(Fc).any(1) | np.isnan(L).any(1))
print('rows %d (complete %d), columns %d and %d' % (len(U), ok.sum(), Fc.shape[1], L.shape[1]))
def show(name, M):
    sv = np.linalg.svd(M[ok], compute_uv=False); sv = sv / sv[0]
    r = lambda t: int((sv > t).sum())
    print('%-34s rank at 1e-4/1e-6/1e-8/1e-10: %3d %3d %3d %3d  σ every 5th: %s' % (name, r(1e-4), r(1e-6), r(1e-8), r(1e-10),
          ' '.join('%.0e' % x for x in sv[:80:5])))
show('parent F_card(u·s)', Fc)
show('child factor λ(u·s; r2) q2^4', L)
for r2 in (Fraction(1, 2), Fraction(1, 3), Fraction(1, 7)):
    j = R2.index(r2)
    show('  one child r2 = %s' % r2, L[:, j::len(R2)])
