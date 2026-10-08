"""Is the multiplicative cascade low-rank?  Rows (r1, u) (r1 ≤ 1/2 with q1 ≤ 8, u with q_u ≤ 4), compared over:
  F(r1; u·s)            within-level normalized area (baseline), columns s (q_s ≤ 12)
  G2(r1, u·s)           = A(r1, u·s) q1^4 q2^4, the level-2 cumulative area without q's
  G3(r1, u·s, r3)       = A(r1, u·s, r3) q1^4 q2^4 q3^4, columns (s, r3), r3 with q3 ≤ 8
  G3 / G2               the level-3 factor alone (its context dependence)
If rank(G3) ≈ rank(G2) the cascade can be summed by one automaton state; if it multiplies, products need tensors."""
import os, sys
import numpy as np
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from store import Store, val
from gen import words, rationals
st = Store()
R1 = [x for x in rationals(8) if x <= Fraction(1, 2)]
U4 = [()] + [u for u in words(4, False) if val(list(u)) != 1]
S12 = words(12, True)
R3 = rationals(8)
rows = [(r1, u) for r1 in R1 for u in U4]
v2 = lambda r, s: val(list(r[1]) + list(s))
st.request([(r[0], v2(r, s)) for r in rows for s in S12] + [(r[0], v2(r, s), r3) for r in rows for s in S12 for r3 in R3])
n = sum(len(v) for v in st.wanted.values()); print('computing', n, 'bulbs under', len(st.wanted), 'parents'); sys.stdout.flush()
st.compute()
q = lambda x: x.denominator ** 4.0
nz = lambda v: np.nan if v is None else v
F = np.array([[nz(st.F((r[0], v2(r, s)))) for s in S12] for r in rows], dtype=float)
G2 = np.array([[nz(st.A((r[0], v2(r, s)))) * q(r[0]) * q(v2(r, s)) for s in S12] for r in rows], dtype=float)
G3 = np.array([[nz(st.A((r[0], v2(r, s), r3))) * q(r[0]) * q(v2(r, s)) * q(r3) for s in S12 for r3 in R3] for r in rows], dtype=float)
G2rep = np.repeat(G2, len(R3), axis=1)
print('rows %d, columns %d (level 2) and %d (level 3); missing %d %d %d' % (len(rows), len(S12), G3.shape[1],
      np.isnan(F).sum(), np.isnan(G2).sum(), np.isnan(G3).sum()))
okr = ~(np.isnan(F).any(1) | np.isnan(G2).any(1) | np.isnan(G3).any(1))
def spec(M):
    s = np.linalg.svd(M[okr], compute_uv=False); return s / s[0]
def rank(sv, tol): return int((sv > tol).sum())
for name, M in (('F (within level)', F), ('G2 = A q1^4 q2^4', G2), ('G3 = A q1^4 q2^4 q3^4', G3), ('G3 / G2 (level-3 factor)', G3 / G2rep)):
    sv = spec(np.nan_to_num(M))
    print('%-26s rank at 1e-4/1e-6/1e-8/1e-10: %3d %3d %3d %3d   σ: %s' % (name, rank(sv, 1e-4), rank(sv, 1e-6), rank(sv, 1e-8),
          rank(sv, 1e-10), ' '.join('%.0e' % x for x in sv[:40:4])))
