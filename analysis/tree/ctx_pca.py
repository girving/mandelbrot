"""How many directions do contexts span?  Rows F(r1; s) over digit suffixes s for all depth-1 contexts r1 (q1 ≤ 20),
minus the cardioid's F(s): singular values of the context deviations."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from store import Store, val
from gen import S40, rationals
st = Store()
ctx = rationals(20)
M = np.array([[st.F((r1, val(list(s)))) for s in S40] for r1 in ctx], dtype=float)
card = np.array([st.F((val(list(s)),)) for s in S40], dtype=float)
ok = ~np.isnan(M).any(0) & ~np.isnan(card)
M, card = M[:, ok], card[ok]
okr = ~np.isnan(M).any(1); M = M[okr]
print('%d contexts x %d suffixes' % M.shape)
for name, X in (('F(r1; s) - F_card(s)', M - card), ('F(r1; s) - mean over contexts', M - M.mean(0))):
    s = np.linalg.svd(X, compute_uv=False)
    print('%-30s max |X| %.1e; singular values: %s' % (name, np.abs(X).max(), ' '.join('%.0e' % x for x in s[:30])))
# Weighted by importance in the sum: suffix s weighted by q_s^-2 (its share of the bulb sum)
w = np.array([val(list(s)).denominator ** -2.0 for s, k in zip(S40, ok) if k])
s = np.linalg.svd((M - M.mean(0)) * w, compute_uv=False)
print('%-30s singular values: %s' % ('weighted by q_s^-2', ' '.join('%.0e' % x for x in s[:30])))
