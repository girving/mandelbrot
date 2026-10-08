"""Does the cascade's rank saturate with depth?  Rows u (prefixes of a cardioid-level word, q_u ≤ 16); columns
(s) [F_card(u·s)], (s, r2) [G2 = A(u·s, r2) q_{us}^4 q_2^4], (s, r2, r3) [G3 = A(u·s, r2, r3) q_{us}^4 q_2^4 q_3^4],
with s q ≤ 6, r2 and r3 q ≤ 5.  If rank(G3) ≈ rank(G2), one automaton state carries the cascade at any depth."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from store import Store, val
from gen import words, rationals
st = Store()
U = [()] + [u for u in words(16, False) if val(list(u)) != 1]
S = words(6, True); R = rationals(5)
v = lambda u, s: val(list(u) + list(s))
st.request([(v(u, s), r2, r3) for u in U for s in S for r2 in R for r3 in R])
n = sum(len(x) for x in st.wanted.values()); print('computing', n, 'bulbs under', len(st.wanted), 'parents'); sys.stdout.flush()
st.compute()
nz = lambda x: np.nan if x is None else x
q4 = lambda x: x.denominator ** 4.0
F1 = np.array([[nz(st.F((v(u, s),))) for s in S] for u in U])
G2 = np.array([[nz(st.A((v(u, s), r2))) * q4(v(u, s)) * q4(r2) for s in S for r2 in R] for u in U])
G3 = np.array([[nz(st.A((v(u, s), r2, r3))) * q4(v(u, s)) * q4(r2) * q4(r3) for s in S for r2 in R for r3 in R] for u in U])
ok = ~(np.isnan(F1).any(1) | np.isnan(G2).any(1) | np.isnan(G3).any(1))
print('rows %d (complete %d); columns %d, %d, %d' % (len(U), ok.sum(), F1.shape[1], G2.shape[1], G3.shape[1]))
for name, M in (('F_card(u·s)', F1), ('G2 (u·s, r2)', G2), ('G3 (u·s, r2, r3)', G3)):
    sv = np.linalg.svd(M[ok], compute_uv=False); sv = sv / sv[0]
    r = lambda t: int((sv > t).sum())
    print('%-18s rank at 1e-4/1e-6/1e-8/1e-10: %3d %3d %3d %3d  σ every 5th: %s' % (name, r(1e-4), r(1e-6), r(1e-8), r(1e-10),
          ' '.join('%.0e' % x for x in sv[:120:5])))
