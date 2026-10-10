"""Separator learning from many contexts.  Digits as in learn3 (root + context rows, suffix closure); each context's
start vector α(ctx) = F(ctx; ·) Q⁺ from its row over the digit suffixes.  Separator models, fitted on 70% of the
depth-1 contexts (q1 ≤ 20) and tested on the rest and on depth-2 contexts:
  linear:   α(r1) = s(r1)ᵀ Σ,                      s(r1) = α_rootᵀ N(digits of r1)
  analytic: α(r1) = Σ_k T_k(2 x_{r1} - 1) s(r1)ᵀ Σ_k
  BULB_DATA=dir python3 learn4.py K [K ...]"""
import os, sys
import numpy as np
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
exec(open('learn3.py').read().split("T = os.environ")[0].replace("rng = np.random.default_rng(7)\ntest = rng.random(len(cols)) < 0.2", "test = np.zeros(len(cols), bool)"))
from gen import rationals
ctx1 = [(r1,) for r1 in rationals(20)]
ctx2 = [(r1, r2) for r1 in (Fraction(1, 2), Fraction(1, 3), Fraction(2, 5), Fraction(1, 4), Fraction(2, 7)) for r2 in rationals(6)]
def ctx_row(c):
    return np.array([st.F(c + (val(list(s)),)) for s in ctr], dtype=float)
R1rows = {c: ctx_row(c) for c in ctx1}; R2rows = {c: ctx_row(c) for c in ctx2}
R1rows = {c: v for c, v in R1rows.items() if not np.isnan(v).any()}; R2rows = {c: v for c, v in R2rows.items() if not np.isnan(v).any()}
print('%d depth-1 contexts, %d depth-2 contexts with complete rows' % (len(R1rows), len(R2rows)))
rng = np.random.default_rng(11)
train = {c: rng.random() < 0.7 for c in R1rows}
for K in map(int, sys.argv[1:]):
    P = Uu[:, :K] * Ss[:K]; Q = Vt[:K]; Pp = np.linalg.pinv(P); Qp = np.linalg.pinv(Q)
    N = {}
    for b in range(1, 21):
        pr = [(cidx[s], cidx[(b,) + s]) for s in ctr if (b,) + s in cidx]
        if len(pr) >= 3: N[b] = Pp @ Htr[:, [p[1] for p in pr]] @ np.linalg.pinv(Q[:, [p[0] for p in pr]])
    al_root = P[ridx[('root', ())]]
    alpha = {c: v @ Qp for c, v in list(R1rows.items()) + list(R2rows.items())}
    recon = np.median([np.abs(alpha[c] @ Q - R1rows[c]).max() for c in R1rows])
    def state(v, w):
        for a in w:
            if a not in N: return None
            v = v @ N[a]
        return v
    X, Y, xs = [], [], []
    for c in R1rows:
        s = state(al_root, cf(c[0]))
        if s is None: continue
        X.append(s); Y.append(alpha[c]); xs.append(float(c[0]))
    X, Y, xs = np.array(X), np.array(Y), np.array(xs)
    tr = np.array([train[c] for c in R1rows if state(al_root, cf(c[0])) is not None])
    out = ['row reconstruction %.1e' % recon]
    for name, deg in (('linear', 0), ('analytic deg 2', 2), ('analytic deg 4', 4)):
        Z = np.hstack([X * np.polynomial.chebyshev.chebval(2 * xs - 1, np.eye(deg + 1)[k])[:, None] for k in range(deg + 1)])
        W = np.linalg.lstsq(Z[tr], Y[tr], rcond=None)[0]
        def sep(s, x): return np.hstack([s * np.polynomial.chebyshev.chebval(2 * x - 1, np.eye(deg + 1)[k]) for k in range(deg + 1)]) @ W
        e_tr = np.median(np.abs((Z[tr] @ W) @ Q - np.array([R1rows[c] for c, t in zip([c for c in R1rows if state(al_root, cf(c[0])) is not None], tr) if t])).max(1))
        e_te = np.median(np.abs((Z[~tr] @ W) @ Q - np.array([R1rows[c] for c, t in zip([c for c in R1rows if state(al_root, cf(c[0])) is not None], tr) if not t])).max(1))
        e2 = []
        for c, row in R2rows.items():
            a1 = alpha.get((c[0],))
            if a1 is None: continue
            s2 = state(a1, cf(c[1]))
            if s2 is None: continue
            e2.append(np.abs(sep(s2, float(c[1])) @ Q - row).max())
        out.append('%s: train %.1e test %.1e depth-2 contexts %.1e' % (name, e_tr, e_te, np.median(e2)))
    print('K=%d (%d train, %d test contexts): %s' % (K, tr.sum(), (~tr).sum(), ' | '.join(out)))
    sys.stdout.flush()
