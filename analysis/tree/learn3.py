"""Two-step automaton for the satellite tree: digits from a block with root and context rows, then the separator.

  BULB_DATA=dir python3 learn3.py K [K ...]

Rows: root prefixes u (q_u ≤ 20) and context prefixes (r1, u) (r1 ≤ 1/2, q1 ≤ 8; q_u ≤ 8); columns: digit suffixes
s (q_s ≤ 40), F(row·s) = F_card(u·s) or F(r1; u·s).  N(b) by suffix closure; start vectors are rows of P: α_root =
P[()], α(r1) = P[(r1, ())].  Separator: Σ = argmin Σ_r1 |α_rootᵀ N(digits of r1) Σ - α(r1)ᵀ|² over the training
contexts.  Tests: held-out suffix columns (20%), depth-2 data with unseen contexts (tree/d2, q1 ≤ 16), depth 3."""
import os, sys, glob
import numpy as np
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from store import Store, val, read as read_out
from gen import U20, U8, S40, R1, rationals

def cf(x):
    a = []; p, q = x.numerator, x.denominator
    while q: a.append(p // q); p, q = q, p % q
    return tuple(a[1:])

st = Store()
rows = [('root', u) for u in U20 if val(list(u)) != 1] + [('ctx', r1, u) for r1 in R1 for u in U8 if val(list(u)) != 1]
def F(row, s):
    if row[0] == 'root': return st.F((val(list(row[1]) + list(s)),))
    return st.F((row[1], val(list(row[2]) + list(s))))
H = np.array([[F(r, s) for s in S40] for r in rows], dtype=float)
okc = ~np.isnan(H).any(0); okr = ~np.isnan(H[:, okc]).any(1)
cols = [s for s, g in zip(S40, okc) if g]; rows = [r for r, g in zip(rows, okr) if g]; H = H[okr][:, okc]
print('block %d x %d' % H.shape)
rng = np.random.default_rng(7)
test = rng.random(len(cols)) < 0.2
ridx = {r: i for i, r in enumerate(rows)}
Htr = H[:, ~test]; ctr = [c for c, t in zip(cols, test) if not t]; cidx = {c: j for j, c in enumerate(ctr)}
Uu, Ss, Vt = np.linalg.svd(Htr, full_matrices=False)
print('singular values:', ' '.join('%.1e' % x for x in Ss[:70:5]), '(every fifth)')
T = os.environ['BULB_DATA'] + '/tree/'
d2 = []
for f in glob.glob(T + 'd2/*.out'):
    r1 = Fraction(*map(int, os.path.basename(f)[:-4].split('-')))
    for x, v in read_out(f).items(): d2.append((r1, x, v[0]))
d3 = []
for f in glob.glob(T + 'd3/*.out'):
    a, b = os.path.basename(f)[:-4].split('.')
    r1 = Fraction(*map(int, a.split('-'))); r2 = Fraction(*map(int, b.split('-')))
    for x, v in read_out(f).items(): d3.append((r1, r2, x, v[0]))
for K in map(int, sys.argv[1:]):
    P = Uu[:, :K] * Ss[:K]; Q = Vt[:K]; Pp = np.linalg.pinv(P)
    N = {}
    for b in range(1, 13):
        pr = [(cidx[s], cidx[(b,) + s]) for s in ctr if (b,) + s in cidx]
        if len(pr) >= 3: N[b] = Pp @ Htr[:, [p[1] for p in pr]] @ np.linalg.pinv(Q[:, [p[0] for p in pr]])
    # β for any suffix: Q column for training suffixes; else from the closure: β(s) = N(s1) ... β(last)
    be = {c[0]: Q[:, cidx[c]] for c in ctr if len(c) == 1}
    def col(s):
        if s in cidx: return Q[:, cidx[s]]
        if s[-1] not in be or any(a not in N for a in s[:-1]): return None
        v = be[s[-1]]
        for a in reversed(s[:-1]): v = N[a] @ v
        return v
    # Held-out columns, each row its own start state
    e = []
    for j, (s, t) in enumerate(zip(cols, test)):
        if not t: continue
        v = col(s)
        if v is not None: e.append(np.abs(P @ v - H[:, j]).max())
    # Separator from training contexts
    al_root = P[ridx[('root', ())]]
    X, Y = [], []
    for r1 in R1:
        if ('ctx', r1, ()) not in ridx or any(a not in N for a in cf(r1)): continue
        v = al_root.copy()
        for a in cf(r1): v = v @ N[a]
        X.append(v); Y.append(P[ridx[('ctx', r1, ())]])
    Sig = np.linalg.lstsq(np.array(X), np.array(Y), rcond=None)[0]
    def pred(addr):
        v = al_root.copy()
        for li, r in enumerate(addr):
            w = cf(r)
            if li < len(addr) - 1:
                if any(a not in N for a in w): return None
                for a in w: v = v @ N[a]
                v = v @ Sig
            else:
                c = col(w)
                return None if c is None else v @ c
    e2 = [abs(pred((r1, x)) - v) for r1, x, v in d2 if r1 not in R1 and pred((r1, x)) is not None]
    e2s = [abs(pred((r1, x)) - v) for r1, x, v in d2 if r1 in R1 and pred((r1, x)) is not None]
    e3 = [abs(pred((r1, r2, x)) - v) for r1, r2, x, v in d3 if pred((r1, r2, x)) is not None]
    q = lambda e: 'n=%d median %.1e 90%% %.1e' % (len(e), np.median(e), np.quantile(e, .9)) if e else 'n=0'
    print('K=%d: held-out columns (max over rows) %s | depth 2 seen contexts %s | unseen contexts %s | depth 3 %s'
          % (K, q(e), q(e2s), q(e2), q(e3)))
    sys.stdout.flush()
