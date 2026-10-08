"""Stable digits + many-context separator.  Block: rows root prefixes u (q_u ≤ 20) and context rows (r1, u) (training
r1 with q1 ≤ 20, q_u ≤ 4); columns digit suffixes only (q_s ≤ 40).  N(b) by suffix closure (stable, as learn3);
N(#) = argmin |P[root r1] N - P[(r1, ())]| over training contexts.  Tests as learn5.
  BULB_DATA=dir python3 learn6.py K [K ...]"""
import os, sys, glob
import numpy as np
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from store import Store, val, read as read_out
from gen import S40, U20, rationals, words

def cf(x):
    a = []; p, q = x.numerator, x.denominator
    while q: a.append(p // q); p, q = q, p % q
    return tuple(a[1:])

st = Store()
U4 = [()] + [u for u in words(4, False) if val(list(u)) != 1]
ctxs = rationals(20)
rng = np.random.default_rng(5)
held = {r1 for r1 in ctxs if rng.random() < 0.15}
rows = [('root', u) for u in U20 if not u or val(list(u)) != 1] + \
       [('ctx', r1, u) for r1 in ctxs if r1 not in held for u in U4]
def addr(row, s):
    if row[0] == 'root': return (val(list(row[1]) + list(s)),)
    return (row[1], val(list(row[2]) + list(s)))
H = np.array([[st.F(addr(r, s)) for s in S40] for r in rows], dtype=float)
okc = ~np.isnan(H).any(0); H = H[:, okc]; cols = [s for s, k in zip(S40, okc) if k]
okr = ~np.isnan(H).any(1); H = H[okr]; rows = [r for r, k in zip(rows, okr) if k]
ridx = {r: i for i, r in enumerate(rows)}; cidx = {s: j for j, s in enumerate(cols)}
Uu, Ss, Vt = np.linalg.svd(H, full_matrices=False)
print('block %d x %d; singular values (every fifth): %s' % (*H.shape, ' '.join('%.1e' % x for x in Ss[:100:5])))
T = os.environ['BULB_DATA'] + '/tree/'
test_held = [((r1, x), st.F((r1, x))) for r1 in held for x in rationals(30)]
test_held = [(a, v) for a, v in test_held if v is not None]
test_seen = [((r1, x), st.F((r1, x))) for r1 in list(set(ctxs) - held)[:20] for x in rationals(30) if x.denominator > 10]
test_seen = [(a, v) for a, v in test_seen if v is not None]
test_d3 = []
for r1 in (Fraction(1, 2), Fraction(1, 3), Fraction(2, 5), Fraction(1, 4), Fraction(2, 7)):
    for r2 in rationals(6):
        for s in S40:
            a = (r1, r2, val(list(s))); v = st.F(a)
            if v is not None: test_d3.append((a, v))
print('tests: %d seen-context children, %d held-out-context children, %d depth-3 nodes' % (len(test_seen), len(test_held), len(test_d3)))
for K in map(int, sys.argv[1:]):
    P = Uu[:, :K] * Ss[:K]; Q = Vt[:K]; Pp = np.linalg.pinv(P)
    N = {}
    for b in range(1, 13):
        pr = [(cidx[s], cidx[(b,) + s]) for s in cols if (b,) + s in cidx]
        if len(pr) >= 3: N[b] = Pp @ H[:, [p[1] for p in pr]] @ np.linalg.pinv(Q[:, [p[0] for p in pr]])
    pr = [(ridx[('root', cf(r1))], ridx[('ctx', r1, ())]) for r1 in ctxs if ('root', cf(r1)) in ridx and ('ctx', r1, ()) in ridx]
    N['#'] = np.linalg.lstsq(P[[p[0] for p in pr]], P[[p[1] for p in pr]], rcond=None)[0]
    sep_res = np.abs(P[[p[0] for p in pr]] @ N['#'] - P[[p[1] for p in pr]]).max()
    al = P[ridx[('root', ())]]
    be = {s[0]: Q[:, cidx[s]] for s in cols if len(s) == 1}
    def Fm(a):
        v = al
        for i, r in enumerate(a):
            w = cf(r)
            if i < len(a) - 1:
                if any(c not in N for c in w): return None
                for c in w: v = v @ N[c]
                v = v @ N['#']
            else:
                if w[-1] not in be or any(c not in N for c in w[:-1]): return None
                for c in w[:-1]: v = v @ N[c]
                return v @ be[w[-1]]
    out = []
    for nm, data in (('seen contexts', test_seen), ('held-out contexts', test_held), ('depth 3', test_d3)):
        e = [abs(Fm(a) - v) for a, v in data if Fm(a) is not None]
        out.append('%s n=%d median %.1e 90%% %.1e' % (nm, len(e), np.median(e), np.quantile(e, .9)) if e else nm + ' n=0')
    print('K=%d: separator residual %.1e | %s' % (K, sep_res, ' | '.join(out)))
    sys.stdout.flush()
