"""Joint automaton for the satellite tree from a context-rich block, by row closure.

  BULB_DATA=dir python3 learn5.py K [K ...]

Rows: root prefixes u (q_u ≤ 20) and context prefixes (r1, u) (all r1 with q1 ≤ 20, q_u ≤ 4, u may be empty), minus
the context rows of a held-out 15% of the r1.  Columns: digit suffixes s (q_s ≤ 40) and separator futures #s'
(q_s' ≤ 20); invalid words (an empty level) are 0.  With H ≈ P Q (rank K):
  N(b) = argmin |P[rows u] N - P[rows u·b]|  over root and context rows with u·b also a row,
  N(#) = argmin |P[root r1] N - P[(r1, ())]| over the training contexts,
  α = P[root ()], β = P⁺ H[:, ε] (each row's own value F(prefix)).
F(w) = αᵀ N(w_1) ⋯ N(w_n) β.  Tests: children of held-out contexts (tree/d2 and the block), depth-3 nodes (r1, r2)
with q2 ≥ 5 (not in the block), and depth-3 data tree/d3."""
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
cols = [('dig', s) for s in S40] + [('sep', s2) for s2 in rationals(20)] + [('eps',)]

def addr(row, col):
    if row[0] == 'root':
        u = row[1]
        if col[0] == 'dig': return (val(list(u) + list(col[1])),)
        if col[0] == 'sep': return (val(list(u)), col[1]) if u else None
        return (val(list(u)),) if u else None
    r1, u = row[1], row[2]
    if col[0] == 'dig': return (r1, val(list(u) + list(col[1])))
    if col[0] == 'sep': return (r1, val(list(u)), col[1]) if u else None
    return (r1, val(list(u))) if u else None

def value(a):
    return 0.0 if a is None else st.F(a)

H = np.array([[value(addr(r, c)) for c in cols] for r in rows], dtype=float)
bad = np.isnan(H)
print('block %d x %d, missing %d' % (*H.shape, bad.sum()))
keep_r = ~bad.any(1); rows = [r for r, k in zip(rows, keep_r) if k]; H = H[keep_r]
ridx = {r: i for i, r in enumerate(rows)}
ie = len(cols) - 1                       # the empty-suffix column
Uu, Ss, Vt = np.linalg.svd(H[:, :ie], full_matrices=False)
print('using %d rows; singular values (every fifth): %s' % (len(rows), ' '.join('%.1e' % x for x in Ss[:100:5])))

# Tests
T = os.environ['BULB_DATA'] + '/tree/'
test_held = [((r1, x), st.F((r1, x))) for r1 in held for x in rationals(30)]
test_held = [(a, v) for a, v in test_held if v is not None]
test_d3 = []
for r1 in (Fraction(1, 2), Fraction(1, 3), Fraction(2, 5), Fraction(1, 4), Fraction(2, 7)):
    for r2 in rationals(6):
        if r2.denominator < 5: continue
        for s in S40:
            a = (r1, r2, val(list(s))); v = st.F(a)
            if v is not None: test_d3.append((a, v))
for f in glob.glob(T + 'd3/*.out'):
    a_, b_ = os.path.basename(f)[:-4].split('.')
    r1 = Fraction(*map(int, a_.split('-'))); r2 = Fraction(*map(int, b_.split('-')))
    if r2.denominator < 5: continue
    for x, v in read_out(f).items(): test_d3.append(((r1, r2, x), v[0]))
print('tests: %d children of %d held-out contexts, %d depth-3 nodes' % (len(test_held), len(held), len(test_d3)))

for K in map(int, sys.argv[1:]):
    P = Uu[:, :K] * Ss[:K]
    N = {}
    for b in range(1, 21):
        pr = [(ridx[r], ridx[('root', r[1] + (b,))]) for r in rows if r[0] == 'root' and ('root', r[1] + (b,)) in ridx] + \
             [(ridx[r], ridx[('ctx', r[1], r[2] + (b,))]) for r in rows if r[0] == 'ctx' and ('ctx', r[1], r[2] + (b,)) in ridx]
        if len(pr) >= K: N[b] = np.linalg.lstsq(P[[p[0] for p in pr]], P[[p[1] for p in pr]], rcond=None)[0]
    pr = [(ridx[('root', cf(r1))], ridx[('ctx', r1, ())]) for r1 in ctxs
          if ('root', cf(r1)) in ridx and ('ctx', r1, ()) in ridx]
    N['#'] = np.linalg.lstsq(P[[p[0] for p in pr]], P[[p[1] for p in pr]], rcond=None)[0]
    al = P[ridx[('root', ())]]; be = np.linalg.pinv(P) @ H[:, ie]
    def Fm(a):
        w = []
        for i, r in enumerate(a):
            if i: w.append('#')
            w += list(cf(r))
        if any(s not in N for s in w): return None
        v = al
        for s in w: v = v @ N[s]
        return v @ be
    out = []
    for nm, data in (('held-out contexts', test_held), ('depth 3', test_d3)):
        e = [abs(Fm(a) - v) for a, v in data if Fm(a) is not None]
        out.append('%s n=%d median %.1e 90%% %.1e max %.1e' % (nm, len(e), np.median(e), np.quantile(e, .9), max(e)) if e else nm + ' n=0')
    sep_res = np.abs(P[[p[0] for p in pr]] @ N['#'] - P[[p[1] for p in pr]]).max()
    print('K=%d: separator closure residual %.1e, digits %s | %s' % (K, sep_res, sorted(k for k in N if k != '#')[-1], ' | '.join(out)))
    sys.stdout.flush()
