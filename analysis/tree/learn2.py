"""Joint automaton over digits and '#' by the basis method: F(w) = αᵀ N(w_1) ⋯ N(w_n) β.

  BULB_DATA=dir python3 learn2.py K [K ...]

Candidate prefixes: root digit words (q ≤ 14) and context words w#u (w complete, q_w ≤ 8; q_u ≤ 6, u may be
empty).  Candidate suffixes: complete digit words (q ≤ 20), #s' (q ≤ 12), s#s' (q_s ≤ 5, q_s' ≤ 8).  From their
Hankel block pick a pivoted basis of prefixes and suffixes, compute the shifted blocks F(p σ s) for σ ∈ {1..12, #}
and the empty prefix/suffix, and set N(σ) = P⁺ H_σ Q⁺, α = H[ε, S] Q⁺, β = P⁺ H[P, ε].  Invalid words (an empty
level) have F = 0.  Validation on depth-2 data (tree/d2) and depth-3 data (tree/d3)."""
import os, sys, glob, pickle
import numpy as np
from fractions import Fraction
from math import gcd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from store import Store, parse, val, read as read_out
from gen import words, cont, rationals

def cf(x):
    a = []; p, q = x.numerator, x.denominator
    while q: a.append(p // q); p, q = q, p % q
    return tuple(a[1:])

store = Store()
complete8 = [cf(x) for x in rationals(8)]
roots = [()] + words(14, False)
ctx = [w + ('#',) + u for w in complete8 for u in [()] + words(6, False)]
Pc = roots + ctx
Sc = [cf(x) for x in rationals(20)] + [('#',) + cf(x) for x in rationals(12)] + \
     [cf(a) + ('#',) + cf(b) for a in rationals(5) for b in rationals(8)]
def fill(rows, cols):
    store.request([parse(r + c) for r in rows for c in cols])
    store.compute()
    return np.array([[store.F(parse(r + c)) for c in cols] for r in rows], dtype=float)
H = fill(Pc, Sc)
print('Hankel %d x %d, missing %d' % (*H.shape, np.isnan(H).sum()))
H = np.nan_to_num(H)
Uu, Ss, Vt = np.linalg.svd(H, full_matrices=False)
print('singular values:', ' '.join('%.1e' % x for x in Ss[:80:4]), '(every fourth)')

def pivot(M, m):
    R = M.copy(); sel = []
    for _ in range(m):
        j = int(np.argmax(np.linalg.norm(R, axis=0))); sel.append(j)
        nrm = np.linalg.norm(R[:, j])
        if nrm == 0: break
        v = R[:, j] / nrm; R = R - np.outer(v, v @ R); R[:, j] = 0
    return sel

Kmax = max(map(int, sys.argv[1:]))
m = Kmax + 10
Psel = [Pc[i] for i in pivot(Uu[:, :Kmax].T, m)]
Ssel = [Sc[j] for j in pivot(Vt[:Kmax], m)]
if () not in Psel: Psel = [()] + Psel
syms = list(range(1, 13)) + ['#']
Hb = fill(Psel, Ssel)
Hs = {s: fill([p + (s,) for p in Psel], Ssel) for s in syms}
He = fill([()], Ssel)[0]                 # empty prefix
Hp = fill(Psel, [()])[:, 0]              # empty suffix: F(p) itself
print('basis %d x %d; shifted blocks done' % Hb.shape)

# Validation data
T = os.environ['BULB_DATA'] + '/tree/'
d2 = []
for f in glob.glob(T + 'd2/*.out'):
    r1 = Fraction(*map(int, os.path.basename(f)[:-4].split('-')))
    for x, v in read_out(f).items(): d2.append((cf(r1) + ('#',) + cf(x), v[0]))
d3 = []
for f in glob.glob(T + 'd3/*.out'):
    a, b = os.path.basename(f)[:-4].split('.')
    r1 = Fraction(*map(int, a.split('-'))); r2 = Fraction(*map(int, b.split('-')))
    for x, v in read_out(f).items(): d3.append((cf(r1) + ('#',) + cf(r2) + ('#',) + cf(x), v[0]))
root = [(cf(x), store.F((x,))) for x in rationals(60) if x.denominator > 20]

for K in map(int, sys.argv[1:]):
    U, S, Vh = np.linalg.svd(np.nan_to_num(Hb), full_matrices=False)
    P = U[:, :K] * S[:K]; Q = Vh[:K]; Pp = np.linalg.pinv(P); Qp = np.linalg.pinv(Q)
    N = {s: Pp @ np.nan_to_num(Hs[s]) @ Qp for s in syms}
    al = np.nan_to_num(He) @ Qp; be = Pp @ np.nan_to_num(Hp)
    def Fm(w):
        if any(s not in N for s in w): return None
        v = al
        for s in w: v = v @ N[s]
        return v @ be
    out = []
    for nm, data in (('root q 21..60', root), ('depth 2', d2), ('depth 3', d3)):
        e = [abs(Fm(w) - v) for w, v in data if Fm(w) is not None and v is not None]
        out.append('%s n=%d median %.1e 90%% %.1e max %.1e' % (nm, len(e), np.median(e), np.quantile(e, .9), max(e)))
    print('K=%d: %s' % (K, ' | '.join(out)))
    sys.stdout.flush()
