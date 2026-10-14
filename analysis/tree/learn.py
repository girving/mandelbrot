"""Joint automaton over digits and a separator '#' for the satellite tree's normalized areas F(address).

  BULB_DATA=dir python3 learn.py K [K ...]

F(r1 # r2 # ... # rm) = αᵀ N(digits of r1) N(#) N(digits of r2) N(#) ... N(digits of rm but the last) β(last digit).
Hankel rows: root prefixes u (q_u ≤ 20) and context prefixes r1#u (r1 ≤ 1/2, q1 ≤ 8; q_u ≤ 8); columns: digit
suffixes s (q_s ≤ 40) and separator columns #s' (s' rationals, q ≤ 20).  N(b) by column closure on digit suffixes,
N(#) by column closure on separator columns, α from the empty prefix's row, β from single-digit columns.
Validation on depth-2 data (40 parents q1 ≤ 16, children q2 ≤ 30) and depth-3 data (tree/d3)."""
import os, sys, pickle, glob
import numpy as np
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen import val, cont, words, rationals, U20, U8, S40, S20, R1, D

def cf(x):
    a = []; p, q = x.numerator, x.denominator
    while q: a.append(p // q); p, q = q, p % q
    return tuple(a[1:])

def read(path):
    F = {}
    for line in open(path):
        t = line.split()
        if t[2] != 'failed':
            p, q = int(t[0]), int(t[1]); F[Fraction(p, q)] = float(t[5]) + (float(t[10]) if len(t) > 10 else 0.0)
    return F

def card_table():
    F = {}
    for n in ('all_e.out', 'stage2.out', 'stage3.out'):
        for x, f in read(os.environ['BULB_DATA'] + '/' + n).items(): F[x] = F[1 - x] = f
    return F

def node_F(addr, Fcard, ext):
    """F of address (tuple of rationals); ext: {file key: table}"""
    if len(addr) == 1: return Fcard.get(addr[0])
    if len(addr) == 2:
        t = ext.get('A_%d-%d' % (addr[0].numerator, addr[0].denominator)) if addr[1].denominator <= 20 else None
        v = t.get(addr[1]) if t else None
        if v is None:
            t = ext.get('B_%d-%d' % (addr[0].numerator, addr[0].denominator))
            v = t.get(addr[1]) if t else None
        return v
    if len(addr) == 3:
        t = ext.get('C_%d-%d_%d-%d' % (addr[0].numerator, addr[0].denominator, addr[1].numerator, addr[1].denominator))
        return t.get(addr[2]) if t else None

def build_hankel():
    Fcard = card_table()
    ext = {os.path.basename(f)[:-4]: read(f) for f in glob.glob(D + 'ext/*.out')}
    rows = [('root', u) for u in U20 if u and val(list(u)) != 1] + \
           [('ctx', r1, u) for r1 in R1 for u in U8 if u and val(list(u)) != 1]
    cols = [('dig', s) for s in S40] + [('sep', s2) for s2 in S20]
    def addr(row, col):
        if row[0] == 'root':
            if col[0] == 'dig': return (val(list(row[1]) + list(col[1])),)
            return (val(list(row[1])), col[1])
        if col[0] == 'dig': return (row[1], val(list(row[2]) + list(col[1])))
        return (row[1], val(list(row[2])), col[1])
    H = np.full((len(rows), len(cols)), np.nan)
    for i, r in enumerate(rows):
        for j, c in enumerate(cols):
            v = node_F(addr(r, c), Fcard, ext)
            if v is not None: H[i, j] = v
    return rows, cols, H, Fcard, ext

if __name__ == '__main__':
    rows, cols, H, Fcard, ext = build_hankel()
    good_r = ~np.isnan(H).any(1); good_c = ~np.isnan(H).any(0)
    print('Hankel %d x %d; complete rows %d, complete cols %d' % (*H.shape, good_r.sum(), good_c.sum()))
    rows = [r for r, g in zip(rows, good_r) if g]; H = H[good_r]
    keep_c = ~np.isnan(H).any(0); cols = [c for c, g in zip(cols, keep_c) if g]; H = H[:, keep_c]
    print('using %d x %d' % H.shape)
    Uu, Ss, Vt = np.linalg.svd(H, full_matrices=False)
    print('singular values:', ' '.join('%.1e' % x for x in Ss[:60:3]), '(every third)')
    cidx = {c: j for j, c in enumerate(cols)}
    # Validation sets: depth 2 (tree/d2: parents r1 q1 ≤ 16), depth 3 (tree/d3)
    T = D
    d2 = []
    for f in glob.glob(T + 'd2/*.out'):
        p1, q1 = map(int, os.path.basename(f)[:-4].split('-')); r1 = Fraction(p1, q1)
        for x, v in read(f).items():
            if max(cf(x)) <= 12 and max(cf(r1)) <= 12: d2.append(((r1, x), v))
    d3 = []
    for f in glob.glob(T + 'd3/*.out'):
        a, b = os.path.basename(f)[:-4].split('.')
        r1 = Fraction(*map(int, a.split('-'))); r2 = Fraction(*map(int, b.split('-')))
        for x, v in read(f).items():
            if max(cf(x)) <= 12: d3.append(((r1, r2, x), v))
    for K in map(int, sys.argv[1:]):
        P = Uu[:, :K] * Ss[:K]; Q = Vt[:K]; Pp = np.linalg.pinv(P)
        N = {}
        for b in range(1, 13):
            pairs = [(cidx[('dig', s)], cidx[('dig', (b,) + s)]) for s in S40 if ('dig', s) in cidx and ('dig', (b,) + s) in cidx]
            if len(pairs) < 3: continue
            N[b] = Pp @ H[:, [p[1] for p in pairs]] @ np.linalg.pinv(Q[:, [p[0] for p in pairs]])
        pairs = [(cidx[('dig', cf(s2))], cidx[('sep', s2)]) for s2 in S20 if ('sep', s2) in cidx and ('dig', cf(s2)) in cidx]
        Nsep = Pp @ H[:, [p[1] for p in pairs]] @ np.linalg.pinv(Q[:, [p[0] for p in pairs]])
        # α from the empty prefix: F(s) for digit columns
        dig = [(j, c[1]) for c, j in cidx.items() if c[0] == 'dig']
        y = np.array([Fcard[val(list(s))] for j, s in dig])
        al = y @ np.linalg.pinv(Q[:, [j for j, s in dig]])
        be = {s[0]: Q[:, cidx[('dig', s)]] for s in S40 if len(s) == 1 and ('dig', s) in cidx}
        def F_model(addr):
            v = al.copy()
            for li, r in enumerate(addr):
                w = cf(r)
                if any(a not in N for a in w[:-1]) or (li == len(addr) - 1 and w[-1] not in be): return None
                if li < len(addr) - 1:
                    if w[-1] not in N: return None
                    for a in w: v = v @ N[a]
                    v = v @ Nsep
                else:
                    for a in w[:-1]: v = v @ N[a]
                    return v @ be[w[-1]]
        out = []
        for name, data in (('depth 2', d2), ('depth 3', d3)):
            e = [abs(F_model(a) - v) for a, v in data if F_model(a) is not None]
            out.append('%s n=%d median %.1e 90%% %.1e max %.1e' % (name, len(e), np.median(e), np.quantile(e, .9), max(e)))
        print('K=%d: %s' % (K, ' | '.join(out)))
