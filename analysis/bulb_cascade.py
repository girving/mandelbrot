#!/usr/bin/env python3
"""Analyses of satellite bulb areas (bulb_areas output) for the renormalization cascade.

  bulb_cascade.py lists DIR            write DIR/{all128,all40,fade}.txt: "p q" request lists for bulb_areas
  bulb_cascade.py fade DIR             memory of the first continued fraction digit (needs DIR/fade.out)
  bulb_cascade.py levels DIR           bulbs of other parents against the cardioid's (needs DIR/card.out, DIR/lv_*.out)
  bulb_cascade.py psi DIR              F(p/q) against y = (p^-1 mod q)/q, and the per-q average F̄(q)
  bulb_cascade.py modes DIR            spectral expansion of the fading (needs DIR/modes.out, from DIR/modes.txt)

bulb_areas prints "p q center_re center_im area F |c_W'|^2 conv root_err"; F = area q^4 / (π |c_W'(λ0)|^2) is
the bulb's area normalized by its parent's multiplier map at the root.
"""
from fractions import Fraction
from math import gcd
import os
import random
import sys

import numpy as np


def cf_value(w):
    x = Fraction(0)
    for a in reversed(w):
        x = 1 / (a + x)
    return x


def fade_suffixes():
    random.seed(5)
    rnd = [random.choice([1, 1, 2, 3, 4]) for _ in range(20)]
    out = {}
    for L in range(13):
        out[('ones', L)] = [1] * L + [2]
    for L in range(8):
        out[('twos', L)] = [2] * L + [2]
    for L in range(11):
        out[('rand', L)] = rnd[:L] + [3]
    return out


FIRST = [1, 2, 3, 5, 10]


# Modes: prefixes u in front of periodic suffix families s_L = (repeated block)^L + [last digit]
MODE_PREFIXES = [[1], [2], [3], [4], [5], [1, 1], [1, 2], [2, 1], [2, 2], [1, 1, 1], [3, 1], [1, 3]]
MODE_FAMILIES = {'1': ([1], 2, 17), '2': ([2], 2, 10), '3': ([3], 2, 7), '12': ([1, 2], 2, 8), '4': ([4], 3, 6)}
MODE_QMAX = 40000


def mode_words():
    """{(family, L, prefix index): word}, keeping only complete rows (every prefix with q ≤ MODE_QMAX)"""
    out = {}
    for name, (block, last, Lmax) in MODE_FAMILIES.items():
        for L in range(Lmax + 1):
            row = {i: u + block * L + [last] for i, u in enumerate(MODE_PREFIXES)}
            if all(cf_value(w).denominator <= MODE_QMAX for w in row.values()):
                for i, w in row.items():
                    out[(name, L, i)] = w
    return out


def read(path):
    F = {}
    for line in open(path):
        t = line.split()
        if t[2] != 'failed':
            F[(int(t[0]), int(t[1]))] = (float(t[4]), float(t[5]))
    return F


def lists(d):
    with open(os.path.join(d, 'all128.txt'), 'w') as f:
        for q in range(2, 129):
            for p in range(1, q // 2 + 1):
                if gcd(p, q) == 1:
                    print(p, q, file=f)
    with open(os.path.join(d, 'all40.txt'), 'w') as f:
        for q in range(2, 41):
            for p in range(1, q):
                if gcd(p, q) == 1:
                    print(p, q, file=f)
    seen = set()
    with open(os.path.join(d, 'modes.txt'), 'w') as f:
        for w in mode_words().values():
            x = cf_value(w)
            if x not in seen:
                seen.add(x)
                # Both p/q and its mirror (q-p)/q: equal areas by symmetry, so their difference measures the error
                print(x.numerator, x.denominator, file=f)
                print(x.denominator - x.numerator, x.denominator, file=f)
    seen = set()
    with open(os.path.join(d, 'fade.txt'), 'w') as f:
        for s in fade_suffixes().values():
            for a in FIRST:
                x = cf_value([a] + s)
                if x.denominator <= 20000 and x not in seen:
                    seen.add(x)
                    print(x.numerator, x.denominator, file=f)


def fade(d):
    F = read(os.path.join(d, 'fade.out'))

    def f(w):
        x = cf_value(w)
        v = F.get((x.numerator, x.denominator)) or F.get(((1 - x).numerator, x.denominator))
        return v[1] if v else None

    fits = {}
    print('suffix family, L: spread of F over first digits %s, mean, shape (F(a)-F(10))/(F(1)-F(10)) for a = 2,3,5'
          % FIRST)
    for (name, L), s in fade_suffixes().items():
        v = [f([a] + s) for a in FIRST]
        if None in v:
            continue
        v = np.array(v)
        shape = (v[1:4] - v[4]) / (v[0] - v[4])
        print('  %-4s L=%2d  spread %.2e  mean %.8f  shape %s' % (name, L, np.ptp(v), v.mean(),
                                                                  ' '.join('%.4f' % x for x in shape)))
        if L >= 2:
            fits.setdefault(name, []).append((np.log(cf_value(s).denominator), np.log(np.ptp(v))))
    for name, pts in fits.items():
        pts = np.array(pts)
        print('  %s: spread ∝ q_suffix^%.3f' % (name, np.polyfit(pts[:, 0], pts[:, 1], 1)[0]))


def levels(d):
    card = read(os.path.join(d, 'card.out'))
    Fc = {}
    for (p, q), (A, F) in card.items():
        Fc[(p, q)] = Fc[(q - p, q)] = F
    print('R = F / F_cardioid for the bulbs of other parents (q ≤ 40):')
    print('  parent (period, center)          mid |R-1| median  max     R(1/n), n = 5,10,20,40       R((n-1)/n)')
    for name in sorted(os.listdir(d)):
        if not name.startswith('lv_'):
            continue
        r = {k: v[1] / Fc[k] for k, v in read(os.path.join(d, name)).items()}
        mid = [abs(v - 1) for (p, q), v in r.items() if min(p, q - p) / q >= 0.25 and q >= 20]
        fmt = lambda v: ' '.join('%.4f' % x if x else '  -   ' for x in v)
        print('  %-32s %.1e         %.1e  %s   %s' % (name[3:-4], np.median(mid), max(mid),
                                                    fmt([r.get((1, n)) for n in (5, 10, 20, 40)]),
                                                    fmt([r.get((n - 1, n)) for n in (5, 10, 20, 40)])))


def psi(d):
    card = read(os.path.join(d, 'card.out'))
    p = np.array([k[0] for k in card]); q = np.array([k[1] for k in card])
    A = np.array([v[0] for v in card.values()]); F = np.array([v[1] for v in card.values()])
    x = np.minimum(p / q, 1 - p / q)
    y = np.array([pow(int(a), -1, int(b)) / b for a, b in zip(p, q)])
    y = np.minimum(y, 1 - y)
    m = q >= 100
    bins = np.linspace(0, 0.5, 51)
    for name, v in [('x = p/q', x), ('y = (p^-1 mod q)/q', y)]:
        idx = np.digitize(v[m], bins)
        sp = [np.ptp(F[m][idx == k]) for k in range(1, 51) if (idx == k).sum() > 5]
        print('median spread of F in 50 bins of %-20s (q ≥ 100): %.4f' % (name, np.median(sp)))
    # Per-q average weighted like the sum: Σ_q (π/q^4) W_q F̄(q), W_q = Σ_p sin^2(πp/q) over both half planes
    w2 = np.where(2 * p == q, 1, 2)
    S = np.sin(np.pi * p / q) ** 2
    W = {Q: (w2[q == Q] * S[q == Q]).sum() for Q in range(2, q.max() + 1)}
    Fb = {Q: (w2[q == Q] * A[q == Q]).sum() * Q ** 4 / np.pi / W[Q] for Q in W}
    print('F̄(q):', ' '.join('%d:%.5f' % (Q, Fb[Q]) for Q in (8, 16, 32, 64, 96, 127, 128) if Q in Fb))
    qs = np.arange(16, 65)
    X = np.vstack([np.ones(len(qs)), qs ** -0.54]).T
    c = np.linalg.lstsq(X, np.array([Fb[Q] for Q in qs]), rcond=None)[0]
    truth = (w2 * A)[q > 64].sum()
    for name, fn in [('constant F̄', lambda Q: np.mean([Fb[Q] for Q in qs])),
                     ('drift %.4f %+.4f q^-0.54' % tuple(c), lambda Q: c[0] + c[1] * Q ** -0.54)]:
        pred = sum(np.pi / Q ** 4 * W[Q] * fn(Q) for Q in range(65, q.max() + 1))
        print('tail Σ_{64<q≤%d} A = %.6e from q ≤ 64 data, %s: rel err %.1e' % (q.max(), truth, name,
                                                                                (pred - truth) / truth))
    for Q in (16, 32, 64, 128):
        print('cardioid bulbs through q = %3d: %.12f' % (Q, (w2 * A)[q <= Q].sum()))


def modes(d):
    """Matrix pencil: D[L, u] = F(u s_L) - F(u_0 s_L) = Σ_k a_k λ_k^L φ_k(u) for each suffix family"""
    F = read(os.path.join(d, 'modes.out'))
    err = max(abs(F[(p, q)][1] - F[(q - p, q)][1]) for (p, q) in F if (q - p, q) in F)
    print('max |F(p/q) - F((q-p)/q)| (accuracy): %.1e' % err)
    words = mode_words()
    for name, (block, last, Lmax) in MODE_FAMILIES.items():
        Ls = sorted({L for (n, L, i) in words if n == name})
        M = np.array([[F[(cf_value(words[(name, L, i)]).numerator, cf_value(words[(name, L, i)]).denominator)][1]
                       for i in range(len(MODE_PREFIXES))] for L in Ls])
        D = M[:, 1:] - M[:, :1]
        sv = np.linalg.svd(D, compute_uv=False)
        print('family %s^L + [%d], L = %d..%d: singular values of D %s' % (block, last, Ls[0], Ls[-1],
                                                                       ' '.join('%.1e' % x for x in sv[:6])))
        print('  row norms |D[L]|: %s' % ' '.join('%.1e' % np.linalg.norm(r) for r in D))
        # Contraction of the Gauss map's inverse branches along the periodic orbit [block]^∞: Π y_i^2 over a block
        y = [0.5] * len(block)
        for _ in range(200):
            for j in reversed(range(len(block))):
                y[j] = 1 / (block[j] + y[(j + 1) % len(block)])
        contraction = np.prod(np.array(y) ** 2)
        best = None
        for K in (1, 2, 3, 4):
            for L0 in (0, 2):
                X0, X1 = D[L0:-1], D[L0 + 1:]
                if len(X0) <= K:
                    continue
                U, S, Vt = np.linalg.svd(X0, full_matrices=False)
                Z = U[:, :K].T @ X1 @ Vt[:K].T @ np.diag(1 / S[:K])
                lam = np.linalg.eigvals(Z)
                lam = lam[np.argsort(-abs(lam))]
                # Residual of the K-mode fit: least squares of D on the modes λ^L
                V = np.array([[l ** L for l in lam] for L in range(L0, len(D))])
                coef = np.linalg.lstsq(V, D[L0:], rcond=None)[0]
                res = np.abs(V @ coef - D[L0:]).max()
                if L0 == 2 and len(D) - L0 > K + 1 and (best is None or res < best[0]):
                    best = (res, K, lam)
                print('  K=%d from L=%d: λ = %s   max residual %.1e' % (
                    K, Ls[L0], '  '.join('%.4f%+.4fi' % (l.real, l.imag) if abs(l.imag) > 1e-6 else '%.4f' % l.real
                                         for l in lam), res))
        res, K, lam = best
        print('  block contraction Π y_i^2 = %.5f; best fit (K=%d from L=2): γ = log|λ| / log(Π y_i^2) = %s' % (
            contraction, K, ' '.join('%.3f' % (np.log(abs(l)) / np.log(contraction)) for l in lam)))


if __name__ == '__main__':
    {'lists': lists, 'fade': fade, 'levels': levels, 'psi': psi, 'modes': modes}[sys.argv[1]](sys.argv[2])
