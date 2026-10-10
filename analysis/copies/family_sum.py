"""Σ_w C_w over all seahorse families (the phase-plane area at −3/4, times π²/16): census plus skeleton classes.

  python3 family_sum.py keys K0 census_constants campaign_constants

census_constants: limb_constants.py output (j i C err ...) for the census (tails ≤ 13 symbols = digit sum ≤ 13).
campaign_constants: family_constants.py output (name C err) for the skeleton classes of skeletons.py.
Beyond the census, each class (skeleton, parity) sums its members with digit sum ≥ 14 whose long digit n is the first
maximal digit: computed values while available (n ≤ 64), then the model C(n) n^6 = Σ_{i<I} a_i n^{-i/2} fitted on
n ≥ N0 (checked against the sparse points to n = 512).  The model's uncertainty is the spread over (I, N0); classes
outside the campaign are estimated from their census mass (the top classes' ratio of beyond-census to census mass)."""
import sys, math
from collections import defaultdict
import numpy as np
from digits import census
from skeletons import skeleton

def class_members(sk, par, Cw):
    """n -> C for the class's computed words"""
    out = {}
    for name, (C, e) in Cw.items():
        d = (name.split('_')[0],) + tuple(int(x) for x in name.split('_')[1:])
        if len(d) != len(sk): continue
        ok = all(a == b for a, b in zip(d, sk) if b != '*')
        if not ok: continue
        p = sk.index('*'); n = d[p]
        if n % 2 != par or skeleton(d) != (sk, par): continue
        out[n] = (C, e)
    return out

def model_sum(ns, Cs, n_from, I, N0, par):
    """Σ over n ≥ n_from (parity par) of the model fitted on the points with n ≥ N0"""
    pts = [(n, c) for n, c in zip(ns, Cs) if n >= N0]
    if len(pts) < I + 2: return None, None
    x = np.array([n for n, _ in pts], float); y = np.array([c * n**6 for n, c in pts])
    A = np.array([[n**(-i / 2) for i in range(I)] for n in x])
    a = np.linalg.lstsq(A, y, rcond=None)[0]
    resid = np.abs(A @ a - y).max() / abs(y).max()
    n = np.arange(n_from, 200000, 2, dtype=float)
    s = float(sum(a[i] * np.sum(n**(-6 - i / 2)) for i in range(I)))
    s += a[0] * 200000.0**-5 / 10  # Σ_{m ≥ N, step 2} m^-6 ≈ N^-5 / 10
    return s, resid

if __name__ == '__main__':
    cs = census(sys.argv[1], int(sys.argv[2]))
    Cc = {(int(r[0]), int(r[1])): float(r[2]) for r in (l.split() for l in open(sys.argv[3]))}
    Cw = {r[0]: (float(r[1]), float(r[2])) for r in (l.split() for l in open(sys.argv[4]))}
    S_census = sum(Cc.values())
    # class masses in the census (digit sum ≥ 9) and the campaign's classes
    mass = defaultdict(float)
    for d, (j, i, lo, hi, t) in cs.items():
        if len(d) >= 2 and sum(d[1:]) >= 9: mass[skeleton(d)] += Cc[(j, i)]
    classes = set()
    for name in Cw:
        d = (name.split('_')[0],) + tuple(int(x) for x in name.split('_')[1:])
        classes.add(skeleton(d))
    print('census Σ C = %.12e (%d families)' % (S_census, len(Cc)))
    total_beyond, budget, cmass, rows = 0.0, 0.0, 0.0, []
    for sk, par in sorted(classes, key=lambda c: -mass.get(c, 0)):
        mem = class_members(sk, par, Cw)
        if not mem: continue
        others = sum(x for x in sk[1:] if x != '*')
        n_min = max(14 - others, 0)
        ns = sorted(mem); Cs = [mem[n][0] for n in ns]
        direct = sum(mem[n][0] for n in ns if n >= n_min and n <= 64)
        last = max(n for n in ns if n <= 64)
        sums = []
        for I in (3, 4, 5):
            for N0 in (16, 24, 32):
                s, r = model_sum(ns, Cs, last + 2, I, N0, par)
                if s is not None: sums.append((s, r, I, N0))
        if not sums or len(ns) < 10:  # a class met only through stray small-n words of other classes: not summed
            continue
        best = sorted(sums, key=lambda x: x[1])[0]
        spread = max(s for s, *_ in sums) - min(s for s, *_ in sums)
        beyond = direct + best[0]
        total_beyond += beyond; budget += spread; cmass += mass.get((sk, par), 0)
        rows.append((sk, par, mass.get((sk, par), 0), direct, best, spread))
    for sk, par, m, direct, best, spread in rows[:12]:
        print('  %-18s par %d  census mass %.3e  beyond: n ≤ 64 %.4e + model %.4e (I=%d N0=%d, fit resid %.0e, spread %.0e)' %
              (' '.join(map(str, sk)), par, m, direct, best[0], best[2], best[3], best[1], spread))
    tot_mass = sum(mass.values())
    rest = total_beyond * (tot_mass - cmass) / cmass
    print('beyond the census: %.6e from %d classes (model spread %.1e), other classes ≈ %.1e (from their census share %.2e)' %
          (total_beyond, len(rows), budget, rest, (tot_mass - cmass) / tot_mass))
    S = S_census + total_beyond + rest
    print('Σ C ≈ %.12e ± %.1e   (phase-plane area 16 Σ/π² = %.10e)' % (S, budget + rest, 16 * S / math.pi**2))
