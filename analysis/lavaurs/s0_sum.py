"""The one-transit sector S_0 = Σ_w C_w of the seahorse families from Lavaurs-model constants.

  python3 s0_sum.py model_census walk_classes model_walk camp_C

Census (digit sum ≤ 13) from the model run on all 16,383 census families; beyond it, per skeleton class (walk_classes,
lavaurs_area --walk output named C<i>_<n>, n = excursion = digit sum), members with n ≥ 14 from the walk, the campaign
constant (M side) where the walk stopped short, and past the last member the fit C(n) n^6 = Σ_{i<I} a_i n^{-i/2} on the
computed members (spread over I and fit ranges = error).  Classes past 800 and small-digit words beyond the census are
the estimates of the M-side sum (6.6e-13, ~1e-13)."""
import sys, gzip
from collections import defaultdict
import numpy as np

def tail(ns, Cs, n_from, I, N0):
    pts = [(n, c) for n, c in zip(ns, Cs) if n >= N0]
    if len(pts) < I + 3: return None
    x = np.array([n for n, _ in pts], float); y = np.array([c * n**6 for n, c in pts])
    A = np.array([[n**(-i / 2) for i in range(I)] for n in x]); a = np.linalg.lstsq(A, y, rcond=None)[0]
    n = np.arange(n_from, 400000, 2, dtype=float)
    return float(sum(a[i] * np.sum(n**(-6 - i / 2)) for i in range(I)) + a[0] * 400000.0**-5 / 10)

if __name__ == '__main__':
    census = sum(float(l.split()[7]) for l in gzip.open(sys.argv[1], 'rt') if l.split()[3] != 'failed')
    classes = {}
    for l in open(sys.argv[2]):
        p = l.split(); classes[p[0]] = (p[1:-1], int(p[-1]))
    walk = defaultdict(dict)
    for l in gzip.open(sys.argv[3], 'rt'):
        r = l.split(); nm, n = r[0].rsplit('_', 1); walk[nm][int(n)] = float(r[7])
    camp = {l.split()[0]: float(l.split()[1]) for l in open(sys.argv[4])}
    beyond = spread_tot = 0.0; used_camp = 0
    for ci, (sk, par) in classes.items():
        if ci not in walk: continue
        others = sum(int(x) for x in sk[1:] if x != '*')
        mem = dict(walk[ci])
        # campaign members (M side) past the walk's reach: class word for digit m = n - others
        last = max(mem)
        for n in range(last + 2, 66 + others, 2):
            word = '_'.join([sk[0]] + [str(n - others) if x == '*' else x for x in sk[1:]])
            if word in camp: mem[n] = camp[word]; used_camp += 1
            else: break
        ns = sorted(mem); Cs = [mem[n] for n in ns]
        direct = sum(c for n, c in zip(ns, Cs) if n >= 14)
        ests = [t for I in (3, 4, 5) for N0 in (16, 24, 32) for t in [tail(ns, Cs, ns[-1] + 2, I, N0)] if t is not None]
        t = sorted(ests)[len(ests) // 2] if ests else 0.0
        spread = (max(ests) - min(ests)) if ests else 0.0
        beyond += direct + t; spread_tot += spread
    S0 = census + beyond + 6.6e-13 + 1e-13
    print('census (model) %.15e' % census)
    print('beyond the census: %.12e from %d classes (%d campaign members; tail spread %.1e)' % (beyond, len(walk), used_camp, spread_tot))
    print('S_0 = %.13e  (+ classes past 800 ≈ 6.6e-13, small-digit words ≈ 1e-13 included)' % S0)
