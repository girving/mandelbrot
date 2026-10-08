"""Family constants C(word) = lim_k k^4 area(k, word) from `bulb_batch --local` output over the k of family_jobs.py.

  python3 family_constants.py area.txt[.gz] > constants   (lines "name C error")

g(k) = k^4 area (sin πt / πt)^6 with t = n/(2k), n the word's largest digit (a run of n kneading symbols is a pass of
~n/2 turns: the gate model's shape divided out), fitted as a polynomial in t over k ≥ 2K (K the word's least k);
the limit is the constant term, and the error estimate the spread over polynomial degrees 4..7."""
import sys, gzip, math
from collections import defaultdict
from decimal import Decimal as D
import warnings
import numpy as np
warnings.simplefilter("ignore", np.exceptions.RankWarning)

def load(path):
    a = defaultdict(dict)
    for l in (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)):
        r = l.split()
        if r[3] == 'failed': continue
        name, k = r[0].rsplit('_', 1)
        a[name][int(k)] = float(D(float(r[5])) + D(float(r[10])))
    return a

def constant(d, n):
    K = min(d)
    ks = sorted(k for k in d if k >= 2 * K)
    if len(ks) < 10: return None, None  # a truncated family (the size engine stopped early)
    t = np.array([n / (2 * k) for k in ks])
    g = np.array([d[k] * k**4 * (math.sin(math.pi * x) / (math.pi * x))**6 for k, x in zip(ks, t)])
    fits = [np.polynomial.polynomial.polyfit(t, g, deg)[0] for deg in range(4, 8)]
    return fits[-2], max(fits) - min(fits)

if __name__ == '__main__':
    for name, d in sorted(load(sys.argv[1]).items()):
        n = max(int(x) for x in name.split('_')[1:])
        C, e = constant(d, n)
        if C is None:
            print('skipped %s: too few large-k points' % name, file=sys.stderr)
            continue
        print('%s %.15e %.1e' % (name, C, abs(e / C)))
