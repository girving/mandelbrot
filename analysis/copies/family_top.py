"""High-precision constants for the heaviest families: C = lim k^4 area by Neville in 1/k (Decimal) on two subsets of
the k grid (the spread is the error estimate).

  python3 family_top.py area.txt[.gz] > constants   (lines "name C error")"""
import sys, gzip
from collections import defaultdict
from decimal import Decimal as D, getcontext
getcontext().prec = 50

def neville(xs, ys):
    p = list(ys)
    for m in range(1, len(xs)):
        p = [(xs[i+m] * p[i] - xs[i] * p[i+1]) / (xs[i+m] - xs[i]) for i in range(len(p) - 1)]
    return p[0]

a = defaultdict(dict)
path = sys.argv[1]
for l in (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)):
    r = l.split()
    if r[3] == 'failed': continue
    name, k = r[0].rsplit('_', 1)
    a[name][int(k)] = D(float(r[5])) + D(float(r[10]))
for name, d in sorted(a.items()):
    ks = sorted(d)
    big = [k for k in ks if k >= 48]
    s1 = big[::max(1, len(big) // 12)][-12:]
    s2 = big[1::max(1, len(big) // 12)][-12:]
    est = [neville([D(1) / k for k in s], [d[k] * D(k)**4 for k in s]) for s in (s1, s2)]
    print('%s %.18e %.1e' % (name, est[0], abs(est[0] - est[1]) / abs(est[0])))
