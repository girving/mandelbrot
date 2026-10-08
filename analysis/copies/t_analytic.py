"""Is Φ(t) = 256^n area(n, w0) analytic in t = 4^-n down to n = 0?  Fit polynomials in t on the deep members and
predict the shallow ones (n = 0, 1, 2)."""
import sys, math, collections
import numpy as np
from load import load, nrp
roots = load(sys.argv[1])
fam = collections.defaultdict(dict)
for r in roots:
    if not nrp(r) or abs(r.c.imag) > 1e-12 or not (-2 < r.c.real < -1.7): continue
    c = r.c.real; p = r.p
    z = [0.0]
    for i in range(p): z.append(z[-1] ** 2 + c)
    bc = (1 + math.sqrt(1 - 4 * c)) / 2
    k = 2; n = 0
    while k < p and abs(z[k] - bc) < 0.25: n += 1; k += 1
    itin = ''.join('R' if z[j] > 0 else 'L' for j in range(k, p))
    fam[itin][n] = r.area
for itin in sorted(fam, key=lambda s: (len(s), s)):
    d = fam[itin]
    if len(d) < 9: continue
    ns = np.array(sorted(d)); phi = np.array([d[n] * 256.0 ** n for n in ns]); t = 4.0 ** -ns
    print('family %-6s (m=%d): n = %d..%d, Φ(n=0) = %.10e, Φ(deepest) = %.10e' % (itin, len(itin), ns[0], ns[-1], phi[0], phi[-1]))
    for lo in (3, 4, 5):
        deep = ns >= lo
        for deg in range(2, min(7, deep.sum() - 1)):
            c = np.polyfit(t[deep], phi[deep], deg)
            pred = np.polyval(c, t[~deep])
            err = (pred - phi[~deep]) / phi[~deep]
            print('   fit on n ≥ %d, degree %d in t: predicted Φ at n = %s rel errors %s' % (lo, deg, list(ns[~deep]), ' '.join('%+.1e' % e for e in err)))
