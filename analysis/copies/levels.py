"""At fixed lingering depth n near -2, is Φ_n(w0*) = 256^n area |(f_a^m)'(w0*)|^4 smooth in w0* across all copies?"""
import sys, math, collections
import numpy as np
from load import load, nrp
roots = load(sys.argv[1])
lev = collections.defaultdict(list)
for r in roots:
    if not nrp(r) or abs(r.c.imag) > 1e-12 or not (-2 < r.c.real < -1.7): continue
    c = r.c.real; p = r.p
    z = [0.0]
    for i in range(p): z.append(z[-1] ** 2 + c)
    bc = (1 + math.sqrt(1 - 4 * c)) / 2
    k = 2; n = 0
    while k < p and abs(z[k] - bc) < 0.25: n += 1; k += 1
    itin = ''.join('R' if z[j] > 0 else 'L' for j in range(k, p))
    ws = 0.0
    for ch in reversed(itin): ws = math.sqrt(ws + 2) * (1 if ch == 'R' else -1)
    der = 1.0; x = ws
    for j in range(len(itin)): der *= 2 * x; x = x * x - 2
    lev[n].append((ws, len(itin), r.area * 256.0 ** n * abs(der) ** 4, z[k]))
for n in sorted(lev):
    v = sorted(lev[n]); xs = np.array([a[0] for a in v]); ys = np.log(np.array([a[2] for a in v]))
    if len(v) < 12: continue
    out = []
    for deg in (4, 8, 12):
        if len(v) > deg + 3:
            cf = np.polynomial.chebyshev.chebfit(xs, ys, deg); res = ys - np.polynomial.chebyshev.chebval(xs, cf)
            out.append('deg %d rms %.1e' % (deg, np.sqrt((res ** 2).mean())))
    print('n=%2d: %4d copies, w0* in [%.3f, %.3f], log Φ_n range [%.2f, %.2f]; smooth fits: %s' % (n, len(v), xs.min(), xs.max(), ys.min(), ys.max(), ', '.join(out)))
