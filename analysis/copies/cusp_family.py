"""Copies near the airplane's cusp c = -1.75 (period-3 saddle-node), on the intermittent side c > -1.75: count the
critical orbit's passes through the period-3 gate (runs of |z_{i+3} - z_i| small) and the area against the
distance to the cusp."""
import sys, math, collections
import numpy as np
from load import load, nrp
roots = load(sys.argv[1])
rows = []
for r in roots:
    if not nrp(r) or abs(r.c.imag) > 1e-12 or not (-1.75 < r.c.real < -1.70): continue
    c = r.c.real; p = r.p
    z = [0.0]
    for i in range(p): z.append(z[-1] ** 2 + c)
    near = [abs(z[i + 3] - z[i]) < 0.02 for i in range(1, p - 2)]
    longest = 0; run = 0
    for v in near:
        run = run + 1 if v else 0; longest = max(longest, run)
    rows.append((c, p, r.area, longest))
rows.sort()
print('%d real NRPs in (-1.75, -1.70); the ten largest with c → -1.75:' % len(rows))
for c, p, a, L in sorted(rows, key=lambda t: -t[2])[:10]:
    print('  c=%.10f p=%2d area %.3e  longest gate run %d   (c+1.75)=%.3e  area/(c+1.75)^3 = %.3e' % (c, p, a, L, c + 1.75, a / (c + 1.75) ** 3))
# Closest-to-cusp family: the largest copy of each period
best = {}
for c, p, a, L in rows:
    if p not in best or a > best[p][2]: best[p] = (c, p, a, L)
print('largest copy per period:')
for p in sorted(best):
    c, p, a, L = best[p]
    print('  p=%2d c=%.12f area %.4e gate run %2d  d=c+1.75 %.4e  log area / log d %.3f' % (p, c, a, L, c + 1.75, math.log(a) / math.log(c + 1.75)))
