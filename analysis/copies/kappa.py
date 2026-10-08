import sys, math, collections
import numpy as np
from load import load, nrp, A_CARD
roots = load(sys.argv[1])
print('%d roots; NRPs %d, satellite letters %d' % (len(roots), sum(nrp(r) for r in roots), sum(r.parent < 0 and r.sat for r in roots)))
# NRP area by period
by = collections.defaultdict(float)
for r in roots:
    if nrp(r): by[r.p] += r.area
print('NRP area by period: ' + ' '.join('%d:%.2e' % (p, by[p]) for p in sorted(by)))
print('NRP total through 16: %.6e;  p^3 * per-period NRP area: %s' % (sum(by.values()), ' '.join('%.3f' % (p ** 3 * by[p]) for p in sorted(by))))
# κ for tunings, by kind of W0
kap = {'primitive W0': [], 'satellite W0': []}
for r in roots:
    if r.parent < 0 or r.v is None: continue
    w0 = roots[r.parent]
    k = r.area * A_CARD / (w0.area * r.v.area)
    kap['satellite W0' if w0.sat else 'primitive W0'].append((k, r, w0))
for name, v in kap.items():
    ks = np.array([x[0] for x in v])
    print('%s: %d tunings; κ = area(W0⋆V) A_card/(area(W0) area(V)): median %.4f, range [%.4f, %.4f]' % (name, len(ks), np.median(ks), ks.min(), ks.max()))
# Primitive copies: δ = κ - 1 by W0 (copy) and by V's kind
print('\nprimitive copies: δ = κ - 1 per copy W0 (period, center), over its components V:')
per = collections.defaultdict(list)
for k, r, w0 in kap['primitive W0']: per[w0.i].append((k - 1, r.v))
for i in sorted(per, key=lambda i: (roots[i].p, -roots[i].area))[:12]:
    w0 = roots[i]; d = np.array([x[0] for x in per[i]])
    vs = per[i]
    near_root = [x[0] for x in vs if x[1].sat and x[1].parent < 0]
    print('  W0 p=%2d c=%+.5f%+.5fi area %.2e: %3d V, δ median %+.2e range [%+.2e, %+.2e]; V = cardioid bulbs: median %+.2e'
          % (w0.p, w0.c.real, w0.c.imag, w0.area, len(d), np.median(d), d.min(), d.max(), np.median(near_root) if near_root else float('nan')))
