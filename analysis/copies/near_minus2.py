"""Primitive copies near the Misiurewicz point a = -2: labels (n, w0), size estimate, Ruelle weight.

For each non-renormalizable primitive component X with real center c in (-2, -1.7) (period p ≤ 16):
  orbit z_1 = c, z_{i+1} = z_i² + c, z_p = 0; n = steps the orbit spends near β_c = (1 + sqrt(1 - 4c))/2 ≈ 2;
  w0 = the first point after leaving β; m = steps from w0 to 0.  Size estimate s = 1/(β Λ²) with Λ = Π_{i<p} 2 z_i,
  β = Σ_{i=1}^{p-1} 1/Π_{j≤i} 2 z_j (the standard small-copy size estimate); its area model is A_card |s|².
Also the same w0 under f_a (a = -2): the precritical point of f_a with the same itinerary, and |(f_a^m)'(w0)|."""
import sys, math, collections
from load import load, nrp, A_CARD
roots = load(sys.argv[1])
rows = []
for r in roots:
    if not nrp(r) or abs(r.c.imag) > 1e-12 or not (-2 < r.c.real < -1.7): continue
    c = r.c.real; p = r.p
    z = [0.0]
    for i in range(p): z.append(z[-1] ** 2 + c)
    lam, beta, prod = 1.0, 0.0, 1.0
    for i in range(1, p):
        prod *= 2 * z[i]; beta += 1 / prod
    lam = prod
    s = 1 / (beta * lam * lam)
    bc = (1 + math.sqrt(1 - 4 * c)) / 2
    k = 2; n = 0
    while k < p and abs(z[k] - bc) < 0.25: n += 1; k += 1
    w0 = z[k] if k < p else float('nan'); m = p - k
    rows.append((p, c, r.area, s, n, w0, m))
rows.sort(key=lambda t: t[1])
print('%d real NRPs in (-2, -1.7)' % len(rows))
print(' p   center              area         area/(A_card s²)   n   w0         m')
for p, c, a, s, n, w0, m in rows[:40]:
    print('%2d  %.12f  %.4e   %.6f         %2d  %+.6f  %2d' % (p, c, a, a / (A_CARD * s * s), n, w0, m))
# Ratio statistics
rat = [a / (A_CARD * s * s) for p, c, a, s, n, w0, m in rows]
print('area/(A_card s²): min %.4f max %.4f' % (min(rat), max(rat)))
print('\nratio area/(A_card s²) by n (steps near β):')
by = collections.defaultdict(list)
for p, c, a, s, n, w0, m in rows: by[n].append(a / (A_CARD * s * s))
for n in sorted(by):
    v = sorted(by[n]); print('  n=%2d count %4d  median %.6f  10%%-90%% [%.6f, %.6f]  min %.4f max %.4f' % (n, len(v), v[len(v)//2], v[len(v)//10], v[9*len(v)//10], v[0], v[-1]))
print('\nlargest NRPs:')
for p, c, a, s, n, w0, m in sorted(rows, key=lambda t: -t[2])[:8]:
    print('  p=%2d c=%.6f area %.3e ratio %.4f n=%d m=%d' % (p, c, a, a / (A_CARD * s * s), n, m))
