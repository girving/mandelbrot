"""Families near -2: copies (n, w0) with the same exit itinerary, as n grows.  G(w0) = lim 256^n area(n, w0);
H(w0) = G(w0) |(f_a^m)'(w0*)|^4 with w0* the precritical point of f_a = z² - 2 the family converges to.  If H is a
smooth function of w0*, Σ_w0 is a Ruelle operator L g(z) = Σ_{f(w)=z} |f'(w)|^-4 g(w) with analytic data."""
import sys, math, collections
from load import load, nrp, A_CARD
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
    # Exit itinerary: signs of the orbit after leaving β, up to 0 (the family key)
    itin = ''.join('R' if z[j] > 0 else 'L' for j in range(k, p))
    fam[itin][n] = (r.area, z[k])
fa = lambda w: w * w - 2
rows = []
for itin, d in fam.items():
    ns = sorted(d)
    if len(ns) < 3 or ns[-1] < 4: continue
    n = ns[-1]; a, w = d[n]
    m = len(itin)
    # Precritical point of f_a with this itinerary: invert f_a along the itinerary from 0 (w_m = 0)
    ws = 0.0  # z_p = 0; invert f_a along the whole exit itinerary z_k, ..., z_{p-1}
    for ch in reversed(itin):
        ws = math.sqrt(ws + 2) * (1 if ch == 'R' else -1)
    der = 1.0; x = ws
    for j in range(m): der *= 2 * x; x = fa(x)
    G = a * 256.0 ** n
    # convergence of 256^n area along the family
    conv = [d[k][0] * 256.0 ** k for k in ns]
    rows.append((ws, m, G, G * abs(der) ** 4, n, conv[-1] / conv[-2] - 1 if len(conv) > 1 else 0))
rows.sort()
import numpy as np
# Smoothness: fit log H as a polynomial in w0* over families converged to < 1e-3, and report the residual
good = [(ws, H) for ws, m, G, H, n, ch in rows if abs(ch) < 1e-3]
xs = np.array([g[0] for g in good]); ys = np.log(np.array([g[1] for g in good]))
dup = collections.Counter(np.round(xs, 6))
for deg in (2, 4, 6):
    cf = np.polyfit(xs, ys, deg); res = ys - np.polyval(cf, xs)
    print('log H(w0*) polynomial fit deg %d over %d converged families: rms %.2e, max %.2e' % (deg, len(xs), np.sqrt((res ** 2).mean()), abs(res).max()))
print('w0* values shared by two families: %d' % sum(1 for v in dup.values() if v > 1))
print('%d families with ≥3 members reaching n ≥ 4' % len(rows))
print('  w0*        m   G = 256^n area   H = G |(f_a^m)\'|^4   (last n, last step change of 256^n area)')
for ws, m, G, H, n, ch in rows:
    print('  %+.6f  %2d  %.6e      %.6e          (n=%d, %+.1e)' % (ws, m, G, H, n, ch))
