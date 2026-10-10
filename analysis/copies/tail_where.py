"""Where does the slowly converging (p^-3) part of the copy area live?  Real NRPs in (-2, -1.7): area by period
split by lingering depth n near β ≈ 2 (and the largest lingering anywhere: the longest run of near-repeats z_{i+k} ≈ z_i
for small k, any repelling cycle)."""
import sys, math, collections
from load import load, nrp
roots = load(sys.argv[1])
tab = collections.defaultdict(lambda: collections.defaultdict(float))
tab2 = collections.defaultdict(lambda: collections.defaultdict(float))
for r in roots:
    if not nrp(r) or abs(r.c.imag) > 1e-12 or not (-2 < r.c.real < -1.7): continue
    c = r.c.real; p = r.p
    z = [0.0]
    for i in range(p): z.append(z[-1] ** 2 + c)
    bc = (1 + math.sqrt(1 - 4 * c)) / 2
    k = 2; n = 0
    while k < p and abs(z[k] - bc) < 0.25: n += 1; k += 1
    tab[p][min(n, 6)] += r.area
    # Longest lingering near any cycle of period ≤ 4: max run of i with |z_{i+per} - z_i| < 0.05
    best = 0
    for per in (1, 2, 3, 4):
        run = 0
        for i in range(1, p - per):
            if abs(z[i + per] - z[i]) < 0.05: run += 1; best = max(best, run * per)
            else: run = 0
    tab2[p][min(best // 2, 6)] += r.area
print('area share by lingering depth near β=2 (n: 0..5, ≥6), per period:')
for p in sorted(tab):
    tot = sum(tab[p].values())
    print('  p=%2d total %.2e: %s' % (p, tot, ' '.join('%d:%.3f' % (n, tab[p][n] / tot) for n in range(7))))
print('area share by longest lingering near any cycle of period ≤ 4 (steps/2 capped at 6), per period:')
for p in sorted(tab2):
    tot = sum(tab2[p].values())
    print('  p=%2d: %s' % (p, ' '.join('%d:%.3f' % (n, tab2[p][n] / tot) for n in range(7))))
