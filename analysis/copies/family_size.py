"""Families' areas against the small-copy size estimate: F = area / (A_card |s|²), s = 1/(β Λ²) from the center's
critical orbit (Λ = Π_{i<p} 2 z_i, β = Σ_{i<p} 1/Π_{j≤i} 2 z_j), and F's limit in k per family.

  python3 family_size.py census.txt[.gz] > F.txt     (lines "name F area"; limits and quantiles on stderr)"""
import sys, gzip, math
from collections import defaultdict
from family_wfa import neville
A_CARD = 3 * math.pi / 8
F = defaultdict(dict)
path = sys.argv[1]
for l in (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)):
    r = l.split()
    if r[3] == 'failed': continue
    c = complex(float(r[3]), float(r[4])); p = int(r[2]); area = float(r[5]) + float(r[10])
    z, prod, beta = 0j, 1 + 0j, 0j
    for i in range(1, p):
        z = z * z + c; prod *= 2 * z; beta += 1 / prod
    f = area / (A_CARD * abs(1 / (beta * prod * prod))**2)
    print('%s %.17g %.17g' % (r[0], f, area))
    name, k = r[0].rsplit('_', 1); F[name][int(k)] = f
lim = sorted(neville([1 / k for k in sorted(d)], [d[k] for k in sorted(d)]) for d in F.values())
q = lambda v, x: v[int(x * (len(v) - 1))]
print('F_inf quantiles (0, 1, 10, 50, 90, 99, 100%%): %s' % ' '.join('%.8f' % q(lim, x) for x in (0, .01, .1, .5, .9, .99, 1)),
      file=sys.stderr)
