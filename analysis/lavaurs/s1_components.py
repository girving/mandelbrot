"""Write the deduplicated, non-doubling two-transit components (as s1_sum.py) with their labels:
"key C n center_re center_im how" (how: label name or source~target|side|j), heaviest first.

  python3 s1_components.py model_census model_diag[,…] model_island[,…] > comps"""
import sys, gzip, math

def opened(p): return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)

fam = {}
for l in opened(sys.argv[1]):
    r = l.split()
    if r[3] != 'failed': fam[r[0]] = (complex(float(r[3]), float(r[4])), math.sqrt(float(r[5]) / math.pi))
fam['bulb'] = (complex(-1.0074583370365449, 0.16135210336429348), math.sqrt(0.084348283804329571 / math.pi))
comp, doubling = {}, set()
def add(key, C, n, c, how):
    if key not in comp: comp[key] = (C, n, c, how)
for p in sys.argv[2].split(','):
    for l in opened(p):
        r = l.split()
        if r[3] == 'failed': continue
        c = complex(float(r[3]), float(r[4])); add((round(c.real, 7), round(c.imag, 7)), float(r[7]), int(r[2]), c, r[0])
for p in sys.argv[3].split(','):
    for l in opened(p):
        r = l.split(); name, side, j = r[0].split('|'); u, t = name.split('~')
        c = complex(float(r[3]), float(r[4])); key = (round(c.real, 7), round(c.imag, 7))
        add(key, float(r[7]), int(r[2]), c, r[0])
        if t in fam and abs(c - fam[t][0]) < 2.5 * fam[t][1]: doubling.add(key)
for key, (C, n, c, how) in sorted(comp.items(), key=lambda x: -x[1][0]):
    if key in doubling: continue
    print('%s %.17g %d %.17g %.17g %s' % ('T%.7f%+.7f' % key, C, n, c.real, c.imag, how))
