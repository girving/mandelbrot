"""The two-transit sector S_1 from Lavaurs-model components: the union of the model runs on the M-side labels
(lavaurs_area on diag seeds) and the island census (lavaurs_area --island), deduplicated by center, without the
doublings (a single-transit family's period-doubling satellite is tuned; it is attached to the family, so its center
lies within 2.5 radii of the family's center).

  python3 s1_sum.py model_census model_diag[,…] model_island[,…]

Prints S_1, the masses by source (labels only, islands only, both), and how the island-found mass grows with the
source's and target's rank (convergence of the enumeration)."""
import sys, gzip, math
from collections import defaultdict

def opened(p): return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)

if __name__ == '__main__':
    fam = {}   # single-transit families: name -> (C, center, radius in σ)
    for l in opened(sys.argv[1]):
        r = l.split()
        if r[3] == 'failed': continue
        area = float(r[5]); fam[r[0]] = (float(r[7]), complex(float(r[3]), float(r[4])), math.sqrt(area / math.pi))
    rank = {nm: i for i, nm in enumerate(sorted(fam, key=lambda x: -fam[x][0]))}
    rank['bulb'] = -1
    comp = {}  # key -> [C, center, sources set, best (source rank, target rank)]
    def add(key, C, c, src, ranks=None):
        e = comp.setdefault(key, [C, c, set(), (10**9, 10**9)])
        e[2].add(src)
        if ranks and ranks < e[3]: e[3] = ranks
    for p in sys.argv[2].split(','):
        for l in opened(p):
            r = l.split()
            if r[3] == 'failed': continue
            c = complex(float(r[3]), float(r[4])); add((round(c.real, 7), round(c.imag, 7)), float(r[7]), c, 'label')
    doubling = set()
    for p in sys.argv[3].split(','):
        for l in opened(p):
            r = l.split(); name, side, j = r[0].split('|'); u, t = name.split('~')
            c = complex(float(r[3]), float(r[4])); key = (round(c.real, 7), round(c.imag, 7))
            add(key, float(r[7]), c, 'island', (max(rank.get(u, 10**9), 0), rank.get(t, 10**9)))
            # a doubling of the target family t: attached to t (found from any source region)
            if t in fam and abs(c - fam[t][1]) < 2.5 * fam[t][2]: doubling.add(key)
    nr = {k: v for k, v in comp.items() if k not in doubling}
    S1 = sum(v[0] for v in nr.values())
    by = defaultdict(float)
    for v in nr.values(): by['+'.join(sorted(v[2]))] += v[0]
    print('two-transit components: %d distinct, %d doublings removed (mass %.4e)' % (len(comp), len(doubling), sum(comp[k][0] for k in doubling)))
    print('S_1 = %.12e  (by source: %s)' % (S1, ', '.join('%s %.4e' % kv for kv in sorted(by.items()))))
    # convergence of the island census in the source and target ranks (components found from islands only)
    isl = [v for v in nr.values() if 'island' in v[2]]
    for cap in (10, 30, 100, 150, 1000, 10**9):
        m = sum(v[0] for v in isl if v[3][0] < cap and v[3][1] < cap)
        print('   island-found mass with source and target rank < %s: %.8e' % (cap if cap < 10**9 else '∞', m))
