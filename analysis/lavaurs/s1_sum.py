"""The two-transit sector S_1 from Lavaurs-model components, classified by the cusp measure.

Input: lavaurs_area (default mode) run on every distinct two-transit center found so far (M-side labels and island
census), lines "name r n cre cim area_hi area_lo C conv cusp [address]"; island names are "source~target|side|j".
A component is primitive (a non-renormalizable parameter of M, given two transits) if its multiplier map has a cusp
at μ = 1, cusp = |σ'(1)|/|a_1| ≈ 0, and a satellite (a doubling, the bulb's children, a cardioid bulb) if cusp = O(1).
The σ-plane is a cylinder (period 1 in Re σ: the same M family at limb k and k + 1), so components are deduplicated
by center mod 1.

  python3 s1_sum.py model_s1all.txt.gz [census.txt.gz]

Prints S_1, the satellite mass and kinds, the masses by provenance, and the cusp histogram."""
import sys, gzip, math
from collections import defaultdict

def opened(p): return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)

def key(c): return (round(c.real % 1.0, 8) % 1.0, round(c.imag, 8))

if __name__ == '__main__':
    comp = {}  # key -> [C, center, cusp, names, conv]
    failed = 0
    for l in opened(sys.argv[1]):
        r = l.split()
        if r[3] == 'failed': failed += 1; continue
        c = complex(float(r[3]), float(r[4]))
        e = comp.setdefault(key(c), [float(r[7]), c, float(r[9]), [], float(r[8])])
        e[3].append(r[0])
        if abs(e[0] - float(r[7])) > 1e-9 * e[0]: print('mismatch at %s: %s vs %s' % (key(c), e[0], r[7]))
    hist = defaultdict(lambda: [0, 0.0])
    for v in comp.values():
        b = math.floor(math.log10(max(v[2], 1e-30)))
        hist[b][0] += 1; hist[b][1] += v[0]
    print('%d distinct components mod 1 (%d failed)' % (len(comp), failed))
    print('cusp histogram (decade: count, mass):')
    for b in sorted(hist): print('   1e%+03d: %6d %.4e' % (b, hist[b][0], hist[b][1]))
    prim = {k: v for k, v in comp.items() if v[2] < 1e-8}
    sat = {k: v for k, v in comp.items() if v[2] >= 1e-8}
    S1 = sum(v[0] for v in prim.values())
    print('S_1 = %.13e over %d primitive components;  satellites %d, mass %.6e' % (
        S1, len(prim), len(sat), sum(v[0] for v in sat.values())))
    by = defaultdict(float)
    for v in prim.values():
        src = {'island' if '~' in nm else 'label' for nm in v[3]}
        by['+'.join(sorted(src))] += v[0]
    print('   by provenance: ' + ', '.join('%s %.6e' % kv for kv in sorted(by.items())))
    print('   worst conv among primitives: %.2e;  mass-weighted conv error %.2e' % (
        max(v[4] for v in prim.values()), sum(v[0] * v[4] for v in prim.values())))
    print('largest satellites:')
    for v in sorted(sat.values(), key=lambda v: -v[0])[:12]:
        print('   C %.6e cusp %.3f σ %.10f%+.10fi  %s' % (v[0], v[2], v[1].real, v[1].imag, ' '.join(v[3][:3])))
    print('largest primitives:')
    for v in sorted(prim.values(), key=lambda v: -v[0])[:12]:
        print('   C %.6e cusp %.1e σ %.10f%+.10fi  %s' % (v[0], v[2], v[1].real, v[1].imag, ' '.join(v[3][:3])))
