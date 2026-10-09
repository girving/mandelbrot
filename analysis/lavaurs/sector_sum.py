"""An r-transit sector of the 1/2-root NRP layer from Lavaurs-model components (r ≥ 3; s1_sum.py is r = 2).

Inputs: lavaurs_area --island outputs with LAVAURS_R = r ("name|side|j r n cre cim area_hi area_lo C conv uaddr caddr
side own|other isl_re isl_im cusp"), and lavaurs_area outputs with LAVAURS_TUNE (lines "U*X r n cre cim area_hi area_lo
C conv cusp ratio").  A component counts if it is primitive (cusp ≈ 0) and not a tuning U*X (U a component with
fewer transits, X primitive in M): tunings are matched by center.  Centers are deduplicated mod 1 (the σ cylinder).

  python3 sector_sum.py r island[,…] tune[,…]

Prints the sector sum, the tuned and satellite masses removed, and the convergence of the island census in the
source and target ranks."""
import sys, gzip, math
from collections import defaultdict

def opened(p): return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)

def key(c): return (round(c.real % 1.0, 8) % 1.0, round(c.imag, 8))

if __name__ == '__main__':
    r = int(sys.argv[1])
    tuned = {}  # key -> label
    for p in sys.argv[3].split(','):
        for l in opened(p):
            f = l.split()
            if '*' not in f[0] or f[3] == 'failed' or int(f[1]) != r: continue
            tuned[key(complex(float(f[3]), float(f[4])))] = f[0]
    comp = {}  # key -> [C, center, cusp, names]
    failed = 0
    for p in sys.argv[2].split(','):
        for l in opened(p):
            f = l.split()
            if len(f) < 10: failed += 1; continue
            assert int(f[1]) == r, l
            c = complex(float(f[3]), float(f[4]))
            e = comp.setdefault(key(c), [float(f[7]), c, float(f[-1]), []])
            e[3].append(f[0])
    kinds = defaultdict(lambda: [0, 0.0])
    for k, v in comp.items():
        kind = 'satellite' if v[2] >= 1e-8 else 'tuned' if k in tuned else 'NRP'
        v.append(kind)
        kinds[kind][0] += 1; kinds[kind][1] += v[0]
    print('r = %d: %d distinct components mod 1 (%d failed lines); tunings known: %d' % (r, len(comp), failed, len(tuned)))
    for kind in ('NRP', 'tuned', 'satellite'):
        print('   %-9s %6d  mass %.12e' % (kind, kinds[kind][0], kinds[kind][1]))
    for kind in ('NRP', 'tuned', 'satellite'):
        print('largest %s:' % kind)
        for v in sorted((v for v in comp.values() if v[4] == kind), key=lambda v: -v[0])[:8]:
            lab = tuned.get(key(v[1]), '')
            print('   C %.6e cusp %.1e σ %.10f%+.10fi  %s %s' % (v[0], v[2], v[1].real, v[1].imag, lab, ' '.join(v[3][:2])))
    # Census convergence: NRP mass by the rank of the best (source, target) pair that found it
    rank_src, rank_tgt = {}, {}
    for v in comp.values():
        for nm in v[3]:
            st = nm.split('|')[0]
            s, t = st.split('~')
            rank_src.setdefault(s, len(rank_src)); rank_tgt.setdefault(t, len(rank_tgt))
    for cap in (10, 30, 100, 300, 10**9):
        m = sum(v[0] for v in comp.values() if v[4] == 'NRP' and any(
            rank_src[nm.split('|')[0].split('~')[0]] < cap and rank_tgt[nm.split('|')[0].split('~')[1]] < cap for nm in v[3]))
        print('   NRP mass from source and target rank < %s (input order): %.10e' % (cap if cap < 10**9 else '∞', m))
