"""Fast census of the σ strip's components, level by level, from centers alone.

An r-transit component is a preimage of a single-transit target t under Θ_r (Θ_r(σ) = σ_t + j/2) near an
(r-1)-transit source, and its family constant is C ≈ C_t |Θ_r' Π H'|^-2 (checked to ~0.2% median on 3000 three-
transit components; large targets need exact areas).  So a level needs only Newton solves (lavaurs_area --children,
which also verifies each child is a center), no boundary tracing.  Children of a source carry ~0.02 of its mass
(W_U ∝ C_U), so sources below a cutoff can be dropped with a controlled loss.

Per level: sources (components of the previous level above --src-min, plus given satellite sources) × targets (the
bulb and the --targets heaviest families) × shifts j in [-(n_U+1) - JD, -(n_U+1) + JU]; children deduplicated by
center mod 1; heavy children (C > --exact-min) recomputed exactly with lavaurs_area (area, cusp) — satellites
dropped by cusp.

  python3 fast_census.py --census model_census.txt.gz --levels 2 [--sources N] [--targets N] [--out prefix]"""
import argparse, gzip, math, os, subprocess, sys
from collections import defaultdict

BIN = os.path.join(os.path.dirname(__file__), '../../build/release/lavaurs_area')
BULB = ('bulb', 1, complex(-1.0074583370365449, 0.16135210336429348), 0.20812104826488634)

def key(c): return (round(c.real % 1.0, 8) % 1.0, round(c.imag, 8))

def radius(C): return math.sqrt(4 * C / math.pi ** 3)   # σ radius of a component with constant C (C = π²/4 area)

def run(args, text, threads, env=None):
    e = dict(os.environ); e.update(env or {})
    p = subprocess.run([BIN] + args + [str(threads)], input=text, capture_output=True, text=True, timeout=36000, env=e)
    return p.stdout

def clean(name): return name.replace('~', '^').replace('|', '.')

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--census', required=True)
    ap.add_argument('--levels', type=int, default=2)
    ap.add_argument('--sources', type=int, default=300)
    ap.add_argument('--targets', type=int, default=100)
    ap.add_argument('--src-min', type=float, default=1e-12)
    ap.add_argument('--exact-min', type=float, default=1e-7)
    ap.add_argument('--jd', type=int, default=16)
    ap.add_argument('--ju', type=int, default=4)
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--extra-sources', default='')   # file: "name r n re im C" (e.g. satellites B_m as sources)
    ap.add_argument('--out', default='')
    ap.add_argument('--big-src', type=float, default=1e-7)     # sources above this (and satellites): many starts
    ap.add_argument('--starts', type=int, default=8)
    ap.add_argument('--small-targets', type=int, default=30)
    ap.add_argument('--start', default='')   # file "name r n re im C sat": the starting level's components (sat = 1:
                                             # a satellite, a source but not counted)
    a = ap.parse_args()
    fams = []
    for l in gzip.open(a.census, 'rt'):
        f = l.split()
        if f[3] != 'failed': fams.append((f[0], int(f[2]), complex(float(f[3]), float(f[4])), float(f[7])))
    fams.sort(key=lambda x: -x[3])
    targets = [BULB] + fams[:a.targets]
    # level 1 sources: single-transit families and the bulb (a satellite: a source, not counted)
    level = [(nm, 1, n, c, C, False) for nm, n, c, C in fams[:a.sources]] + [('bulb', 1, 1, BULB[2], BULB[3], True)]
    extra = []
    if a.extra_sources:
        for l in open(a.extra_sources):
            f = l.split(); extra.append((f[0], int(f[1]), int(f[2]), complex(float(f[3]), float(f[4])), float(f[5]), True))
    totals = {}
    r0 = 2
    if a.start:
        level = []
        for l in open(a.start):
            f = l.split(); level.append((f[0], int(f[1]), int(f[2]), complex(float(f[3]), float(f[4])), float(f[5]), f[6] == '1'))
        r0 = level[0][1] + 1
    for r in range(r0, a.levels + 1):
        srcs = [s for s in level if s[1] == r - 1] + [s for s in extra if s[1] == r - 1]
        big, small = [], []
        for nm, rs, n, c, C, sat in srcs:
            jc = -(n + 1)
            heavy_src = sat or C > a.big_src
            for tn, nt, ct, Ct in (targets if heavy_src else targets[:a.small_targets + 1]):
                (big if heavy_src else small).append('%s~%s %d %d %.17g %.17g %.6g %d %.17g %.17g %d %d' % (
                    clean(nm), tn, r, n, c.real, c.imag, radius(C), nt, ct.real, ct.imag, jc - a.jd, jc + a.ju))
        print('level %d: %d sources (%d heavy or satellite) × %d targets' % (
            r, len(srcs), sum(1 for s_ in srcs if s_[5] or s_[4] > a.big_src), len(targets)), file=sys.stderr)
        out = run(['--children'], '\n'.join(big) + '\n', a.threads, {'CHILDREN_STARTS': str(a.starts)})
        out += run(['--children'], '\n'.join(small) + '\n', a.threads, {'CHILDREN_STARTS': '0'})
        if a.out:
            with open('%s_raw_r%d.txt' % (a.out, r), 'w') as fo: fo.write(out)
        Ct_of = {t[0]: t[3] for t in targets}
        nt_of = {t[0]: t[1] for t in targets}
        kids = {}
        for l in out.splitlines():
            f = l.split()
            nm, br, j = f[0].rsplit('|', 2); src, t = nm.rsplit('~', 1)
            c = complex(float(f[1]), float(f[2])); w = float(f[3]); satf = int(f[4])
            if satf or w >= 1: continue   # w ≥ 1: a degenerate solution (a single-transit center itself: H' = 0)
            k = key(c)
            C = Ct_of[t] * w
            if k not in kids: kids[k] = [C, c, nt_of[t] - int(j), clean(f[0]), False]
        # exact areas (and cusp) for the heavy children
        heavy = [(k, v) for k, v in kids.items() if v[0] > a.exact_min]
        if heavy:
            text = ''.join('H%d %d %d %.17g %.17g\n' % (i, r, v[2], v[1].real, v[1].imag) for i, (k, v) in enumerate(heavy))
            res = run([], text, a.threads)
            for l in res.splitlines():
                f = l.split()
                i = int(f[0][1:]); k, v = heavy[i]
                if f[3] == 'failed': continue
                v[0] = float(f[7])
                if float(f[9]) >= 1e-8: v[4] = True               # satellite: a source, not counted
        prim = {k: v for k, v in kids.items() if v[0] > 0 and not v[4]}
        sats = {k: v for k, v in kids.items() if v[4]}
        tot = sum(v[0] for v in prim.values())
        totals[r] = tot
        print('level %d: %d children (%d heavy, exact), primitive mass %.6e; satellites %d (mass %.3e)' % (
            r, len(prim), len(heavy), tot, len(sats), sum(v[0] for v in sats.values())))
        top = sorted(prim.values(), key=lambda v: -v[0])[:6]
        print('   heaviest: ' + ', '.join('%.3e @%.4f%+.4fi' % (v[0], v[1].real, v[1].imag) for v in top))
        if a.out:
            with open('%s_r%d.txt' % (a.out, r), 'w') as fo:
                for k, v in kids.items():
                    if v[0] > 0: fo.write('%s %d %d %.17g %.17g %.6e %d\n' % (v[3], r, v[2], v[1].real, v[1].imag, v[0], int(v[4])))
        sys.stdout.flush()
        level = [(v[3], r, v[2], v[1], v[0], False) for v in prim.values() if v[0] > a.src_min] + \
                [(v[3], r, v[2], v[1], v[0], True) for v in sats.values()]
    return totals

if __name__ == '__main__':
    main()
