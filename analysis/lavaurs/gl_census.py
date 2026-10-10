"""Multi-transit census at the p/q root in the general model (fast_census.py for any root), level by level.

An r-transit component is a preimage of a single-transit target t under Θ_r (Θ_r(σ) = σ_t + j/q, excursion n_t - j)
near an (r-1)-transit source U (Θ_r(σ_U) = σ_U - (n_U + 1)/q, so shifts are centred on j = -(n_U + 1)), with family
constant C ≈ C_t |Θ_r' Π H'|^-2.  Children come from `glavaurs children` (Newton from U's local quadratic plus circle
starts for heavy sources, each verified as a center), deduplicated by center mod 1 (g_{σ+1}^r = f^{rq} g_σ^r); heavy
children get exact areas and cusps (`glavaurs area`), satellites (cusp ≥ 0.01) are kept as sources but not counted.
One gate's σ-plane holds the single-transit components of both sides of the root (the two bulbs, one per side), but
only one half is the gate's own: its multi-transit children are the physical ones (at q = 2, 3 checked by S_1/S_0 =
0.155 and the bulb doubling's ratio 0.072; the other half is the other gate's single transits, whose children here are
spurious).  The halves are split at the mean Im σ of the two bulbs; the gate's half is the one whose bulb has excursion
n = q - 1 (σ ≈ -1).  Sources and targets come from that half only, and children are counted only there.

  python3 gl_census.py p q side single.txt --levels 3 [--sources N] [--targets N] [--out prefix]
single.txt: gl_single.py --out ("n re im C exact|est" rows, "sat<n> re im C cusp" rows)."""
import argparse, math, os, re, subprocess, sys

BIN = os.path.join(os.path.dirname(__file__), '../../build/release/glavaurs')

def key(c): return (round(c.real % 1.0, 8) % 1.0, round(c.imag, 8))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('p', type=int); ap.add_argument('q', type=int); ap.add_argument('side', type=int)
    ap.add_argument('single')
    ap.add_argument('--levels', type=int, default=2)
    ap.add_argument('--sources', type=int, default=300)
    ap.add_argument('--targets', type=int, default=100)
    ap.add_argument('--src-min', type=float, default=1e-12)
    ap.add_argument('--exact-min', type=float, default=1e-8)
    ap.add_argument('--jd', type=float, default=8)     # shift window below / above the centre, in units of 1 in Re σ
    ap.add_argument('--ju', type=float, default=2)
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--big-src', type=float, default=1e-7)
    ap.add_argument('--starts', type=int, default=8)
    ap.add_argument('--small-targets', type=int, default=30)
    ap.add_argument('--out', default='')
    a = ap.parse_args()
    q = a.q
    env = dict(os.environ, GL_SIDE=str(a.side))
    def run(args, text):
        return subprocess.run([BIN, str(a.p), str(q)] + args, input=text, capture_output=True, text=True, env=env,
                              timeout=36000).stdout
    # K_q: area_σ → C, and the side cut
    head = subprocess.run([BIN, str(a.p), str(q), 'consist'], capture_output=True, text=True, env=env).stderr
    a_K = float(re.search(r'K = ([0-9.e+-]+)', head).group(1))
    fams, sats = [], []
    for l in open(a.single):
        f = l.split()
        if f[0].startswith('sat'): sats.append(('B%d' % len(sats), int(f[0][3:]), complex(float(f[1]), float(f[2])), float(f[3])))
        else: fams.append(('F%d' % len(fams), int(f[0]), complex(float(f[1]), float(f[2])), float(f[3])))
    big2 = sorted(sats, key=lambda x: -x[3])[:2]
    cut = (big2[0][2].imag + big2[1][2].imag) / 2
    own = [b for b in big2 if b[1] == q - 1]
    assert len(own) == 1, big2
    upper = own[0][2].imag > cut
    inh = lambda c: (c.imag >= cut) == upper
    fams = [f for f in fams if inh(f[2])]; sats = [f for f in sats if inh(f[2])]
    print('%s half (cut Im σ = %.4f): %d single-transit families, S_0 = %.6e; satellites %s' % (
        'upper' if upper else 'lower', cut, len(fams), sum(f[3] for f in fams), ' '.join('%.6e' % f[3] for f in sats)))
    fams.sort(key=lambda x: -x[3])
    targets = sats + fams[:a.targets]
    radius = lambda C: math.sqrt(C / (a_K * math.pi))   # σ radius of a disk of constant C
    level = [(nm, 1, n, c, C, False) for nm, n, c, C in fams[:a.sources]] + [(nm, 1, n, c, C, True) for nm, n, c, C in sats]
    totals = {}
    jd, ju = int(round(a.jd * q)), int(round(a.ju * q))
    for r in range(2, a.levels + 1):
        srcs = [s for s in level if s[1] == r - 1]
        big, small = [], []
        for nm, rs, n, c, C, sat in srcs:
            jc = -(n + 1)
            heavy = sat or C > a.big_src
            for tn, nt, ct, Ct in (targets if heavy else targets[:len(sats) + a.small_targets]):
                (big if heavy else small).append('%s~%s %d %d %.17g %.17g %.6g %d %.17g %.17g %d %d' % (
                    nm, tn, r, n, c.real, c.imag, radius(C), nt, ct.real, ct.imag, jc - jd, jc + ju))
        print('level %d: %d sources (%d heavy or satellite) x %d targets' % (
            r, len(srcs), sum(1 for s_ in srcs if s_[5] or s_[4] > a.big_src), len(targets)), file=sys.stderr)
        env['CHILDREN_STARTS'] = str(a.starts)
        out = run(['children', str(a.threads)], '\n'.join(big) + '\n')
        env['CHILDREN_STARTS'] = '0'
        out += run(['children', str(a.threads)], '\n'.join(small) + '\n')
        if a.out:
            with open('%s_raw_r%d.txt' % (a.out, r), 'w') as fo: fo.write(out)
        Ct_of = {t[0]: t[3] for t in targets}
        nt_of = {t[0]: t[1] for t in targets}
        kids = {}
        for l in out.splitlines():
            f = l.split()
            nm, br, j = f[0].rsplit('|', 2); src, t = nm.rsplit('~', 1)
            c = complex(float(f[1]), float(f[2])); w = float(f[3])
            if int(f[4]) or w >= 1: continue
            k = key(c)
            if k not in kids: kids[k] = [Ct_of[t] * w, c, nt_of[t] - int(j), 'L%d_%d' % (r, len(kids)), False]
        heavy = [(k, v) for k, v in kids.items() if v[0] > a.exact_min]
        if heavy:
            text = ''.join('H%d %d %d %.17g %.17g\n' % (i, r, v[2], v[1].real, v[1].imag) for i, (k, v) in enumerate(heavy))
            for l in run(['area', str(a.threads)], text).splitlines():
                f = l.split()
                if f[3] == 'failed': continue
                k, v = heavy[int(f[0][1:])]
                v[0] = float(f[6])
                if float(f[8]) >= 0.01: v[4] = True
        out_half = sum(v[0] for v in kids.values() if not inh(v[1]))
        kids = {k: v for k, v in kids.items() if inh(v[1])}
        prim = {k: v for k, v in kids.items() if not v[4]}
        sat = {k: v for k, v in kids.items() if v[4]}
        tot = sum(v[0] for v in prim.values())
        totals[r] = tot
        print('level %d: %d children (%d exact), primitive mass %.6e; satellites %d (%s); outside the half %.3e' % (
            r, len(prim), len(heavy), tot, len(sat), ' '.join('%.4e' % v[0] for v in sorted(sat.values(), key=lambda v: -v[0])[:6]), out_half))
        top = sorted(prim.values(), key=lambda v: -v[0])[:6]
        print('   heaviest: ' + ', '.join('%.3e @%.4f%+.4fi' % (v[0], v[1].real, v[1].imag) for v in top))
        if a.out:
            with open('%s_r%d.txt' % (a.out, r), 'w') as fo:
                for k, v in kids.items(): fo.write('%s %d %d %.17g %.17g %.6e %d\n' % (v[3], r, v[2], v[1].real, v[1].imag, v[0], int(v[4])))
        sys.stdout.flush()
        level = [(v[3], r, v[2], v[1], v[0], False) for v in prim.values() if v[0] > a.src_min] + \
                [(v[3], r, v[2], v[1], v[0], True) for v in sat.values()]
    return totals

if __name__ == '__main__':
    main()
