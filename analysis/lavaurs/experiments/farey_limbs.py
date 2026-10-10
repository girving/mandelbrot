"""κ of the Farey limbs t = 1/m of the q = 2 strip (CF [0; 2, k, m], k → ∞) as m grows: around each Farey bulb B_m
(farey_edge.py), depth 0 = the level-m children of B_{m-1} (B_m's siblings), depth 1 = the level-(m+1) children of
B_m and of the heaviest depth-0 components; all over the satellites + the --targets heaviest single transits,
satellites (cusp ≥ 0.01) and tunings by the big bulb (block test) removed, each classified by kneading
(glavaurs classify, limb b/m' with b = (offset + 1)/2 mod m').  Prints per m: C(B_m), the limb-1/m mass at depth 0
and 1 over C(B_m), counts, the heaviest components, and how much lands in other limbs.

  python3 farey_limbs.py single.txt edge.txt [--targets 60] [--top 10] [--threads 2]"""
import argparse, math, os, subprocess, sys
from collections import defaultdict
from fractions import Fraction
B = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../build/release/glavaurs')
K = 2.4674011
ap = argparse.ArgumentParser(); ap.add_argument('single'); ap.add_argument('edge')
ap.add_argument('--targets', type=int, default=60); ap.add_argument('--top', type=int, default=10)
ap.add_argument('--threads', type=int, default=2); ap.add_argument('--mmin', type=int, default=2); ap.add_argument('--mmax', type=int, default=99)
ap.add_argument('--verbose', action='store_true'); ap.add_argument('--jd', type=int, default=16); ap.add_argument('--jdm', type=int, default=6)   # shift window below jc: jd + jdm m
a = ap.parse_args()
env = dict(os.environ, GL_SIDE='-1', CHILDREN_STARTS='26')
def run(args, text): return subprocess.run([B, '1', '2'] + args, input=text, capture_output=True, text=True, env=env).stdout
comps = []
for l in open(a.single):
    f = l.split(); sat = f[0].startswith('sat'); n = int(f[0][3:]) if sat else int(f[0])
    comps.append(('B%d' % len(comps) if sat else 'F%d' % len(comps), n, complex(float(f[1]), float(f[2])), float(f[3]), float(f[5]), sat))
cut = (-0.16135210336424927 - 4.1583377953217475) / 2
own = [c for c in comps if c[2].imag > cut]
targets = [c for c in own if c[5]] + sorted([c for c in own if not c[5]], key=lambda c: -c[3])[:a.targets]
Ct = {t[0]: t[4] for t in targets}; nt = {t[0]: t[1] for t in targets}
bigB = max((c for c in own if c[5]), key=lambda c: c[3])
rad = lambda C: math.sqrt(C / (K * math.pi))
edge = {1: (1, bigB[1], bigB[2], bigB[3])}
for l in open(a.edge):
    f = l.split(); edge[int(f[0])] = (int(f[2]), int(f[3]), complex(float(f[4]), float(f[5])), float(f[6]))
def children(srcs, r):
    lines = []
    jcs = {}
    for nm, n, s, C in srcs:
        jc = -(n + 1); jcs[nm] = jc
        for tn, ntt, ct, _, _, _ in targets:
            lines.append('%s~%s %d %d %.17g %.17g %.6g %d %.17g %.17g %d %d' % (nm, tn, r, n, s.real, s.imag, rad(C), ntt, ct.real, ct.imag, jc - a.jd - a.jdm * r, jc + 4))
    kids = {}
    for l in run(['children', str(a.threads)], '\n'.join(lines) + '\n').splitlines():
        f = l.split(); nm, br, j = f[0].rsplit('|', 2); src, t = nm.rsplit('~', 1)
        c = complex(float(f[1]), float(f[2])); w = float(f[3])
        if int(f[4]) or w >= 1 or c.imag <= cut: continue
        k = (round(c.real % 1, 8), round(c.imag, 8))
        if k not in kids: kids[k] = ['W%d_%d' % (r, len(kids)), nt[t] - int(j), c, Ct[t] * w, int(j) - jcs[src], t]
    kids = list(kids.values())
    heavy = [k for k in kids if k[3] > 1e-12]
    drop = set()
    for l in run(['area', str(a.threads)], ''.join('%s %d %d %.17g %.17g\n' % (k[0], r, k[1], k[2].real, k[2].imag) for k in heavy)).splitlines():
        f = l.split()
        if f[3] == 'failed': continue
        k = next(k for k in heavy if k[0] == f[0])
        k[3] = float(f[6])
        if float(f[8]) >= 0.01: drop.add(f[0])
    kids = [k for k in kids if k[0] not in drop]
    # tunings by the big bulb (r = 1, n_B): W ~ (σ + j, n - 2 r j) in B's representative, n + 1 = 2 r ... (block test)
    nB, sB = bigB[1], bigB[2]
    pairs = []
    for nm, n, s, C, _, _ in kids:
        j = round((sB - s).real); n2 = n - j * r * 2
        if n2 + 1 == r * (nB + 1): pairs.append('%s %d %d %.17g %.17g B 1 %d %.17g %.17g' % (nm, r, n2, s.real + j, s.imag, nB, sB.real, sB.imag))
    tuned = set()
    for l in run(['blocktune', '0.6666666666666666', str(a.threads)], '\n'.join(pairs) + '\n').splitlines():
        f = l.split()
        if len(f) > 3 and f[2] == '0' and int(f[3]) > 0: tuned.add(f[0])
    return [tuple(k) for k in kids if k[0] not in tuned], len(tuned)
def limb_of(r, ks):
    out = {}
    res = run(['classify', '0.6666666666666666', str(a.threads)], ''.join('%s %d %d %.17g %.17g\n' % (k[0], r, k[1], k[2].real, k[2].imag) for k in ks))
    for l in res.splitlines():
        g = l.split()
        if len(g) == 3 and g[1] not in ('bulb', 'pre', 'failed'):
            mm, off = int(g[1]), int(g[2])
            out[g[0]] = Fraction(((off + 1) // 2) % mm, mm) if (off + 1) % 2 == 0 and mm > 0 else None
        else: out[g[0]] = g[1]
    return out
print('m  C(B_m)       m^4C    depth-0 limb/C  n0   depth-1 limb/C  n1   tuned  other-limb mass/C (top)  heaviest/C', flush=True)
for m in range(max(2, a.mmin), min(max(edge), a.mmax) + 1):
    rB, nBm, sBm, CBm = edge[m]
    prev = edge[m - 1]
    k0, t0 = children([('P', prev[1], prev[2], prev[3])], m)
    L0 = limb_of(m, k0)
    t_m = Fraction(1, m)
    d0 = [k for k in k0 if L0.get(k[0]) == t_m]
    other = defaultdict(float)
    for k in k0:
        if L0.get(k[0]) != t_m: other[str(L0.get(k[0]))] += k[3]
    srcs = [('Bm', nBm, sBm, CBm)] + [(k[0], k[1], k[2], k[3]) for k in sorted(d0, key=lambda k: -k[3])[:a.top]]
    h0 = max(d0, key=lambda k: k[3], default=None)
    k1, t1 = children(srcs, m + 1)
    L1 = limb_of(m + 1, k1)
    d1 = [k for k in k1 if L1.get(k[0]) == t_m]
    for k in k1:
        if L1.get(k[0]) != t_m: other[str(L1.get(k[0]))] += k[3]
    top_other = sorted(other.items(), key=lambda kv: -kv[1])[:2]
    if a.verbose:
        raw = run(['classify', '0.6666666666666666', str(a.threads)], ''.join('%s %d %d %.17g %.17g\n' % (k[0], m + 1, k[1], k[2].real, k[2].imag) for k in k1))
        rawd = {l.split()[0]: ' '.join(l.split()[1:]) for l in raw.splitlines()}
        for k in sorted(k1, key=lambda k: -k[3])[:8]:
            print('    depth-1 %s n %d C/C_Bm %.3e σ %.6f%+.6fi j-jc %d target %s classify: %s -> %s' % (k[0], k[1], k[3] / CBm, k[2].real, k[2].imag, k[4], k[5], rawd.get(k[0]), L1.get(k[0])))
    hv = max((k[3] for k in d0 + d1), default=0)
    h1 = max(d1, key=lambda k: k[3], default=None)
    hj = ' heaviest j-jc: d0 %s (%s) d1 %s (%s)' % (h0[4] if h0 else '-', h0[5] if h0 else '', h1[4] if h1 else '-', h1[5] if h1 else '')
    print('%2d %.4e %.4f  %.4e  %5d  %.4e  %5d  %3d+%-3d %s  %.3e' % (m, CBm, m ** 4 * CBm, sum(k[3] for k in d0) / CBm, len(d0), sum(k[3] for k in d1) / CBm, len(d1),
          t0, t1, ' '.join('%s:%.2e' % (t, v / CBm) for t, v in top_other), hv / CBm) + hj, flush=True)
