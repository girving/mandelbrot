"""Is the island decomposition complete?  Level-2 primitive mass from island children (ichildren) of every own-half
level-1 source vs census children (children, Newton from wide starts), same targets, deduplicated (q = 2, gate -1).

  python3 island_complete.py single.txt [--targets 20]"""
import argparse, math, os, subprocess
B = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../build/release/glavaurs')
K = 2.4674011
ap = argparse.ArgumentParser(); ap.add_argument('single'); ap.add_argument('--targets', type=int, default=20)
ap.add_argument('--threads', type=int, default=2); a = ap.parse_args()
env = dict(os.environ, GL_SIDE='-1', CHILDREN_STARTS='26')
comps = []
for l in open(a.single):
    f = l.split(); sat = f[0].startswith('sat'); n = int(f[0][3:]) if sat else int(f[0])
    comps.append(('B%d' % len(comps) if sat else 'F%d' % len(comps), n, complex(float(f[1]), float(f[2])), float(f[3]), float(f[5]), sat))
cut = (-0.16135210336424927 - 4.1583377953217475) / 2
own = [c for c in comps if c[2].imag > cut]
targets = [c for c in own if c[5]] + sorted([c for c in own if not c[5]], key=lambda c: -c[3])[:a.targets]
Ct = {t[0]: t[4] for t in targets}; nt = {t[0]: t[1] for t in targets}
srcs = [c for c in own if not c[5]] + [c for c in own if c[5]]
def collect(mode):
    lines = []
    for nm, n, s, C, _, _ in srcs:
        jc = -(n + 1)
        for tn, ntt, ct, _, _, _ in targets:
            if mode == 'children':
                lines.append('%s~%s 2 %d %.17g %.17g %.6g %d %.17g %.17g %d %d' % (nm, tn, n, s.real, s.imag, math.sqrt(C / (K * math.pi)), ntt, ct.real, ct.imag, jc - 16, jc + 4))
            else:
                lines.append('%s~%s 2 %.17g %.17g %d %.17g %.17g %d %d' % (nm, tn, s.real, s.imag, ntt, ct.real, ct.imag, jc - 16, jc + 4))
    out = subprocess.run([B, '1', '2', mode, str(a.threads)], input='\n'.join(lines) + '\n', capture_output=True, text=True, env=env).stdout
    kids, srcs_of = {}, {}
    for l in out.splitlines():
        f = l.split(); nm, br, j = f[0].rsplit('|', 2); src, t = nm.rsplit('~', 1)
        c = complex(float(f[1]), float(f[2])); w = float(f[3])
        if (mode == 'children' and int(f[4])) or w >= 1 or c.imag <= cut: continue
        k = (round(c.real % 1, 7), round(c.imag, 7))
        kids.setdefault(k, (Ct[t] * w, c, nt[t] - int(j)))
        srcs_of.setdefault(k, set()).add(src)
    multi = [k for k, v in srcs_of.items() if len(v) > 1]
    print('%s: %d components found from more than one source (mass %.3e of %.3e)' % (mode, len(multi), sum(kids[k][0] for k in multi), sum(v[0] for v in kids.values())), flush=True)
    # exact areas; drop satellites
    heavy = [(k, v) for k, v in kids.items() if v[0] > 1e-10]
    res = subprocess.run([B, '1', '2', 'area', str(a.threads)], input=''.join('H%d 2 %d %.17g %.17g\n' % (i, v[2], v[1].real, v[1].imag) for i, (k, v) in enumerate(heavy)),
                         capture_output=True, text=True, env=env).stdout
    prim = dict(kids)
    for l in res.splitlines():
        f = l.split()
        if f[3] == 'failed': continue
        k, v = heavy[int(f[0][1:])]
        if float(f[8]) >= 0.01: prim.pop(k, None)
        else: prim[k] = (float(f[6]),) + v[1:]
    return prim
ic = collect('ichildren'); cc = collect('children')
both = set(ic) & set(cc)
m = lambda d, ks: sum(d[k][0] for k in ks)
print('island children: %d components, primitive mass %.6e' % (len(ic), m(ic, ic)))
print('census children: %d components, primitive mass %.6e' % (len(cc), m(cc, cc)))
print('in both: %d (mass %.6e); island only %d (%.3e); census only %d (%.3e)' % (len(both), m(cc, both), len(set(ic) - both), m(ic, set(ic) - both), len(set(cc) - both), m(cc, set(cc) - both)))
for k in sorted(set(cc) - both, key=lambda k: -cc[k][0])[:8]: print('  census only: C %.3e at %.6f%+.6fi' % (cc[k][0], cc[k][1].real, cc[k][1].imag))
for k in sorted(set(ic) - both, key=lambda k: -ic[k][0])[:5]: print('  island only: C %.3e at %.6f%+.6fi' % (ic[k][0], ic[k][1].real, ic[k][1].imag))
