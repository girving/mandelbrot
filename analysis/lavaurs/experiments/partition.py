"""The island partition at level 2 (q = 2, gate -1): every census child (all own-half level-1 sources × targets,
deduplicated, satellites dropped by exact cusp) gets its canonical parent from `glavaurs parent` (pieces of the horn
map: Voronoi cells of its critical values).  Reports the mass with a parent among the level-1 components, with a
parent outside the list, orphans (Misiurewicz sheets), failures, and how often the census's source is the parent.

  python3 partition.py single.txt [--targets 20]"""
import argparse, math, os, subprocess
from collections import defaultdict
B = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../build/release/glavaurs')
K = 2.4674011
ap = argparse.ArgumentParser(); ap.add_argument('single'); ap.add_argument('--targets', type=int, default=20)
ap.add_argument('--threads', type=int, default=2); ap.add_argument('--save', default='')
a = ap.parse_args()
env = dict(os.environ, GL_SIDE='-1', CHILDREN_STARTS='26')
comps = []
for l in open(a.single):
    f = l.split(); sat = f[0].startswith('sat'); n = int(f[0][3:]) if sat else int(f[0])
    comps.append(('B%d' % len(comps) if sat else 'F%d' % len(comps), n, complex(float(f[1]), float(f[2])), float(f[3]), float(f[5]), sat))
cut = (-0.16135210336424927 - 4.1583377953217475) / 2
own = [c for c in comps if c[2].imag > cut]
targets = [c for c in own if c[5]] + sorted([c for c in own if not c[5]], key=lambda c: -c[3])[:a.targets]
Ct = {t[0]: t[4] for t in targets}; nt = {t[0]: t[1] for t in targets}
lines = []
for nm, n, s, C, _, _ in own:
    jc = -(n + 1)
    for tn, ntt, ct, _, _, _ in targets:
        lines.append('%s~%s 2 %d %.17g %.17g %.6g %d %.17g %.17g %d %d' % (nm, tn, n, s.real, s.imag, math.sqrt(C / (K * math.pi)), ntt, ct.real, ct.imag, jc - 16, jc + 4))
out = subprocess.run([B, '1', '2', 'children', str(a.threads)], input='\n'.join(lines) + '\n', capture_output=True, text=True, env=env).stdout
kids = {}
for l in out.splitlines():
    f = l.split(); nm, br, j = f[0].rsplit('|', 2); src, t = nm.rsplit('~', 1)
    c = complex(float(f[1]), float(f[2])); w = float(f[3])
    if int(f[4]) or w >= 1 or c.imag <= cut: continue
    k = (round(c.real % 1, 7), round(c.imag, 7))
    if k not in kids: kids[k] = [Ct[t] * w, c, nt[t] - int(j), {src}]
    else: kids[k][3].add(src)
heavy = [(k, v) for k, v in kids.items() if v[0] > 1e-10]
res = subprocess.run([B, '1', '2', 'area', str(a.threads)], input=''.join('H%d 2 %d %.17g %.17g\n' % (i, v[2], v[1].real, v[1].imag) for i, (k, v) in enumerate(heavy)),
                     capture_output=True, text=True, env=env).stdout
for l in res.splitlines():
    f = l.split()
    if f[3] == 'failed': continue
    k, v = heavy[int(f[0][1:])]
    if float(f[8]) >= 0.01: kids.pop(k, None)
    else: v[0] = float(f[6])
keys = list(kids)
pr = subprocess.run([B, '1', '2', 'parent', str(a.threads)], input=''.join('K%d 2 %.17g %.17g\n' % (i, kids[k][1].real, kids[k][1].imag) for i, k in enumerate(keys)),
                    capture_output=True, text=True, env=env).stdout
lev1 = {(round(c[2].real % 1, 6), round(c[2].imag, 6)): c for c in comps}
tot = sum(v[0] for v in kids.values())
by = defaultdict(float); cnt = defaultdict(int); src_ok = 0.0; per_parent = defaultdict(float); orphans = []
for l in pr.splitlines():
    f = l.split(); k = keys[int(f[0][1:])]; v = kids[k]; kind = f[1]
    if kind == 'parent':
        ps = complex(float(f[2]), float(f[3])); pk = (round(ps.real % 1, 6), round(ps.imag, 6))
        if pk in lev1:
            kind = 'parent (listed)'; per_parent[lev1[pk][0]] += v[0]
            if lev1[pk][0] in v[3]: src_ok += v[0]
        else: kind = 'parent (not listed)'
    elif kind == 'orphan': orphans.append((v[0], v[1], complex(float(f[2]), float(f[3]))))
    by[kind] += v[0]; cnt[kind] += 1
print('level-2 primitive census children: %d, mass %.6e' % (len(kids), tot))
for kd in sorted(by): print('  %-20s %7d components, mass %.6e (%.4f%%)' % (kd, cnt[kd], by[kd], 100 * by[kd] / tot))
print('  the parent is among the census sources that found it: mass %.6e' % src_ok)
for o in sorted(orphans, key=lambda o: -o[0])[:6]: print('  orphan C %.3e at %.6f%+.6fi (sheet end %.6f%+.6fi)' % (o[0], o[1].real, o[1].imag, o[2].real, o[2].imag))
top = sorted(per_parent.items(), key=lambda kv: -kv[1])[:6]
print('  heaviest islands: ' + ', '.join('%s %.3e' % kv for kv in top))
if a.save:
    with open(a.save, 'w') as fo:
        for l in pr.splitlines():
            f = l.split(); k = keys[int(f[0][1:])]; v = kids[k]
            fo.write('%.17g %.17g %d %.6e %s %s %s\n' % (v[1].real, v[1].imag, v[2], v[0], f[1], f[2], f[3]))
