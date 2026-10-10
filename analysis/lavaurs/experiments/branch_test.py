"""Does the δ-series of a child converge with ratio |σ_W - σ_X| / D (D = the frozen chain's closest approach to the
critical value)?  Levels 2 and 3 from a few heavy sources at q = 2 (gate -1)."""
import math, os, subprocess, sys
import numpy as np
B = '/Users/irving/mandelbrot/build/release/glavaurs'
env = dict(os.environ, GL_SIDE='-1')
K = 2.4674011
def run(args, text, starts=26):
    e = dict(env, CHILDREN_STARTS=str(starts))
    return subprocess.run([B, '1', '2'] + args, input=text, capture_output=True, text=True, env=e).stdout
single = sys.argv[1]; out = sys.argv[2]
comps = []
for l in open(single):
    f = l.split(); sat = f[0].startswith('sat'); n = int(f[0][3:]) if sat else int(f[0])
    comps.append(('B%d' % len(comps) if sat else 'F%d' % len(comps), n, complex(float(f[1]), float(f[2])), float(f[3]), float(f[5]), sat))
cut = (-0.16135210336424927 - 4.1583377953217475) / 2
own = [c for c in comps if c[2].imag > cut]
targets = [c for c in own if c[5]] + sorted([c for c in own if not c[5]], key=lambda c: -c[3])[:60]
Ct = {t[0]: t[4] for t in targets}; nt = {t[0]: t[1] for t in targets}
rad = lambda C: math.sqrt(C / (K * math.pi))
def children(srcs, r):
    lines = []
    for nm, n, s, C in srcs:
        jc = -(n + 1)
        for tn, ntt, ct, _, _, _ in targets:
            lines.append('%s~%s %d %d %.17g %.17g %.6g %d %.17g %.17g %d %d' % (nm, tn, r, n, s.real, s.imag, rad(C), ntt, ct.real, ct.imag, jc - 16, jc + 4))
    kids = {}
    for l in run(['children', '2'], '\n'.join(lines) + '\n').splitlines():
        f = l.split(); nm, br, j = f[0].rsplit('|', 2); src, t = nm.rsplit('~', 1)
        c = complex(float(f[1]), float(f[2])); w = float(f[3])
        if int(f[4]) or w >= 1 or c.imag <= cut: continue
        k = (round(c.real % 1, 8), round(c.imag, 8))
        if k not in kids: kids[k] = (Ct[t] * w, c, nt[t] - int(j), src, f[0])
    return list(kids.values())
src1 = [(c[0], c[1], c[2], c[3]) for c in sorted([c for c in own if not c[5]], key=lambda c: -c[3])[:6]]
lev2 = children(src1, 2)
print('level 2: %d children of %d sources' % (len(lev2), len(src1)), file=sys.stderr)
src2 = [('S%d' % i, n, s, C) for i, (C, s, n, _, _) in enumerate(sorted(lev2, key=lambda v: -v[0])[:8])]
lev3 = children(src2, 3)
print('level 3: %d children of %d sources' % (len(lev3), len(src2)), file=sys.stderr)
sx_of = {nm: s for nm, n, s, C in src2}
# the exact target point y = Θ_3(σ_W) (locate), then frefine
loc = run(['locate'], ''.join('W%d 3 %.17g %.17g\n' % (i, v[1].real, v[1].imag) for i, v in enumerate(lev3)))
ys = {}
for l in loc.splitlines():
    f = l.split()
    if len(f) > 3 and f[2] != 'failed': ys[f[0]] = complex(float(f[2]), float(f[3]))
text = ''
for i, v in enumerate(lev3):
    nm = 'W%d' % i
    if nm not in ys: continue
    y = ys[nm]; sx = sx_of[v[3]]
    text += '%s 3 %.17g %.17g %.17g %.17g %.17g %.17g\n' % (nm, y.real, y.imag, sx.real, sx.imag, v[1].real, v[1].imag)
res = run(['frefine', '2'], text)
rows = []
for l in res.splitlines():
    f = l.split()
    if f[1] == 'failed': continue
    i = int(f[0][1:]); C, sw, n, src, _ = lev3[i]; sx = sx_of[src]
    sf = complex(float(f[1]), float(f[2])); s1 = complex(float(f[6]), float(f[7])); s2 = complex(float(f[8]), float(f[9])); D = float(f[10])
    d = abs(sw - sx)
    y = ys['W%d' % i]; tau = 2j * math.pi * 1.375 / 2
    Dt = min(abs((lambda z: z - round(z.real * 2) / 2)(y - sx - k * tau)) for k in range(-2, 3))
    rows.append((C, d, D, abs(sf - sw) / d, abs(s1 - sw) / d, abs(s2 - sw) / d, Dt, abs(sx - sw)))
np.save(out, np.array(rows))
print('%d children refined' % len(rows), file=sys.stderr)
