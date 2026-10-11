"""Depth saturation of the model's cascade rank: rows = level-1 parents U, columns = (child label, grandchild label)
with child label (t1, j1 - j_c(U)) and grandchild label (t2, j2 - j_c(W)), sheets summed; entries = grandchild masses.
If the cascade is a linear map on a small context, the mass-weighted rank of the two-level matrix grows additively
over the one-level rank, not multiplicatively.

  python3 cascade_rank2.py single.txt ich_level1.txt [--parents 60] [--kids 12] [--targets 16]"""
import argparse, math, os, subprocess
import numpy as np
B = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../build/release/glavaurs')
ap = argparse.ArgumentParser(); ap.add_argument('single'); ap.add_argument('ich')
ap.add_argument('--parents', type=int, default=60); ap.add_argument('--kids', type=int, default=12); ap.add_argument('--targets', type=int, default=16)
ap.add_argument('--jlo', type=int, default=-14); ap.add_argument('--jhi', type=int, default=3); ap.add_argument('--threads', type=int, default=2)
ap.add_argument('--out', default=None); ap.add_argument('--load', action='store_true'); a = ap.parse_args()
env = dict(os.environ, GL_SIDE='-1')
comps = []
for l in open(a.single):
    f = l.split(); sat = f[0].startswith('sat'); n = int(f[0][3:]) if sat else int(f[0])
    comps.append(('B%d' % len(comps) if sat else 'F%d' % len(comps), n, complex(float(f[1]), float(f[2])), float(f[3]), float(f[5]), sat))
info = {c[0]: c for c in comps}
cut = (-0.16135210336424927 - 4.1583377953217475) / 2
own = [c for c in comps if c[2].imag > cut]
targets = [c for c in own if c[5]] + sorted([c for c in own if not c[5]], key=lambda c: -c[3])[:a.targets]
# level-1 island children: (parent, t1, j1rel, sheet) -> (σ, mass, n)
kids = {}
for l in open(a.ich):
    f = l.split(); nm, br, j = f[0].rsplit('|', 2); src, t = nm.rsplit('~', 1); w = float(f[3])
    if w >= 1: continue
    c = complex(float(f[1]), float(f[2]))
    if c.imag <= cut: continue
    kids.setdefault(src, []).append(((t, int(j) + info[src][1] + 1, br), c, info[t][4] * w, info[t][1] - int(j)))
parents = sorted(kids, key=lambda s: -info[s][3])[:a.parents]
lines = []
for p in parents:
    for lab, c, m, n in sorted(kids[p], key=lambda x: -x[2])[:a.kids]:
        jc = -(n + 1)
        for tn, ntt, ct, _, _, _ in targets:
            lines.append('%s#%s#%d#%s~%s 2 %.17g %.17g %d %.17g %.17g %d %d' % (p, lab[0], lab[1], lab[2], tn, c.real, c.imag, ntt, ct.real, ct.imag, jc + a.jlo, jc + a.jhi))
if a.load: out = open(a.out).read()
else:
    out = subprocess.run([B, '1', '2', 'ichildren', str(a.threads)], input='\n'.join(lines) + '\n', capture_output=True, text=True, env=env).stdout
    if a.out: open(a.out, 'w').write(out)
nkid = {}
for p in parents:
    for lab, c, m, n in kids[p]: nkid[(p, lab)] = (m, n)
one, two = {}, {}
for p in parents:
    for lab, c, m, n in kids[p]: one[(p, (lab[0], lab[1]))] = one.get((p, (lab[0], lab[1])), 0) + m
for l in out.splitlines():
    f = l.split(); nm, br, j = f[0].rsplit('|', 2); src, t2 = nm.rsplit('~', 1); w = float(f[3])
    if w >= 1: continue
    p, t1, j1, s1 = src.split('#'); mk, nk = nkid[(p, (t1, int(j1), s1))]
    key = (p, (t1, int(j1), t2, int(j) + nk + 1))
    two[key] = two.get(key, 0) + info[t2][4] * w      # grandchild mass relative to the child's C (C_t2 w)
    # absolute: child's mass × relative; w here is |Θ'ΠH'|^-2 for the level-2 island, already absolute in σ
def spectrum(d, name):
    rows = {p: i for i, p in enumerate(parents)}; cols = {}
    for (p, c) in d: cols.setdefault(c, len(cols))
    A = np.zeros((len(rows), len(cols)))
    for (p, c), v in d.items(): A[rows[p], cols[c]] += v
    s = np.linalg.svd(A / A.sum(), compute_uv=False); s /= s[0]
    rk = lambda e: int((s > e).sum())
    print('%s: %d x %d, mass-weighted rank at 1e-2/1e-4/1e-6/1e-8/1e-10: %d %d %d %d %d' % (name, A.shape[0], A.shape[1], rk(1e-2), rk(1e-4), rk(1e-6), rk(1e-8), rk(1e-10)))
    print('   ', ' '.join('%.1e' % x for x in s[:48:3]))
spectrum(one, 'one level (children)')
spectrum(two, 'two levels (grandchildren)')
