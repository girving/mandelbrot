"""Hankel rank of the model's cascade on true island lineages (approach C).  Rows: islands at levels 1 and 2 (prefixes);
columns: child labels l = (target, shift relative to the island's j_c), sheets summed, and two-step suffixes (l1, l2)
with l1 in a fixed label set L1; entries: descendant masses (island children from ichildren, the near-field
tuning-like child at the island's critical value dropped).  A linear cascade with a K-dimensional state has Hankel
rank ≤ K whatever the depth.  Stages are cached in --dir.

  python3 cascade_hankel.py single.txt ich_level1.txt --dir D [--rows1 150] [--L1 8] [--rows2 200] [--targets 10]"""
import argparse, os, re, subprocess
import numpy as np
B = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../build/release/glavaurs')
ap = argparse.ArgumentParser(); ap.add_argument('single'); ap.add_argument('ich'); ap.add_argument('--dir', required=True)
ap.add_argument('--rows1', type=int, default=150); ap.add_argument('--L1', type=int, default=8); ap.add_argument('--rows2', type=int, default=200)
ap.add_argument('--targets', type=int, default=10); ap.add_argument('--jlo', type=int, default=-11); ap.add_argument('--jhi', type=int, default=2)
ap.add_argument('--threads', type=int, default=2); a = ap.parse_args()
os.makedirs(a.dir, exist_ok=True)
env = dict(os.environ, GL_SIDE='-1', ICH_DEBUG='1')
comps = []
for l in open(a.single):
    f = l.split(); sat = f[0].startswith('sat'); n = int(f[0][3:]) if sat else int(f[0])
    comps.append(('B%d' % len(comps) if sat else 'F%d' % len(comps), n, complex(float(f[1]), float(f[2])), float(f[3]), float(f[5]), sat))
info = {c[0]: c for c in comps}
cut = (-0.16135210336424927 - 4.1583377953217475) / 2
own = [c for c in comps if c[2].imag > cut]
targets = [c for c in own if c[5]] + sorted([c for c in own if not c[5]], key=lambda c: -c[3])[:a.targets]
tset = {t[0] for t in targets}
def isl_of(err):
    isl = {}
    for l in err.splitlines():
        m = re.match(r"(\S+)~\S+: σ_c (\S+?)([+-][0-9.e+-]+)i .* κ'' (\S+)", l)
        if m: isl[m.group(1)] = (complex(float(m.group(2)), float(m.group(3))), float(m.group(4)))
    return isl
def parse(out, err, nof):
    """children per source: list of (label (t, jrel), sheet, σ, mass, n); near field dropped"""
    isl = isl_of(err); kids = {}
    for l in out.splitlines():
        f = l.split(); nm, br, j = f[0].rsplit('|', 2); src, t = nm.rsplit('~', 1); w = float(f[3])
        if w >= 1 or src not in isl or t not in tset: continue
        c = complex(float(f[1]), float(f[2])); sc, k2 = isl[src]
        if c.imag <= cut or 0.5 * k2 * abs(c - sc) ** 2 < 0.03: continue
        kids.setdefault(src, []).append(((t, int(j) + nof(src) + 1), br, c, info[t][4] * w, info[t][1] - int(j)))
    return kids
def run(name, sources):   # sources: list of (name, r, σ, n); cached
    o = os.path.join(a.dir, name)
    if not os.path.exists(o):
        lines = []
        for nm, r, s, n in sources:
            jc = -(n + 1)
            for tn, ntt, ct, _, _, _ in targets:
                lines.append('%s~%s %d %.17g %.17g %d %.17g %.17g %d %d' % (nm, tn, r, s.real, s.imag, ntt, ct.real, ct.imag, jc + a.jlo, jc + a.jhi))
        res = subprocess.run([B, '1', '2', 'ichildren', str(a.threads)], input='\n'.join(lines) + '\n', capture_output=True, text=True, env=env)
        open(o, 'w').write(res.stdout); open(o + '.err', 'w').write(res.stderr)
    return open(o).read(), open(o + '.err').read()
# level 1: existing island children of the level-1 sources
k1 = parse(open(a.ich).read(), open(a.ich + '.err').read(), lambda s: info[s][1])
rows1 = sorted(k1, key=lambda s: -info[s][3])[:a.rows1]
tot = {}
for s in rows1:
    for lab, br, c, m, n in k1[s]: tot[lab] = tot.get(lab, 0) + m
L1 = sorted(tot, key=lambda l: -tot[l])[:a.L1]
print('L1 labels:', L1, flush=True)
# level 2 islands: children of rows1 with labels in L1 (both sheets)
src2 = []; n2 = {}
for s in rows1:
    for lab, br, c, m, n in k1[s]:
        if lab in L1:
            nm = '%s#%s#%d#%s' % (s, lab[0], lab[1], br); src2.append((nm, 3, c, n)); n2[nm] = (n, m)   # ichildren's r = the island map's transit count = level + 1
print('%d level-2 islands' % len(src2), flush=True)
o2, e2 = run('lev2.txt', src2)
k2 = parse(o2, e2, lambda s: n2[s][0])
# level 3 islands: children (labels in L1) of the heaviest rows2 level-2 islands
rows2 = sorted([s for s in n2 if s in k2], key=lambda s: -n2[s][1])[:a.rows2]
src3 = []; n3 = {}
for s in rows2:
    for lab, br, c, m, n in k2[s]:
        if lab in L1:
            nm = '%s@%s#%d#%s' % (s, lab[0], lab[1], br); src3.append((nm, 4, c, n)); n3[nm] = (n, m)
print('%d level-3 islands' % len(src3), flush=True)
o3, e3 = run('lev3.txt', src3)
k3 = parse(o3, e3, lambda s: n3[s][0])
# Hankel blocks.  Row X (an island); column l: Σ masses of X's children with label l; column (l1, l2): Σ over X's
# children Y with label l1 (both sheets) of Σ masses of Y's children with label l2.
def kids_of(X): return k1.get(X) or k2.get(X) or k3.get(X) or []
def name_of_child(X, lab, br): return ('%s#%s#%d#%s' if '#' not in X else '%s@%s#%d#%s') % (X, lab[0], lab[1], br)
rows = rows1 + rows2
cols1 = sorted({lab for X in rows for lab, *_ in kids_of(X)})
cols2 = sorted({(l1, lab) for X in rows for l1 in L1 for br in '+-' for lab, *_ in kids_of(name_of_child(X, l1, br))})
ci = {c: i for i, c in enumerate(cols1 + cols2)}
H = np.zeros((len(rows), len(ci)))
for i, X in enumerate(rows):
    for lab, br, c, m, n in kids_of(X):
        H[i, ci[lab]] += m
        for lab2, br2, c2, m2, n2_ in kids_of(name_of_child(X, lab, br)):
            if (lab, lab2) in ci: H[i, ci[(lab, lab2)]] += m2
print('Hankel %d rows (%d level 1, %d level 2) x %d columns (%d one-step, %d two-step)' % (len(rows), len(rows1), len(rows2), H.shape[1], len(cols1), len(cols2)))
def ranks(M, lab):
    s = np.linalg.svd(M, compute_uv=False); s = s / s[0]
    print('  %-44s rank at 1e-2/1e-4/1e-6/1e-8/1e-10: %s' % (lab, ' '.join(str(int((s > e).sum())) for e in (1e-2, 1e-4, 1e-6, 1e-8, 1e-10))))
# row normalization by the row's own one-step total (each island's children relative to itself), then weights
rn = H[:, :len(cols1)].sum(1, keepdims=True); ok = rn[:, 0] > 0
Hn = H[ok] / rn[ok]
n1 = sum(1 for X in rows1 if X in rows and ok[rows.index(X)])
ranks(Hn[:n1, :len(cols1)], 'level-1 rows, one-step')
ranks(Hn[:n1], 'level-1 rows, one- and two-step')
ranks(Hn[n1:], 'level-2 rows, one- and two-step')
ranks(Hn, 'all rows, all columns')
# mass weighting (what the area sum sees): absolute descendant masses, no per-row normalization
np.save(os.path.join(a.dir, 'H.npy'), H)
print('mass-weighted (absolute masses / block total):')
n1 = len(rows1); c1 = len(cols1)
ranks(H[:n1, :c1] / H[:n1, :c1].sum(), 'level-1 rows, one-step suffixes')
ranks(H[:n1, c1:] / H[:n1, c1:].sum(), 'level-1 rows, two-step suffixes')
ranks(H[n1:, :c1] / H[n1:, :c1].sum(), 'level-2 rows, one-step suffixes')
ranks(H[n1:, c1:] / H[n1:, c1:].sum(), 'level-2 rows, two-step suffixes')
B = np.vstack([H[:n1] / H[:n1].sum(), H[n1:] / H[n1:].sum()])
ranks(B, 'levels 1+2 rows (each block to unit mass), all suffixes')
# per-row normalized by the row's total mass, then rows weighted by sqrt(mass share) (a middle ground)
rt = H.sum(1, keepdims=True); okr = rt[:, 0] > 0
ranks(H[okr] / rt[okr] * np.sqrt(rt[okr] / rt[okr].sum()), 'all rows, row-normalized x sqrt(mass)')
# depth saturation of the row space with one-step suffixes: level-1, -2, -3 rows (each block to unit mass)
rows3 = sorted([s for s in n3 if s in k3], key=lambda s: -n3[s][1])[:300]
H3 = np.zeros((len(rows3), len(cols1))); c1i = {c: i for i, c in enumerate(cols1)}
for i, X in enumerate(rows3):
    for lab, br, c, m, n in k3[X]:
        if lab in c1i: H3[i, c1i[lab]] += m
B1, B2 = H[:n1, :c1] / H[:n1, :c1].sum(), H[n1:, :c1] / H[n1:, :c1].sum(); B3 = H3 / H3.sum()
print('one-step row spaces by depth (mass-weighted; trust ≥ 1e-6, weights printed to 7 digits):')
ranks(B3, 'level-3 rows (%d)' % len(rows3))
ranks(np.vstack([B1, B2]), 'levels 1+2 rows')
ranks(np.vstack([B1, B2, B3]), 'levels 1+2+3 rows')
ranks(np.vstack([B2, B3]), 'levels 2+3 rows')
