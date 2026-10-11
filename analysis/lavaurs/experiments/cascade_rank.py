"""Is the model's cascade low-rank?  Rows: parents U (level-1 islands), columns: child labels (target t, shift j
relative to the parent's j_c = -(n_U + 1)); entries: island-children weights w_+ + w_- (ichildren, both sheets, so the
sheet orientation convention does not matter), normalized per row.  Geometric singular values = a linear cascade with
a small context (as for the satellite tree); algebraic = the frozen kernel's singularity at the targets dominates.

  python3 cascade_rank.py single.txt [--sources 200] [--targets 12] [--jlo -10] [--jhi 2]"""
import argparse, math, os, subprocess
import numpy as np
B = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../build/release/glavaurs')
ap = argparse.ArgumentParser(); ap.add_argument('single'); ap.add_argument('--sources', type=int, default=200)
ap.add_argument('--targets', type=int, default=12); ap.add_argument('--jlo', type=int, default=-10); ap.add_argument('--jhi', type=int, default=2)
ap.add_argument('--threads', type=int, default=2); ap.add_argument('--out', default=None); a = ap.parse_args()
env = dict(os.environ, GL_SIDE='-1', ICH_DEBUG='1')
comps = []
for l in open(a.single):
    f = l.split(); sat = f[0].startswith('sat'); n = int(f[0][3:]) if sat else int(f[0])
    comps.append(('B%d' % len(comps) if sat else 'F%d' % len(comps), n, complex(float(f[1]), float(f[2])), float(f[3]), float(f[5]), sat))
cut = (-0.16135210336424927 - 4.1583377953217475) / 2
own = [c for c in comps if c[2].imag > cut]
targets = [c for c in own if c[5]] + sorted([c for c in own if not c[5]], key=lambda c: -c[3])[:a.targets]
srcs = sorted([c for c in own if not c[5]], key=lambda c: -c[3])[:a.sources]
lines = []
for nm, n, s, C, _, _ in srcs:
    jc = -(n + 1)
    for tn, ntt, ct, _, _, _ in targets:
        lines.append('%s~%s 2 %.17g %.17g %d %.17g %.17g %d %d' % (nm, tn, s.real, s.imag, ntt, ct.real, ct.imag, jc + a.jlo, jc + a.jhi))
res = subprocess.run([B, '1', '2', 'ichildren', str(a.threads)], input='\n'.join(lines) + '\n', capture_output=True, text=True, env=env)
out = res.stdout
if a.out: open(a.out, 'w').write(out); open(a.out + '.err', 'w').write(res.stderr)
# island data per source: σ_c and |κ''| (ICH_DEBUG lines "name~t: σ_c re+imi (...) κ_c ... κ'' x")
import re
isl = {}
for l in res.stderr.splitlines():
    m = re.match(r"(\S+)~\S+: σ_c (\S+?)([+-][0-9.e+-]+)i .* κ'' (\S+)", l)
    if m: isl[m.group(1)] = (complex(float(m.group(2)), float(m.group(3))), float(m.group(4)))
Ct = {t[0]: t[4] for t in targets}; nU = {s[0]: s[1] for s in srcs}
cols = {}; rows = {s[0]: i for i, s in enumerate(srcs)}
M = {}; R = {}
for l in out.splitlines():
    f = l.split(); nm, br, j = f[0].rsplit('|', 2); src, t = nm.rsplit('~', 1)
    w = float(f[3])
    if w >= 1: continue
    jr = int(j) + nU[src] + 1
    key = (t, jr); cols.setdefault(key, len(cols))
    M[(rows[src], cols[key])] = M.get((rows[src], cols[key]), 0) + w * Ct[t]
    if src in isl:
        sc, k2 = isl[src]; sg = complex(float(f[1]), float(f[2]))
        R[(rows[src], cols[key])] = R.get((rows[src], cols[key]), 0) + 0.5 * w * (k2 * abs(sg - sc)) ** 2   # w |κ'_quad|^2, both sheets averaged
A = np.zeros((len(rows), len(cols)))
for (i, j), v in M.items(): A[i, j] = v
keep = A.sum(1) > 0; A = A[keep]
print('%d parents x %d child labels, fill %.2f' % (A.shape[0], A.shape[1], (A > 0).mean()))
An = A / np.linalg.norm(A, axis=1, keepdims=True)
s = np.linalg.svd(An, compute_uv=False)
print('row-normalized singular values / s0:', ' '.join('%.1e' % x for x in (s / s[0])[:40]))
# mass-weighted: rows by their total (what the area sum sees)
s2 = np.linalg.svd(A / A.sum(), compute_uv=False)
print('mass-weighted:', ' '.join('%.1e' % x for x in (s2 / s2[0])[:40]))

# the explicit kernel divided out: R = w |κ'_quadratic(σ_child)|^2 should be a column factor (the landing |ΠH'|^-2)
Rm = np.zeros((len(rows), len(cols)))
for (i, j), v in R.items(): Rm[i, j] = v
Rm = Rm[keep]
cf = (Rm > 0).mean(0) >= 0.9          # well-populated child labels, then the parents complete on them
Rm = Rm[:, cf]; ok = (Rm > 0).all(1); Rm = Rm[ok]
print('kernel-normalized block: %d parents x %d labels' % Rm.shape)
colf = np.median(Rm, axis=0)
Q = Rm / colf
print('kernel-normalized: %d full rows; |log(R / column median)|: median %.2e, 90%% %.2e, max %.2e' % (
    len(Q), np.median(abs(np.log(Q))), np.quantile(abs(np.log(Q)), 0.9), abs(np.log(Q)).max()))
s3 = np.linalg.svd(Q, compute_uv=False)
print('kernel-normalized singular values / s0:', ' '.join('%.1e' % x for x in (s3 / s3[0])[:40]))
