"""Residual of the frozen-σ operator against |ε| at level 3 of a q = 2 gl_census run (glavaurs frefine): every
census-kept child of the top sources refined under σ frozen at its source; per source the median relative position
error and |log w| against ε_X = 1/|Θ'_2(σ_X)|.

  python3 gl_eps_scaling.py dir nsources   (dir holds q2m.txt and q2u26_r2/_r3/_raw_r3/_tuned_r2/_tuned_r3.txt)"""
import subprocess, os, sys, math
from collections import defaultdict
S, BIN = sys.argv[1], os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../build/release/glavaurs')
p, q, gate, single, prefix, r, nsrc = 1, 2, -1, S + '/q2m.txt', S + '/q2u26', 3, int(sys.argv[2])
env = dict(os.environ, GL_SIDE=str(gate))
T = {}; nf = ns = 0
for l in open(single):
    f = l.split()
    if f[0].startswith('sat'): T['B%d' % ns] = (complex(float(f[1]), float(f[2])), float(f[3])); ns += 1
    else: T['F%d' % nf] = (complex(float(f[1]), float(f[2])), float(f[3])); nf += 1
src = {}
src_C = {}
for l in open('%s_r%d.txt' % (prefix, r - 1)):
    f = l.split(); src[f[0]] = (complex(float(f[3]), float(f[4])), float(f[5]), f[6] == '1'); src_C[f[0]] = float(f[5])
tuned = {l.split()[0] for l in open('%s_tuned_r%d.txt' % (prefix, r - 1))}
top = sorted([k for k, v in src.items() if not v[2] and k not in tuned], key=lambda k: -src[k][1])[:nsrc]
tops = set(top)
final, tun3 = {}, {l.split()[0] for l in open('%s_tuned_r%d.txt' % (prefix, r))}
for l in open('%s_r%d.txt' % (prefix, r)):
    f = l.split(); c = complex(float(f[3]), float(f[4]))
    final[(round(c.real % 1, 7) % 1, round(c.imag, 7))] = (f[6] == '1' or f[0] in tun3)
lines, meta, seen = [], [], set()
for l in open('%s_raw_r%d.txt' % (prefix, r)):
    f = l.split(); k = f[0].rsplit('|', 2); s_, t = k[0].split('~'); j = int(k[2])
    if s_ not in tops or int(f[4]) or float(f[3]) >= 1: continue
    c = complex(float(f[1]), float(f[2])); kk = (round(c.real % 1, 7) % 1, round(c.imag, 7))
    if final.get(kk, True) or (s_, kk) in seen: continue
    seen.add((s_, kk))
    y = T[t][0] + j / q; sx = src[s_][0]
    lines.append('c%d %d %.17g %.17g %.17g %.17g %.17g %.17g' % (len(meta), r, y.real, y.imag, sx.real, sx.imag, c.real, c.imag))
    meta.append((s_, T[t][1], c))
out = subprocess.run([BIN, str(p), str(q), 'frefine', '8'], input='\n'.join(lines) + '\n', capture_output=True, text=True, env=env).stdout
loc = subprocess.run([BIN, str(p), str(q), 'locate'], input=''.join('%s %d %.17g %.17g\n' % (k, r - 1, src[k][0].real, src[k][0].imag) for k in top), capture_output=True, text=True, env=env).stdout
eps = {l.split()[0]: 1 / abs(complex(float(l.split()[4]), float(l.split()[5]))) for l in loc.splitlines()}
per = defaultdict(list)
for l in out.splitlines():
    f = l.split(); i = int(f[0][1:]); s_, Ct, c = meta[i]
    sx = src[s_][0]
    if f[1] == 'failed': continue
    sf = complex(float(f[1]), float(f[2])); we, wT, wm = float(f[3]), float(f[4]), float(f[5])
    m = Ct * we
    if not (we > 0 and wT > 0 and wm > 0): continue
    dpos = abs(sf - c) / abs(c - sx)
    per[s_].append((dpos, abs(math.log(we / wT)), abs(math.log(we / wm))))
fails = sum(1 for l in out.splitlines() if l.split()[1] == 'failed')
print('%d children of %d sources refined (%d failed)' % (len(lines), len(top), fails))
def qt(v, f): v = sorted(v); return v[min(len(v) - 1, int(f * len(v)))]
print('source   C_X       |eps_X|  children   dpos/dist median,90%    |log w/wT| median,90%   |log w/wmult| median')
rows = sorted(per, key=lambda k: eps[k])
xs, ys, zs = [], [], []
for k in rows:
    v = per[k]; d = [x[0] for x in v]; a = [x[1] for x in v]; b = [x[2] for x in v]
    print('%-9s %.3e %.4f %7d   %.2e %.2e       %.2e %.2e         %.2e' % (k, src_C[k], eps[k], len(v), qt(d, .5), qt(d, .9), qt(a, .5), qt(a, .9), qt(b, .5)))
    xs.append(math.log(eps[k])); ys.append(math.log(qt(a, .5))); zs.append(math.log(qt(d, .5)))
def slope(x, y):
    n = len(x); mx = sum(x) / n; my = sum(y) / n
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / sum((a - mx) ** 2 for a in x)
print('log-log slope of the medians vs |eps|: weight error (T kept) %.2f, position error %.2f' % (slope(xs, ys), slope(xs, zs)))
