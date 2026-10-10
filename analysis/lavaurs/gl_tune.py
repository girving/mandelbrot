"""Tunings U*X among gl_census components, by the little-Julia-set test (glavaurs tuned): W (r transits) is in U's
copy of M (r_U | r, r_U < r) iff at σ_W the critical orbit of U's return map stays within U's little Julia set through
the r/r_U returns.  Candidate U: components of lower levels (primitive or satellite) whose center is within --reach
radii of σ_W (mod 1, W shifted into U's representative: (σ + j, n) ~ (σ, n + j r q)).

  python3 gl_tune.py p q gate single.txt prefix --level r [--reach 12]
reads single.txt (gl_single --out) and prefix_r{2..r}.txt (gl_census --out); prints the tuned components of level r
and writes prefix_tuned_r{r}.txt ("name U ratio")."""
import argparse, math, os, re, subprocess, sys
from collections import defaultdict

BIN = os.path.join(os.path.dirname(__file__), '../../build/release/glavaurs')

def load(a):
    comps = []
    for l in open(a.single):
        f = l.split()
        if f[0].startswith('sat'): comps.append(('B%d' % len(comps), 1, int(f[0][3:]), complex(float(f[1]), float(f[2])), float(f[3]), 1))
        else: comps.append(('F%d' % len(comps), 1, int(f[0]), complex(float(f[1]), float(f[2])), float(f[3]), 0))
    for r in range(2, a.level + 1):
        for l in open('%s_r%d.txt' % (a.prefix, r)):
            f = l.split(); comps.append((f[0], r, int(f[2]), complex(float(f[3]), float(f[4])), float(f[5]), int(f[6])))
    return comps

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('p', type=int); ap.add_argument('q', type=int); ap.add_argument('gate', type=int)
    ap.add_argument('single'); ap.add_argument('prefix')
    ap.add_argument('--level', type=int, required=True)
    ap.add_argument('--reach', type=float, default=12)
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--umin', type=float, default=1e-9)   # candidate tuners above this constant
    ap.add_argument('--wmin', type=float, default=0)
    a = ap.parse_args()
    env = dict(os.environ, GL_SIDE=str(a.gate))
    head = subprocess.run([BIN, str(a.p), str(a.q), 'consist'], capture_output=True, text=True, env=env).stderr
    K = float(re.search(r'K = ([0-9.e+-]+)', head).group(1))
    rad = lambda C: math.sqrt(C / (K * math.pi))
    comps = load(a)
    r = a.level
    W = [c for c in comps if c[1] == r and c[4] > a.wmin]
    U = [c for c in comps if c[1] < r and r % c[1] == 0 and c[4] > a.umin]
    # grid over Re σ mod 1 and Im σ
    G = 0.05
    grid = defaultdict(list)
    for u in U: grid[(int(math.floor((u[3].real % 1) / G)), int(math.floor(u[3].imag / G)))].append(u)
    lines = []
    for w in W:
        sw = w[3]
        x0, y0 = (sw.real % 1) / G, sw.imag / G
        # neighbours within reach radii: scan cells up to the largest tuner reach (bounded by 2 in σ)
        span = 2.0
        for dx in range(-int(span / G) - 1, int(span / G) + 2):
            for dy in range(-int(span / G) - 1, int(span / G) + 2):
                for u in grid.get(((int(math.floor(x0)) + dx) % int(1 / G), int(math.floor(y0)) + dy), []):
                    su = u[3]
                    j = round(su.real - sw.real)
                    s2 = sw + j
                    if abs(s2 - su) > a.reach * rad(u[4]): continue
                    n2 = w[2] - j * r * a.q
                    if n2 + 1 != (r // u[1]) * (u[2] + 1): continue
                    lines.append('%s %d %d %.17g %.17g %s %d %d %.17g %.17g' % (w[0], r, n2, s2.real, s2.imag, u[0], u[1], u[2], su.real, su.imag))
    print('level %d: %d components, %d tuner candidates, %d pairs' % (r, len(W), len(U), len(lines)), file=sys.stderr)
    out = subprocess.run([BIN, str(a.p), str(a.q), 'tuned', str(a.threads)], input='\n'.join(lines) + '\n',
                         capture_output=True, text=True, env=env).stdout
    best = {}
    for l in out.splitlines():
        f = l.split(); ratio = float(f[2])
        if f[0] not in best or ratio < best[f[0]][1]: best[f[0]] = (f[1], ratio)
    Cof = {c[0]: (c[4], c[5]) for c in comps}
    tuned = {k: v for k, v in best.items() if v[1] <= 1.5}
    hist = defaultdict(int)
    for k, (u, ratio) in best.items(): hist[min(8, int(math.floor(math.log10(ratio)))) if ratio > 0 and ratio < float('inf') else 'inf'] += 1
    print('ratio histogram (log10 floor):', dict(hist))
    tot = sum(Cof[k][0] for k in tuned if not Cof[k][1])
    print('tuned: %d components (primitive mass %.6e)' % (len(tuned), tot))
    for k, (u, ratio) in sorted(tuned.items(), key=lambda kv: -Cof[kv[0]][0])[:15]:
        print('  %s C %.4e sat %d  by %s (C %.3e)  ratio %.3f' % (k, Cof[k][0], Cof[k][1], u, Cof[u][0], ratio))
    with open('%s_tuned_r%d.txt' % (a.prefix, r), 'w') as fo:
        for k, (u, ratio) in tuned.items(): fo.write('%s %s %.4g\n' % (k, u, ratio))

if __name__ == '__main__':
    main()
