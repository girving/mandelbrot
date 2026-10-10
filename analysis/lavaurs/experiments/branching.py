"""Branching factor and scaling of the near-critical-value refinement (q = 2, gate -1).

For an island U (level r, centre σ_U, ε_U = 1/|Θ_r'(σ_U)|), its grandchildren G (one transit after U's island) are
atoms of the frozen centre measure whose δ-series around σ_U has radius D = |y_t - v_U| (y_t the target, v_U = ζ0 + σ_U,
mod 1/q and τ).  With ρ = |σ_G - σ_U| / D, "bad" grandchildren (ρ > 1/2) need their own island (their parent child Y)
as expansion base.  Primitive components only (satellites, cusp ≥ 0.01 from exact areas, are dropped and not
expanded).  Prints per island: C_U, ε_U, children and grandchildren mass / C_U, the bad share (ρ > 1/2, > 1), and
B90, B99 = the number of children carrying 90%, 99% of the bad mass (the sub-islands to process).  Islands: the
heaviest level-1 components, then down the heaviest-child lineage to --depth.  Children are island children
(glavaurs ichildren: continuation from the island's branch point), or with --census the census's (Newton from wide
starts, which also finds solutions in other islands).

  python3 branching.py single.txt [--a 2] [--per 3] [--depth 4] [--y 20] [--targets 60]"""
import argparse, math, os, subprocess, sys
from collections import defaultdict
B = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../build/release/glavaurs')
env = dict(os.environ, GL_SIDE='-1')
K = 2.4674011
TAU = 2j * math.pi * 1.375 / 2
def run(args, text, starts=26):
    return subprocess.run([B, '1', '2'] + args, input=text, capture_output=True, text=True,
                          env=dict(env, CHILDREN_STARTS=str(starts))).stdout
def lat(z):   # distance to the deck lattice 1/2 Z + τ Z (q = 2)
    return min(abs((lambda w: w - round(w.real * 2) / 2)(z - k * TAU)) for k in range(-2, 3))
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('single'); ap.add_argument('--a', type=int, default=2); ap.add_argument('--per', type=int, default=3)
    ap.add_argument('--depth', type=int, default=4); ap.add_argument('--y', type=int, default=20)
    ap.add_argument('--targets', type=int, default=60); ap.add_argument('--threads', type=int, default=2)
    ap.add_argument('--workdir', default='branching_work'); ap.add_argument('--verbose', action='store_true')
    ap.add_argument('--census', action='store_true')   # census children (Newton from wide starts) instead of island children
    a = ap.parse_args()
    comps = []
    for l in open(a.single):
        f = l.split(); sat = f[0].startswith('sat'); n = int(f[0][3:]) if sat else int(f[0])
        comps.append(('B%d' % len(comps) if sat else 'F%d' % len(comps), n, complex(float(f[1]), float(f[2])), float(f[3]), float(f[5]), sat))
    cut = (-0.16135210336424927 - 4.1583377953217475) / 2
    own = [c for c in comps if c[2].imag > cut]
    targets = [c for c in own if c[5]] + sorted([c for c in own if not c[5]], key=lambda c: -c[3])[:a.targets]
    Ct = {t[0]: t[4] for t in targets}; nt = {t[0]: t[1] for t in targets}
    rad = lambda C: math.sqrt(C / (K * math.pi))
    def children(srcs, r):
        """srcs: [(name, n, σ, C)] at level r-1 -> {name: [(name, n, σ, C)]}, primitive, deduplicated across srcs"""
        seen = set()
        lines = []
        for nm, n, s, C in srcs:
            jc = -(n + 1)
            for tn, ntt, ct, _, _, _ in targets:
                if a.census:
                    lines.append('%s~%s %d %d %.17g %.17g %.6g %d %.17g %.17g %d %d' % (nm, tn, r, n, s.real, s.imag, rad(C), ntt, ct.real, ct.imag, jc - 16, jc + 4))
                else:
                    lines.append('%s~%s %d %.17g %.17g %d %.17g %.17g %d %d' % (nm, tn, r, s.real, s.imag, ntt, ct.real, ct.imag, jc - 16, jc + 4))
        raw = defaultdict(list)
        for l in run(['children' if a.census else 'ichildren', str(a.threads)], '\n'.join(lines) + '\n').splitlines():
            f = l.split(); nm, br, j = f[0].rsplit('|', 2); src, t = nm.rsplit('~', 1)
            c = complex(float(f[1]), float(f[2])); w = float(f[3])
            if (a.census and int(f[4])) or w >= 1 or c.imag <= cut: continue
            k = (r, round(c.real % 1, 8), round(c.imag, 8))
            if k in seen: continue
            seen.add(k); raw[src].append(['%s/%d' % (src, len(raw[src])), nt[t] - int(j), c, Ct[t] * w])
        heavy = [k for v in raw.values() for k in v if k[3] > 1e-10]
        res = run(['area', str(a.threads)], ''.join('H%d %d %d %.17g %.17g\n' % (i, r, k[1], k[2].real, k[2].imag) for i, k in enumerate(heavy)))
        drop = set()
        for l in res.splitlines():
            f = l.split()
            if f[3] == 'failed': continue
            k = heavy[int(f[0][1:])]
            k[3] = float(f[6])
            if float(f[8]) >= 0.01: drop.add(k[0])
        res = {s: [tuple(k) for k in v if k[0] not in drop] for s, v in raw.items()}
        for v in res.values():
            for k in v: known[r][k[0]] = k
        return res
    known = defaultdict(dict)
    wd = a.workdir
    os.makedirs(wd, exist_ok=True)
    def tuned(R):
        """names at level R tuned by known components of lower levels (gl_tune --rule block)"""
        if R < 2: return set()
        for k in range(2, R + 1):
            with open('%s/c_r%d.txt' % (wd, k), 'w') as fo:
                for nm, n, sg, C in known[k].values(): fo.write('%s %d %d %.17g %.17g %.6e 0\n' % (nm, k, n, sg.real, sg.imag, C))
        subprocess.run(['python3', os.path.join(os.path.dirname(os.path.abspath(__file__)), '../gl_tune.py'), '1', '2', '-1',
                        a.single, wd + '/c', '--level', str(R), '--threads', str(a.threads), '--rule', 'block'],
                       capture_output=True, text=True, env=env)
        return {l.split()[0] for l in open('%s/c_tuned_r%d.txt' % (wd, R))}
    def locate(r, sigmas):
        res = run(['locate'], ''.join('P%d %d %.17g %.17g\n' % (i, r, s.real, s.imag) for i, s in enumerate(sigmas)))
        th, d = {}, {}
        for l in res.splitlines():
            f = l.split()
            if len(f) > 5: th[int(f[0][1:])] = complex(float(f[2]), float(f[3])); d[int(f[0][1:])] = complex(float(f[4]), float(f[5]))
        return th, d
    def measure(r, U):
        kids = children([U], r + 1).get(U[0], [])
        t1 = tuned(r + 1)
        kids = [k for k in kids if k[0] not in t1]
        top = sorted(kids, key=lambda k: -k[3])[:a.y]
        gk = children(top, r + 2)
        t2 = tuned(r + 2)
        gk = {y: [g for g in v if g[0] not in t2] for y, v in gk.items()}
        _, dth = locate(r, [U[2]]); eps = 1 / abs(dth[0]) if 0 in dth else float('nan')
        G = [(g, y[0]) for y in top for g in gk.get(y[0], [])]
        th, _ = locate(r + 2, [g[2] for g, _ in G])
        tot = bad5 = bad1 = 0.0; badby = defaultdict(float)
        for i, (g, yn) in enumerate(G):
            if i not in th: continue
            D = lat(th[i] - U[2]); dd = g[2] - U[2]; dd = complex((dd.real + 0.5) % 1 - 0.5, dd.imag)
            rho = abs(dd) / D
            tot += g[3]
            if a.verbose and g[3] > 1e-3 * U[3]:
                print('    G %-22s C/C_U %.2e  d %.3e  D %.3e  rho %.2f  y-σ_U %+.4f%+.4fi' % (g[0], g[3] / U[3], abs(dd), D, rho, (th[i] - U[2]).real, (th[i] - U[2]).imag))
            if rho > 0.5: bad5 += g[3]; badby[yn] += g[3]
            if rho > 1: bad1 += g[3]
        bm = sorted(badby.values(), reverse=True); acc = 0; b90 = b99 = 0
        for k, m in enumerate(bm, 1):
            acc += m
            if not b90 and acc >= 0.9 * bad5: b90 = k
            if not b99 and acc >= 0.99 * bad5: b99 = k
        km = sum(k[3] for k in kids); cov = sum(k[3] for k in top) / km if km else 0
        print('L%d %-16s C %.3e eps %.4f  children %4d (mass/C %.3f, top %d = %.0f%%)  grandchildren mass/C %.4f'
              '  bad rho>.5 %.2e rho>1 %.2e  B90 %d B99 %d of %d' % (r, U[0], U[3], eps, len(kids), km / U[3], len(top),
              100 * cov, tot / U[3], bad5 / tot if tot else 0, bad1 / tot if tot else 0, b90, b99, len(bm)), flush=True)
        return top
    A = [(c[0], c[1], c[2], c[3]) for c in sorted([c for c in own if not c[5]], key=lambda c: -c[3])[:a.a]]
    level = A
    for r in range(1, a.depth + 1):
        nxt = []
        for U in level:
            top = measure(r, U)
            nxt += top[:a.per]
        uniq = {}
        for k in nxt: uniq.setdefault((round(k[2].real % 1, 7), round(k[2].imag, 7)), k)
        level = sorted(uniq.values(), key=lambda k: -k[3])[:max(a.per, 2)]
main()
