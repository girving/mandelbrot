"""Does the frozen-σ horn-map dynamics (a transfer operator) reproduce the parameter census?

For the level-r census of gl_census.py (same sources, targets, shifts and starts), run `glavaurs dchildren` (children
as preimages under p ↦ H(p) + σ_X at the source's frozen σ_X) and compare with the parameter children (the census's
raw output, Newton on Θ_r):
  child by child: w_param / w_T (transversality factor kept) and w_param / w_deep (purely multiplicative weight);
  per source: R = Σ_children C_t w / C_X under each weighting, and the scatter of R across sources (the old
  "universal kernel" asked for R = const).

  python3 gl_opdiag.py p q gate single.txt prefix --level r [census options as used]"""
import argparse, math, os, re, subprocess, sys
from collections import defaultdict
import numpy as np

BIN = os.path.join(os.path.dirname(__file__), '../../build/release/glavaurs')

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('p', type=int); ap.add_argument('q', type=int); ap.add_argument('gate', type=int)
    ap.add_argument('single'); ap.add_argument('prefix')
    ap.add_argument('--level', type=int, default=2)
    ap.add_argument('--sources', type=int, default=300)
    ap.add_argument('--targets', type=int, default=100)
    ap.add_argument('--src-min', type=float, default=1e-10)
    ap.add_argument('--jd', type=float, default=8); ap.add_argument('--ju', type=float, default=2)
    ap.add_argument('--big-src', type=float, default=1e-9)
    ap.add_argument('--starts', type=int, default=26)
    ap.add_argument('--small-targets', type=int, default=30)
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--max-sources', type=int, default=10**9)   # diagnose only the heaviest sources
    ap.add_argument('--cache', default='')                      # dchildren output file (reused if present)
    ap.add_argument('--own-param', action='store_true')         # compute the parameter children here (cached) instead of
                                                                # the census raw file; no census verdicts (all children
                                                                # with w < 1 outside 4 radii are compared)
    a = ap.parse_args()
    q, r = a.q, a.level
    env = dict(os.environ, GL_SIDE=str(a.gate))
    head = subprocess.run([BIN, str(a.p), str(q), 'consist'], capture_output=True, text=True, env=env).stderr
    K = float(re.search(r'K = ([0-9.e+-]+)', head).group(1))
    radius = lambda C: math.sqrt(C / (K * math.pi))
    fams, sats = [], []
    for l in open(a.single):
        f = l.split()
        if f[0].startswith('sat'): sats.append(('B%d' % len(sats), int(f[0][3:]), complex(float(f[1]), float(f[2])), float(f[3])))
        else: fams.append(('F%d' % len(fams), int(f[0]), complex(float(f[1]), float(f[2])), float(f[3])))
    big2 = sorted(sats, key=lambda x: -x[3])[:2]
    cut = (big2[0][2].imag + big2[1][2].imag) / 2
    upper = [b for b in big2 if b[1] == q - 1][0][2].imag > cut
    inh = lambda c: (c.imag >= cut) == upper
    fams = sorted([f for f in fams if inh(f[2])], key=lambda x: -x[3]); sats = [f for f in sats if inh(f[2])]
    targets = sats + fams[:a.targets]
    Ct = {t[0]: t[3] for t in targets}
    if r == 2:
        srcs = [(nm, n, c, C, False) for nm, n, c, C in fams[:a.sources]] + [(nm, n, c, C, True) for nm, n, c, C in sats]
    else:
        tuned = set()
        tf = '%s_tuned_r%d.txt' % (a.prefix, r - 1)
        if os.path.exists(tf): tuned = {l.split()[0] for l in open(tf)}
        srcs = []
        for l in open('%s_r%d.txt' % (a.prefix, r - 1)):
            f = l.split(); C = float(f[5])
            if f[0] in tuned or (f[6] == '0' and C <= a.src_min): continue
            srcs.append((f[0], int(f[2]), complex(float(f[3]), float(f[4])), C, f[6] == '1'))
    srcs = sorted(srcs, key=lambda s: -s[3])[:a.max_sources]
    keep = {s[0] for s in srcs}
    Cs = {s[0]: s[3] for s in srcs}
    jd, ju = int(round(a.jd * q)), int(round(a.ju * q))
    big, small = [], []
    for nm, n, c, C, sat in srcs:
        jc = -(n + 1)
        heavy = sat or C > a.big_src
        for tn, nt, ct, C_t in (targets if heavy else targets[:len(sats) + a.small_targets]):
            (big if heavy else small).append('%s~%s %d %d %.17g %.17g %.6g %d %.17g %.17g %d %d' % (
                nm, tn, r, n, c.real, c.imag, radius(C), nt, ct.real, ct.imag, jc - jd, jc + ju))
    def run(lines, starts, use_eps):
        e = dict(env, CHILDREN_STARTS=str(starts), DCH_EPS=str(use_eps))
        return subprocess.run([BIN, str(a.p), str(q), 'dchildren', str(a.threads)], input='\n'.join(lines) + '\n',
                              capture_output=True, text=True, env=e).stdout
    dyns = {}
    for use_eps in (1, 0):
        cf = '%s.eps%d' % (a.cache, use_eps) if a.cache else ''
        if cf and os.path.exists(cf): dyns[use_eps] = open(cf).read()
        else:
            dyns[use_eps] = run(big, a.starts, use_eps) + run(small, 0, use_eps)
            if cf: open(cf, 'w').write(dyns[use_eps])
    # the census's verdicts: primitive, not tuned, with exact or formula C
    final = {}
    tuned = set()
    tf = '%s_tuned_r%d.txt' % (a.prefix, r)
    if os.path.exists(tf): tuned = {l.split()[0] for l in open(tf)}
    if not a.own_param:
        for l in open('%s_r%d.txt' % (a.prefix, r)):
            f = l.split(); c = complex(float(f[3]), float(f[4]))
            final[(round(c.real % 1, 7) % 1, round(c.imag, 7))] = (float(f[5]), f[6] == '1' or f[0] in tuned)
    nkey = lambda c: (round(c.real % 1, 7) % 1, round(c.imag, 7))
    # parameter children from the census's raw file
    par = defaultdict(list)
    if a.own_param:
        pf = '%s.param' % a.cache
        if os.path.exists(pf): ptext = open(pf).read()
        else:
            def prun(lines, starts):
                e = dict(env, CHILDREN_STARTS=str(starts))
                return subprocess.run([BIN, str(a.p), str(q), 'children', str(a.threads)], input='\n'.join(lines) + '\n',
                                      capture_output=True, text=True, env=e).stdout
            ptext = prun(big, a.starts) + prun(small, 0)
            open(pf, 'w').write(ptext)
        plines = ptext.splitlines()
    else: plines = open('%s_raw_r%d.txt' % (a.prefix, r))
    for l in plines:
        f = l.split(); key = f[0].rsplit('|', 2); src = key[0].split('~')[0]
        if src not in keep: continue
        w = float(f[3]); c = complex(float(f[1]), float(f[2]))
        if int(f[4]) or w >= 1: continue
        fk = final.get(nkey(c)) if not a.own_param else (float('nan'), False)
        if fk is None or fk[1]: continue   # not kept by the census (or a satellite / tuning)
        par[(key[0], key[2])].append((c, w, fk[0]))
    srcpos = {x[0]: x[2] for x in srcs}
    def wq(x, w, qs=(0.1, 0.5, 0.9)):
        o = np.argsort(x); cw = np.cumsum(w[o]) / w.sum()
        return [x[o][min(len(x) - 1, np.searchsorted(cw, qq))] for qq in qs]
    first = True
    for use_eps, label in ((1, 'one-step map with ε (H + ε(p - p*))'), (0, 'frozen σ, multiplicative (ε = 0)')):
        dy = defaultdict(list)
        for l in dyns[use_eps].splitlines():
            f = l.split(); key = f[0].rsplit('|', 2)
            w = float(f[3])
            if int(f[5]) or w >= 1: continue
            dy[(key[0], key[2])].append((complex(float(f[1]), float(f[2])), w))
        rows, unmatched_p, seen, usedd = [], 0.0, defaultdict(set), defaultdict(set)
        for k, P in par.items():
            st, j = k; src, tg = st.split('~'); C_t = Ct[tg]
            D = dy.get(k, [])
            for s_, w, Cx in P:
                kk = nkey(s_)
                if kk in seen[src]: continue
                seen[src].add(kk)
                best = min(range(len(D)), key=lambda i: abs(D[i][0] - s_)) if D else None
                if best is None or abs(D[best][0] - s_) > max(0.05 * abs(s_ - srcpos[src]), 1e-7):
                    unmatched_p += C_t * w; continue
                usedd[k].add(best)
                rows.append((src, C_t * w, C_t * D[best][1], abs(D[best][0] - s_) / abs(s_ - srcpos[src]), Cx))
        unmatched_d = sum(Ct[k[0].split('~')[1]] * D[i][1] for k, D in dy.items() for i in range(len(D)) if i not in usedd[k])
        P_ = np.array([x[1] for x in rows]); M_ = np.array([x[2] for x in rows]); dist = np.array([x[3] for x in rows])
        print('%s:' % label)
        print('  %d matched children (parameter mass %.6e, model %.6e); unmatched: parameter %.3e, model %.3e' % (
            len(rows), P_.sum(), M_.sum(), unmatched_p, unmatched_d))
        print('  child log(w_param/w_model), mass-weighted 10/50/90%%: %s;  relative position error 50/90%%: %s' % (
            ' '.join('%+.2e' % v for v in wq(np.log(P_ / M_), P_)), ' '.join('%.1e' % v for v in wq(dist, P_, (0.5, 0.9)))))
        RP, RM, RX = defaultdict(float), defaultdict(float), defaultdict(float)
        for src, p_, m_, _, cx in rows: RP[src] += p_; RM[src] += m_; RX[src] += cx
        sl = [x for x in srcs if RP[x[0]] > 0]
        wts = np.array([RP[x[0]] for x in sl])
        def spread(vals):
            v = np.log(np.array(vals)); lo, md, hi = wq(v, wts)
            return 'median %.4g, 10-90%% range x%.3f' % (math.exp(md), math.exp(hi - lo))
        if first: print('  per source W/C_X with the census constants (the old universal-kernel test): %s' % spread([RX[x[0]] / x[3] for x in sl]))
        print('  per source W_param / W_model: %s' % spread([RP[x[0]] / RM[x[0]] for x in sl]))
        if use_eps:
            # per source: |ε_X| = 1/|Θ'_{r-1}(σ_X)| against the source's agreement
            loc = subprocess.run([BIN, str(a.p), str(q), 'locate'], input=''.join('%s %d %.17g %.17g\n' % (
                x[0], r - 1, x[2].real, x[2].imag) for x in sl), capture_output=True, text=True, env=env).stdout
            eps = {}
            for l in loc.splitlines():
                f = l.split()
                if len(f) >= 6 and f[2] != 'failed': eps[f[0]] = 1 / abs(complex(float(f[4]), float(f[5])))
            per = defaultdict(list)
            for src, p_, m_, d_, cx in rows: per[src].append((p_, math.log(p_ / m_)))
            print('    source        C_X     |eps_X|   W_param/W_model   child |log ratio| (mass-weighted mean)   children')
            for x in sorted(sl, key=lambda x: -x[3])[:20]:
                v = per[x[0]]; ww = sum(t[0] for t in v)
                print('    %-12s %.3e %.3e   %.6f          %.2e                               %d' % (
                    x[0], x[3], eps.get(x[0], float('nan')), RP[x[0]] / RM[x[0]], sum(t[0] * abs(t[1]) for t in v) / ww, len(v)))
        first = False

if __name__ == '__main__':
    main()
