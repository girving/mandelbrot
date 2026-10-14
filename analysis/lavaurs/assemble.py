"""The strip's NRP sectors S_r (r transits) from every census, merged.

Exact censuses: lavaurs_area default-mode files (cusp in column 10) and --island files (cusp last); fast censuses
(fast_census.py output "name r n re im C sat", areas from the center formula, heavy ones exact).  Components are
deduplicated by center mod 1 (exact values preferred), satellites dropped (cusp ≥ 1e-8 or sat flag), and tunings
U*X removed by center (lavaurs_area LAVAURS_TUNE outputs: "U*X r n re im ... C conv cusp ratio").

  python3 assemble.py --default a,b --island c,d --fast e,f --tune g,h [--top 6]"""
import argparse, gzip
from collections import defaultdict

def opened(p): return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)

def key(c): return (round(c.real % 1.0, 7) % 1.0, round(c.imag, 7))

def main():
    ap = argparse.ArgumentParser()
    for k in ('default', 'island', 'fast', 'tune'): ap.add_argument('--' + k, default='')
    ap.add_argument('--top', type=int, default=6)
    a = ap.parse_args()
    files = lambda s: [p for p in s.split(',') if p]
    tuned = {}
    for p in files(a.tune):
        for l in opened(p):
            f = l.split()
            if '*' in f[0] and len(f) > 7 and f[3] != 'failed':
                tuned[key(complex(float(f[3]), float(f[4])))] = (f[0], int(f[1]), float(f[7]))
    comp = {}   # key -> [r, C, sat, exact, name, center]
    def add(r, c, C, sat, exact, name):
        k = key(c)
        e = comp.get(k)
        if e is None or (exact and not e[3]): comp[k] = [r, C, sat or (e[2] if e else False), exact, name, c]
        elif sat: e[2] = True
    for p in files(a.default):
        for l in opened(p):
            f = l.split()
            if len(f) >= 10 and f[3] != 'failed':
                add(int(f[1]), complex(float(f[3]), float(f[4])), float(f[7]), float(f[9]) >= 1e-8, True, f[0])
    for p in files(a.island):
        for l in opened(p):
            f = l.split()
            if len(f) >= 16:
                add(int(f[1]), complex(float(f[3]), float(f[4])), float(f[7]), float(f[-1]) >= 1e-8, True, f[0])
    for p in files(a.fast):
        for l in opened(p):
            f = l.split()
            if len(f) >= 6:
                add(int(f[1]), complex(float(f[3]), float(f[4])), float(f[5]), len(f) > 6 and f[6] == '1', False, f[0])
    S = defaultdict(float); Tm = defaultdict(float); Sat = defaultdict(float); n = defaultdict(int)
    heavy = defaultdict(list)
    for k, (r, C, sat, exact, name, c) in comp.items():
        if sat: Sat[r] += C; continue
        if k in tuned: Tm[r] += C; continue
        S[r] += C; n[r] += 1
        heavy[r].append((C, name, c, exact))
    print('per transit count r: NRP sector S_{r-1} (primitive, not tuned), tunings removed, satellites')
    for r in sorted(S):
        print('  r = %d: %8d components, S = %.6e   (tunings %.3e, satellites %.3e)' % (r, n[r], S[r], Tm[r], Sat[r]))
        for C, name, c, exact in sorted(heavy[r], key=lambda x: -x[0])[:a.top]:
            print('        %.3e %s σ %.5f%+.5fi  %s' % (C, 'exact' if exact else 'fast ', c.real, c.imag, name[:50]))
    print('tunings known: %d (by r: %s)' % (len(tuned), dict(sorted(defaultdict(int, {}).items()))))

if __name__ == '__main__':
    main()
