"""Farey limbs as copies of the bulb limb, tail by tail.

A component of the limb t = b/m has kneading 1^(q-1) 0 τ; the same tail τ after the bulb limb's first 0 gives a bulb-
limb component (for one-transit families τ is the census word).  This computes τ for every component (limb_class.tail),
pairs components across limbs by τ, and reports per limb the ratios ρ_F(τ) = C_F(τ)/C_bulb(τ): how much of the limb's
mass has a bulb-limb partner, and how ρ depends on the tail (its first run, the number of gate passages).

  python3 farey_tails.py --cache tails.txt --census model_census.txt.gz --tune tune --default ... --island ... [--min 1e-12]"""
import argparse, gzip, os, sys, math
from collections import defaultdict
from multiprocessing import Pool

def opened(p): return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)

def key(c): return (round(c.real % 1.0, 8) % 1.0, round(c.imag, 8))

def _work(args):
    from limb_class import tail
    k, s, n, r = args
    try:
        t = tail(s, n, r)
        return k, None if t is None else (t[0], t[1], ' '.join(t[3]))
    except Exception:
        return k, None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cache', required=True)
    ap.add_argument('--census', required=True)
    ap.add_argument('--tune', default='')
    ap.add_argument('--default', default='')
    ap.add_argument('--island', default='')
    ap.add_argument('--procs', type=int, default=4)
    ap.add_argument('--min', type=float, default=1e-12)
    a = ap.parse_args()
    tuned = set()
    for p in filter(None, a.tune.split(',')):
        for l in opened(p):
            f = l.split()
            if '*' in f[0] and f[3] != 'failed': tuned.add(key(complex(float(f[3]), float(f[4]))))
    comp = {}
    def add(f, cusp):
        if f[3] == 'failed' or cusp >= 1e-8: return
        c = complex(float(f[3]), float(f[4])); k = key(c)
        if k in tuned or k in comp: return
        comp[k] = (float(f[7]), int(f[2]), int(f[1]), c)
    for l in opened(a.census):
        f = l.split()
        if f[3] != 'failed': add(f, 0.0)
    for p in filter(None, a.default.split(',')):
        for l in opened(p):
            f = l.split()
            if len(f) >= 10: add(f, float(f[9]) if f[3] != 'failed' else 1)
    for p in filter(None, a.island.split(',')):
        for l in opened(p):
            f = l.split()
            if len(f) >= 16: add(f, float(f[-1]))
    cache = {}
    if os.path.exists(a.cache):
        for l in open(a.cache):
            f = l.rstrip('\n').split('\t')
            cache[(float(f[0]), float(f[1]))] = None if f[2] == '-' else (int(f[2]), int(f[3]), f[4])
    todo = sorted(((k, v[3], v[1], v[2]) for k, v in comp.items() if k not in cache and v[0] > a.min), key=lambda t: -comp[t[0]][0])
    print('%d components above %g, %d to compute' % (sum(1 for v in comp.values() if v[0] > a.min), a.min, len(todo)), file=sys.stderr)
    with open(a.cache, 'a') as out, Pool(a.procs) as pool:
        for i, (k, res) in enumerate(pool.imap_unordered(_work, todo, chunksize=16)):
            cache[k] = res
            out.write('%r\t%r\t%s\n' % (k[0], k[1], '-' if res is None else '%d\t%d\t%s' % res))
            if i % 5000 == 0: out.flush(); print('  %d / %d' % (i, len(todo)), file=sys.stderr)
    # group by limb and tail
    by = defaultdict(dict)   # (m, b) -> tail -> C (sum, in case of duplicates)
    for k, v in comp.items():
        res = cache.get(k)
        if res is None or v[0] <= a.min: continue
        m, b, t = res
        by[(m, b)][t] = by[(m, b)].get(t, 0.0) + v[0]
    bulb = by[(1, 0)]
    print('bulb limb: %d tails, mass %.6e' % (len(bulb), sum(bulb.values())))
    for (m, b) in sorted(by):
        if (m, b) == (1, 0): continue
        d = by[(m, b)]; tot = sum(d.values())
        pairs = [(t, C, bulb[t]) for t, C in d.items() if t in bulb]
        pm = sum(C for _, C, _ in pairs); pb = sum(B for _, _, B in pairs)
        print('t = %d/%d: %d tails, mass %.6e; paired %.1f%% of mass (bulb partners %.6e), mass ratio %.5f' % (
            b, m, len(d), tot, 100 * pm / tot if tot else 0, pb, pm / pb if pb else float('nan')))
        # ratio by the first token (first run) and by number of gate passages in the tail
        grp = defaultdict(lambda: [0.0, 0.0, 0])
        for t, C, B in pairs:
            gates = sum(int(x.split(':')[0]) for x in t.split() if ':' in x)
            g = grp[('gates', gates)]; g[0] += C; g[1] += B; g[2] += 1
        for kk in sorted(grp):
            C, B, n = grp[kk]
            print('     tails with %d gate passages: %5d pairs, ratio %.5f (limb %.3e, bulb %.3e)' % (kk[1], n, C / B, C, B))
        top = sorted(pairs, key=lambda x: -x[2])[:6]
        print('     heaviest partners: ' + '; '.join('%s: %.4f' % (t, C / B) for t, C, B in top))

if __name__ == '__main__':
    main()
