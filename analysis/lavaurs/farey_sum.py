"""The σ strip's NRP layer summed by limb: classify every primitive, non-tuned component by its Farey limb t = b/m
(limb_class.py) and report the mass per limb and per transit count within each limb.

Inputs (gzipped or plain): lavaurs_area default-mode files ("name r n cre cim ahi alo C conv cusp ...") and --island
files ("... cusp" last), tunings (LAVAURS_TUNE output).  Components are deduplicated mod 1 in Re σ; classification
results are cached in the file named by --cache.

  python3 farey_sum.py --cache cls.txt --tune tune_out.txt --default a.txt,b.txt --island c.txt.gz,... [--procs 4]"""
import argparse, gzip, os, sys
from collections import defaultdict
from multiprocessing import Pool

def opened(p): return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)

def key(c): return (round(c.real % 1.0, 8) % 1.0, round(c.imag, 8))

def _work(args):
    from limb_class import classify
    k, s, n, r = args
    try: return k, classify(s, n, r)
    except Exception: return k, None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cache', required=True)
    ap.add_argument('--tune', default='')
    ap.add_argument('--default', default='')
    ap.add_argument('--island', default='')
    ap.add_argument('--procs', type=int, default=4)
    ap.add_argument('--min', type=float, default=0.0)  # classify only components with C above this
    a = ap.parse_args()
    tuned = set()
    for p in filter(None, a.tune.split(',')):
        for l in opened(p):
            f = l.split()
            if '*' in f[0] and f[3] != 'failed': tuned.add(key(complex(float(f[3]), float(f[4]))))
    comp = {}  # key -> (C, n, r, σ, label?)
    def add(f, cusp):
        if f[3] == 'failed' or cusp >= 1e-8: return
        c = complex(float(f[3]), float(f[4])); k = key(c)
        if k in tuned or k in comp: return
        comp[k] = (float(f[7]), int(f[2]), int(f[1]), c, '~' not in f[0] and '|' not in f[0])
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
            f = l.split()
            cache[(float(f[0]), float(f[1]))] = (f[2], f[3])
    todo = [(k, v[3], v[1], v[2]) for k, v in comp.items() if k not in cache and v[0] > a.min]
    print('%d components (%d by r: %s), %d cached, %d to classify' % (
        len(comp), len(comp), dict(sorted(defaultdict(int, {r: sum(1 for v in comp.values() if v[2] == r) for r in {v[2] for v in comp.values()}}).items())),
        len(cache), len(todo)), file=sys.stderr)
    todo.sort(key=lambda t: -comp[t[0]][0])
    with open(a.cache, 'a') as out, Pool(a.procs) as pool:
        for i, (k, res) in enumerate(pool.imap_unordered(_work, todo, chunksize=16)):
            m, off = ('fail', 'fail') if res is None else (str(res[0]), str(res[1]))
            cache[k] = (m, off)
            out.write('%r %r %s %s\n' % (k[0], k[1], m, off))
            if i % 2000 == 0: out.flush(); print('  %d / %d' % (i, len(todo)), file=sys.stderr)
    # per limb (m, b mod m) and transit count
    limb = defaultdict(lambda: defaultdict(float)); bad = defaultdict(float); labels_off = 0.0
    for k, v in comp.items():
        if k not in cache: bad['unclassified'] += v[0]; continue
        m, off = cache[k]
        if m in ('fail', 'bulb'): bad[m] += v[0]; continue
        m, off = int(m), int(off)
        if (off + 1) % 2: bad['even offset'] += v[0]; continue
        b = ((off + 1) // 2) % m
        if v[4] and m != 1: labels_off += v[0]
        limb[(m, b)][v[2]] += v[0]
    print('per limb (m, b): mass by transit count r, total')
    tot = 0.0
    for (m, b) in sorted(limb):
        row = limb[(m, b)]; t = sum(row.values()); tot += t
        print('  t = %d/%d: %s  total %.6e' % (b, m, '  '.join('r%d %.4e' % (r, row[r]) for r in sorted(row)), t))
    print('classified total %.6e; not counted: %s; labels classified off the bulb limb: %.3e' % (
        tot, ', '.join('%s %.3e' % kv for kv in bad.items()), labels_off))

if __name__ == '__main__':
    main()
