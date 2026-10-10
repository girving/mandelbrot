"""How many model components (= M's k-families) are needed to resolve mass η: N(η) = #{W : C_W > η} over a census
directory (c_r2..c_rL, primitive, deduplicated by centre), per level, with the local slope d log N / d log(1/η).
  python3 count_exponent.py dir [--levels 7]"""
import argparse, collections, math
ap = argparse.ArgumentParser(); ap.add_argument('dir'); ap.add_argument('--levels', type=int, default=7); a = ap.parse_args()
comps = {}
for r in range(2, a.levels + 1):
    for l in open('%s/c_r%d.txt' % (a.dir, r)):
        f = l.split()
        if f[6] != '0': continue
        key = (round(float(f[3]), 9), round(float(f[4]), 9))
        comps[key] = max(comps.get(key, (0, 0)), (float(f[5]), r))
print('%d distinct primitive components' % len(comps))
prev = None
for e in range(8, 25):
    eta = 10 ** (-e / 2); cnt = collections.Counter(r for C, r in comps.values() if C > eta); N = sum(cnt.values())
    sl = '%.3f' % (math.log(N / prev) / math.log(10 ** 0.5)) if prev else ''
    print('  %.0e %6d  %s   %s' % (eta, N, ' '.join('%5d' % cnt[r] for r in range(2, a.levels + 1)), sl)); prev = N
