"""Per-limb NRP/bulb ratios κ in the large-k limit from strip censuses (gl_pipeline.sh logs).

For each root and limb t of its strip: κ_t = Σ_r mass_r(t) / C_bulb(t) through the computed levels, the ratios
between consecutive levels, and a geometric tail from the last ratio (a crude estimate, labelled as such).

  python3 gl_kappa.py pipeline_log ...   (each log: gl_pipeline.sh output for one root)"""
import re, sys
from collections import defaultdict

def parse(path):
    lev, cur, single, bulbs = defaultdict(dict), None, None, {}
    root = None
    for l in open(path):
        m = re.search(r'gl_pipeline.sh (\d+) (\d+) (-?\d+)', l)
        if m: root = '%s/%s' % (m.group(1), m.group(2))
        m = re.search(r'\] level (\d+) limbs', l)
        if m: cur = int(m.group(1))
        m = re.search(r'half \(cut .*\): \d+ single-transit families, S_0 = ([0-9.e+-]+); satellites ([0-9.e+-]+)', l)
        if m and single is None: single = float(m.group(1)); bulbs['0'] = float(m.group(2))
        m = re.match(r'\s+limb (\S+)\s+\d+ components\s+([0-9.e+-]+)', l)
        if m and cur: lev[cur][m.group(1)] = float(m.group(2))
        m = re.match(r'\s+sat ([0-9.e+-]+) bulb of (\S+)', l)
        if m: bulbs.setdefault(m.group(2), float(m.group(1)))
    lev[1]['0'] = single
    return root, lev, bulbs

def main():
    for path in sys.argv[1:]:
        root, lev, bulbs = parse(path)
        print('root %s: bulbs %s' % (root, ' '.join('%s %.4e' % kv for kv in sorted(bulbs.items()))))
        limbs = sorted({t for d in lev.values() for t in d}, key=lambda t: (len(t), t))
        for t in limbs:
            if t not in bulbs: continue
            ms = [lev[r].get(t, 0.0) for r in sorted(lev)]
            tot = sum(ms)
            nz = [m for m in ms if m > 0]
            ratio = nz[-1] / nz[-2] if len(nz) >= 2 else float('nan')
            tail = nz[-1] * ratio / (1 - ratio) if 0 < ratio < 1 else float('nan')
            print('  limb %-4s κ (levels %d..%d) = %.4e   + geometric tail ~%.1e (ratio %.2f)   levels %s' % (
                t, min(lev), max(lev), tot / bulbs[t], tail / bulbs[t], ratio, ' '.join('%.3e' % m for m in ms)))

if __name__ == '__main__':
    main()
