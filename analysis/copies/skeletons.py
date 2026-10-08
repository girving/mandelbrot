"""Skeleton classes of seahorse family words: a word's skeleton is the word with its first maximal digit replaced by
'*', and its class the skeleton with that digit's parity.  Every word belongs to exactly one class, so the sum over
words beyond the census is a sum over classes of sums over the long digit n.

  python3 skeletons.py keys K0 constants N > spec   # the N heaviest classes (by census mass over words with digit
                                                    # sum ≥ 9), as family_jobs.py spec lines with n up to 512"""
import sys
from collections import defaultdict
from digits import census

def skeleton(d):
    p = max(range(1, len(d)), key=lambda q: (d[q], -q))
    return d[:p] + ('*',) + d[p + 1:], d[p] % 2

LARGE = [72, 80, 96, 112, 128, 160, 192, 256, 320, 384, 512]

if __name__ == '__main__':
    cs = census(sys.argv[1], int(sys.argv[2]))
    C = {(int(r[0]), int(r[1])): float(r[2]) for r in (l.split() for l in open(sys.argv[3]))}
    mass = defaultdict(float)
    for d, (j, i, lo, hi, t) in cs.items():
        if len(d) < 2 or sum(d[1:]) < 9: continue
        mass[skeleton(d)] += C[(j, i)]
    top = sorted(mass.items(), key=lambda x: -x[1])[:int(sys.argv[4])]
    tot = sum(mass.values())
    print('# %d classes, %.8f of the census mass at digit sum ≥ 9' % (len(top), sum(m for _, m in top) / tot))
    for (sk, par), m in top:
        lo = 2 if par == 0 else 3
        ns = list(range(lo, 65, 2)) + [n + par for n in LARGE]
        print('%s : %s   # mass %.3e' % (' '.join(map(str, sk)), ' '.join(map(str, ns)), m))
