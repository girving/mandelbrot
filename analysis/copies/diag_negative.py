"""Two-transit labels (u, form, d, v) at negative offsets: for d ≤ -(|u| + |v|) the word u R v (R a run of 2k + d:
1s for form '1', an alternation starting with x for form 'ax') is a stable word in limb k, so keys_for gives its keys.

  python3 diag_negative.py keys K0 shapes d_min ks > lines   (shapes: lines "u|form|v"; names N<i>d<d>; lines for
                                                             limb_families custom: "name lo hi k")"""
import sys
from family_wfa import kneading
from digits import census, digits, keys_for

def run(form, L):
    if form == '1': return '1' * L
    x = form[1]; y = '0' if x == '1' else '1'
    return ((x + y) * L)[:L]

if __name__ == '__main__':
    cs = census(sys.argv[1], int(sys.argv[2]))
    shapes = [l.strip() for l in open(sys.argv[3]) if l.strip()]
    d_min = int(sys.argv[4]); ks = list(map(int, sys.argv[5].split(',')))
    for i, sh in enumerate(shapes):
        u, form, v = sh.split('|'); u = '' if u == '-' else u; v = '' if v == '-' else v
        for d in range(-(len(u) + len(v)) - 1, d_min - 1, -1):
            for k in ks:
                L = 2 * k + d
                if L < 2: continue
                t = u + run(form, L) + v
                try: lo, hi = keys_for(cs, digits(t))
                except (KeyError, IndexError): continue
                print('N%dd%d %s %s %d' % (i, d, lo, hi, k))
