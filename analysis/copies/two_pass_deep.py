"""Lines for `limb_families size`: the two-pass sequences (two_pass.py) at large m, with k dense from m+1 to m+16 and
then geometric (×1.05) to 16m.

  python3 two_pass_deep.py M1,M2,... > lines;  limb_families size threads < lines > out"""
import sys
from two_pass import SEQ

def ks(m):
    out = list(range(m + 1, m + 17)); k = m + 16
    while k < 16 * m:
        k = max(k + 1, round(k * 1.05)); out.append(k)
    return out

if __name__ == '__main__':
    for m in map(int, sys.argv[1].split(',')):
        for name, f in SEQ.items():
            lo, hi = f(m)
            print('%sm%d %s %s %s' % (name, m, lo, hi, ','.join(map(str, ks(m)))))
