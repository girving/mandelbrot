"""Size-engine lines (limb_families size) for seahorse families given by run-length digit words (digits.py).

  python3 family_jobs.py keys K0 SPEC > lines     # SPEC: a file of lines "s d1 d2 ... " with one digit '*', then
                                                  # the n values for it after ':' (e.g. "1 3 * : 2 4 6 8 ... 64")

Each word gets K = max(16, ⌈Σd/2⌉ + 1) (inside the stable range: tail length ≤ 2k) and k dense from K to K + 16, then
geometric (×1.05) to 16K.  Lines are "name lo hi k1,k2,..." with name = s_d1_d2_..., the angle keys from keys_for."""
import sys
from digits import census, keys_for

def ks(K):
    out = list(range(K, K + 17)); k = K + 16
    while k < 16 * K:
        k = max(k + 1, int(k * 1.05 + 0.5)); out.append(k)
    return out

def name(d):
    return '_'.join(str(x) for x in d)

if __name__ == '__main__':
    cs = census(sys.argv[1], int(sys.argv[2]))
    for l in open(sys.argv[3]):
        if not l.strip() or l.startswith('#'): continue
        head, ns = l.split(':')
        sk = head.split()
        for n in map(int, ns.split()):
            d = (sk[0],) + tuple(n if x == '*' else int(x) for x in sk[1:])
            lo, hi = keys_for(cs, d)
            K = max(16, (sum(d[1:]) + 1) // 2 + 1)
            print('%s %s %s %s' % (name(d), lo, hi, ','.join(map(str, ks(K)))))
