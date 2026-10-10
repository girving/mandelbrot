"""Near-diagonal seahorse families: in limb k, the NRP with kneading tail 1^(2k+d) (a second pass about as long as the
first), for fixed offset d.  Its angle keys in limb k+1 are those of limb k with one 2-bit insertion (the same pattern
for every k), so keys for any k follow from two small limbs.

  python3 diagonal.py limbA.txt A limbB.txt B  d1,d2,... k1,k2,... > lines   (B = A + 1; lines for limb_families
                                                                             custom: "D<d>_k lo hi k")"""
import sys
from family_wfa import kneading
from digits import insertion

def ones_keys(path, k):
    head, pre = '1' * (2 * k) + '0', '01' * (k - 1)
    out = {}
    for l in open(path):
        j, lo, hi = l.split(); t = kneading(pre + lo)[len(head):]
        if t and set(t) == {'1'}: out[len(t) - 2 * k] = (lo, hi)
    return out

def keys_at(kA, KA, KB, d, k):
    """Keys of offset d in limb k, from limbs kA and kA + 1"""
    out = []
    for w in (0, 1):
        a, b = KA[d][w], KB[d][w]
        i = insertion(a, b)[0]; ins = b[i:i + 2]
        out.append(a[:i] + ins * (k - kA) + a[i:])
    return tuple(out)

if __name__ == '__main__':
    kA = int(sys.argv[2]); KA = ones_keys(sys.argv[1], kA); KB = ones_keys(sys.argv[3], int(sys.argv[4]))
    for d in map(int, sys.argv[5].split(',')):
        for k in map(int, sys.argv[6].split(',')):
            lo, hi = keys_at(kA, KA, KB, d, k)
            print('D%d %s %s %d' % (d, lo, hi, k))
