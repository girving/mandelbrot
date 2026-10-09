"""Census of the two-transit sector: NRPs of limb k whose kneading tail is u R v with R a run of length 2k + d
(of 1s, or an alternation 10...), u and v short.  Words are labeled (u, form of R, d, v), independent of k; their
keys in limb k+1 are one 2-bit insertion away from limb k (checked between limbs 4, 5 and 6).

  python3 diag_census.py limb4.txt limb5.txt [limb6.txt] > words   (lines "label lo4 hi4 lo5 hi5")"""
import re, sys
from family_wfa import kneading
from digits import insertion

def label(t, k):
    """(u, form, d, v) for the tail's unique run of length ≥ 2k - 2 beyond... or None if not exactly one"""
    found = []
    for m in re.finditer('1+', t):
        if len(m.group()) >= 2 * k - 2: found.append((m.start(), m.end(), '1'))
    for m in re.finditer('(?:10)+1?|(?:01)+0?', t):
        if len(m.group()) >= 2 * k - 2: found.append((m.start(), m.end(), 'a' + m.group()[0]))
    if len(found) != 1: return None
    a, b, form = found[0]
    return (t[:a], form, b - a - 2 * k, t[b:])

def words(path, k):
    head, pre = '1' * (2 * k) + '0', '01' * (k - 1)
    out = {}
    for l in open(path):
        j, lo, hi = l.split(); t = kneading(pre + lo)[len(head):]
        lab = label(t, k)
        if lab is not None: out.setdefault(lab, []).append((lo, hi))
    return out

if __name__ == '__main__':
    W4, W5 = words(sys.argv[1], 4), words(sys.argv[2], 5)
    W6 = words(sys.argv[3], 6) if len(sys.argv) > 3 else {}
    ok = bad = amb = 0; checked = good6 = 0
    for lab, ks4 in sorted(W4.items(), key=str):
        if lab not in W5: continue
        if len(ks4) != 1 or len(W5[lab]) != 1: amb += 1; continue
        (lo4, hi4), (lo5, hi5) = ks4[0], W5[lab][0]
        il, ih = insertion(lo4, lo5), insertion(hi4, hi5)
        if not il or not ih: bad += 1; continue
        ok += 1
        if lab in W6 and len(W6[lab]) == 1:
            checked += 1
            i, j = il[0], ih[0]
            lo6 = lo4[:i] + lo5[i:i + 2] * 2 + lo4[i:]; hi6 = hi4[:j] + hi5[j:j + 2] * 2 + hi4[j:]
            good6 += (lo6, hi6) == W6[lab][0]
        print('%s|%s|%d|%s %s %s %s %s' % (lab[0] or '-', lab[1], lab[2], lab[3] or '-', lo4, hi4, lo5, hi5))
    print('labels in both limbs: %d with a 2-bit insertion, %d without, %d ambiguous; limb 6 predicted %d of %d' %
          (ok, bad, amb, good6, checked), file=sys.stderr)
