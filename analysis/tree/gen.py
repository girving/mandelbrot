"""Data for the extended (digits + separator) Hankel block of the satellite tree.

Nodes are addresses (r1, r2, ...) of rotation numbers; F(address) is bulb_areas' normalized area of the last bulb
relative to its parent's multiplier map.  Stage A: children (q ≤ 20) of every root-level word u (q_u ≤ 20).  Stage B:
children u·s (q_u ≤ 8, q_s ≤ 40) of each r1 (q1 ≤ 8, r1 ≤ 1/2), plus the nodes (r1, u).  Stage C: children (q ≤ 20)
of each node (r1, u).  This script writes run lists; run.sh runs bulb_areas per parent.
  python3 gen.py A|B|C"""
import os, sys, pickle
from fractions import Fraction
from math import gcd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/../wfa')
D = os.environ['BULB_DATA'] + '/tree/'
def val(w):
    x = Fraction(0)
    for a in reversed(w): x = 1 / (a + x)
    return x
def cont(w):
    q0, q1 = 0, 1
    for a in w: q0, q1 = q1, a * q1 + q0
    return q1
def words(maxq, canonical):
    out = []
    def rec(w):
        if w and (not canonical or w[-1] >= 2): out.append(tuple(w))
        a = 1
        while cont(w + [a]) <= maxq:
            rec(w + [a]); a += 1
    rec([])
    return out
def rationals(maxq):
    return [Fraction(p, q) for q in range(2, maxq + 1) for p in range(1, q) if gcd(p, q) == 1]
U20 = [()] + words(20, False); U8 = [()] + words(8, False); S40 = words(40, True); S20 = rationals(20)
R1 = [x for x in rationals(8) if x <= Fraction(1, 2)]
def center_table(path):
    """{child rational: (center_re, center_im)} from a bulb_areas output"""
    out = {}
    for line in open(path):
        t = line.split()
        if t[2] != 'failed': out[Fraction(int(t[0]), int(t[1]))] = (t[2], t[3])
    return out
if __name__ == '__main__':
    stage = sys.argv[1]
    os.makedirs(D + 'ext', exist_ok=True)
    jobs = []   # (name, period, center_re, center_im, children)
    if stage == 'A':
        card = {}
        for n in ('all_e.out', 'stage2.out'):
            card.update(center_table(os.environ['BULB_DATA'] + '/' + n))
        vals = sorted({val(list(u)) for u in U20 if u} - {Fraction(1)})
        for v in vals:
            c = card.get(v) or card.get(1 - v)
            if v not in card:  # mirror: conjugate center
                c = (c[0], repr(-float(c[1])))
            jobs.append(('A_%d-%d' % (v.numerator, v.denominator), v.denominator, c[0], c[1], S20))
    elif stage == 'B':
        card = center_table(D + 'cardioid.out')
        for r1 in R1:
            kids = sorted({val(list(u) + list(s)) for u in U8 for s in S40} | ({val(list(u)) for u in U8 if u} - {Fraction(1)}))
            c = card[r1]
            jobs.append(('B_%d-%d' % (r1.numerator, r1.denominator), r1.denominator, c[0], c[1], kids))
    elif stage == 'C':
        for r1 in R1:
            ct = center_table(D + 'ext/B_%d-%d.out' % (r1.numerator, r1.denominator))
            for v in sorted({val(list(u)) for u in U8 if u} - {Fraction(1)}):
                c = ct[v]
                jobs.append(('C_%d-%d_%d-%d' % (r1.numerator, r1.denominator, v.numerator, v.denominator),
                             r1.denominator * v.denominator, c[0], c[1], S20))
    with open(D + 'ext/jobs_%s.txt' % stage, 'w') as f:
        for name, P, cr, ci, kids in jobs:
            with open(D + 'ext/%s.in' % name, 'w') as g:
                for x in kids: print(x.numerator, x.denominator, file=g)
            print(name, P, cr, ci, file=f)
    print(stage, len(jobs), 'parents,', sum(len(j[4]) for j in jobs), 'children')
