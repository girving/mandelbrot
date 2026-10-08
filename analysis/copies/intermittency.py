"""Intermittency families of primitive copies at the airplane's cusp c = -1.75 (period-3 saddle-node), c > -1.75.

Family A: kneading L (RLL)^k C (period 3k + 2); family B: L (RLL)^k RL C (period 3k + 4).  Centers by Newton on
f_c^p(0) = 0 from an extrapolation d = c + 1.75 ≈ A/(k + B)^2 of the previous members, accepted only if the
itinerary is exactly the family's word.  Writes bulb_batch jobs (P = 0: the component itself)."""
import sys

def itinerary(c, p):
    z, s = 0.0, ''
    for i in range(p):
        z = z * z + c
        s += 'C' if i == p - 1 else ('R' if z > 0 else 'L')
    return s

def center(c, p):
    for _ in range(100):
        z, dz = 0.0, 0.0
        for i in range(p): dz = 2 * z * dz + 1; z = z * z + c
        step = z / dz; c -= step
        if abs(step) < 1e-16 * abs(c): break
    return c

def family(word, kmax):
    out = {}
    seeds = {'A': [(2, -1.7110794700131522), (3, -1.73200627287), (4, -1.739717601451)],
             'B': [(2, -1.721915099589), (3, -1.73574321564), (4, -1.741394537769585)]}[word]
    pat = (lambda k: 'L' + 'RLL' * k + 'C') if word == 'A' else (lambda k: 'L' + 'RLL' * k + 'RLC')
    for k, c in seeds:
        p = len(pat(k)); out[k] = (center(c, p), p)
    for k in range(5, kmax + 1):
        (k1, (c1, _)), (k2, (c2, _)) = sorted(out.items())[-2:]
        # d = A/(k + B)^2: 1/sqrt(d) is linear in k
        s1, s2 = (c1 + 1.75) ** -0.5, (c2 + 1.75) ** -0.5
        s = s2 + (s2 - s1) * (k - k2) / (k2 - k1)
        p = len(pat(k)); c = center(-1.75 + s ** -2, p)
        if itinerary(c, p) != pat(k):
            print('family %s: k=%d Newton left the family' % (word, k), file=sys.stderr); break
        out[k] = (c, p)
    return out

if __name__ == '__main__':
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 120
    for word in 'AB':
        for k, (c, p) in sorted(family(word, kmax).items()):
            print('%s%d 0 %.17g 0 0 %d' % (word, k, c, p))
