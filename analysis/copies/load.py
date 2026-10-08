"""Hyperbolic components through period 16 with angles and tuning (tuning's TUNING_DUMP), and tuning decomposition.

Roots: index, period p, lo and hi angle words (p bits), satellite flag, maximal tuning (root index or -1), center,
area.  A component with maximal tuning W0 (period d) is W0 ⋆ V: decoding W's lo word in d-bit blocks (W0's lo
word → 0, hi word → 1) gives V's lo word."""
import os, math
A_CARD = 3 * math.pi / 8

class Root:
    __slots__ = ('i', 'p', 'lo', 'hi', 'sat', 'parent', 'c', 'area', 'v')

def load(path):
    roots = []
    for line in open(path):
        t = line.split()
        r = Root()
        r.i, r.p, r.lo, r.hi, r.sat, r.parent = int(t[0]), int(t[1]), int(t[2]), int(t[3]), int(t[4]) == 1, int(t[5])
        r.c = complex(float(t[6]), float(t[7])); r.area = float(t[8]); r.v = None
        roots.append(r)
    by_lo = {(r.p, r.lo): r for r in roots}
    for r in roots:
        if r.parent < 0: continue
        w0 = roots[r.parent]; d = w0.p; m = r.p // d
        bits = 0
        for j in range(m):
            blk = (r.lo >> (d * (m - 1 - j))) & ((1 << d) - 1)
            bit = 0 if blk == w0.lo else 1 if blk == w0.hi else None
            assert bit is not None, (r.i, blk, w0.lo, w0.hi)
            bits = 2 * bits + bit
        r.v = by_lo.get((m, bits))
        if r.v is None and m == 1: r.v = None   # W = W0 ⋆ cardioid: not a proper tuning
    return roots

def nrp(r):
    """Non-renormalizable primitive component"""
    return r.parent < 0 and not r.sat
