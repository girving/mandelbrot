"""Primitive hyperbolic components of M of a given period (centers), for the Lavaurs-model tuning filter.

Centers of period p: roots of f_c^p(0) = 0 not of lower period (Newton from numpy's roots of the Gleason polynomial,
refined).  A component is a satellite iff at its root (multiplier 1 on the period-p cycle) some cycle of period d | p,
d < p, is indifferent (|multiplier| = 1); otherwise primitive (a cusp).

  python3 primitives.py 5 6   (prints "p c_re c_im" for the primitive centers)"""
import sys, cmath
import numpy as np

def fp(c, p):
    z, dz = 0j, 0j
    for _ in range(p): dz = 2 * z * dz + 1; z = z * z + c
    return z, dz

def centers(p):
    # polynomial in c: f^p(0) as coefficients
    poly = np.array([0.0])  # z_0 = 0
    zpoly = np.poly1d([0.0])
    cpoly = np.poly1d([1.0, 0.0])
    for _ in range(p): zpoly = zpoly * zpoly + cpoly
    roots = np.roots(zpoly.coeffs)
    out = []
    for c in roots:
        c = complex(c)
        for _ in range(50):
            z, dz = fp(c, p)
            c -= z / dz
        # exact period: smallest d with f^d(0) = 0
        per = next(d for d in range(1, p + 1) if abs(fp(c, d)[0]) < 1e-9)
        if per == p and not any(abs(c - o) < 1e-9 for o in out): out.append(c)
    return out

def cycle_mult(c, d, z0):
    """A period-d cycle near z0 (Newton on f^d(z) = z) and its multiplier"""
    z = z0
    for _ in range(100):
        w, dw = z, 1
        for _ in range(d): dw = 2 * w * dw; w = w * w + c
        step = (w - z) / (dw - 1); z -= step
        if abs(step) < 1e-14: break
    w, dw = z, 1
    for _ in range(d): dw = 2 * w * dw; w = w * w + c
    return z, dw

def root(c0, p):
    """The root of the component with center c0: solve f^p(z) = z, (f^p)'(z) = 1 in (z, c)"""
    c, z = c0, 0j
    # move along the multiplier map μ = 0 → 1 in steps
    for k in range(1, 41):
        mu = k / 40
        for _ in range(40):
            # F(z, c) = (f^p(z) - z, (f^p)'(z) - mu)
            w, dwz, dwc, d2 = z, 1, 0, 0
            for _ in range(p):
                d2 = 2 * (dwz * dwz + w * d2); dwc = 2 * w * dwc + 1; dwz = 2 * w * dwz; w = w * w + c
            # derivative of dwz wrt c: approximate numerically
            h = 1e-7
            w2, dz2 = z, 1
            for _ in range(p): dz2 = 2 * w2 * dz2; w2 = w2 * w2 + (c + h)
            dmu_dc = (dz2 - dwz) / h
            J = np.array([[dwz - 1, dwc], [d2, dmu_dc]], dtype=complex)
            r = np.array([w - z, dwz - mu], dtype=complex)
            dz_, dc_ = np.linalg.solve(J, r)
            z -= dz_; c -= dc_
            if abs(dz_) + abs(dc_) < 1e-14: break
    return c, z

def primitive(c0, p):
    cr, zr = root(c0, p)
    for d in range(1, p):
        if p % d: continue
        # cycles of period d at the root: start Newton from the orbit points of the period-p cycle
        w = zr
        for _ in range(p):
            zc, m = cycle_mult(cr, d, w)
            if abs(abs(m) - 1) < 1e-6: return False
            w = w * w + cr
    return True

if __name__ == '__main__':
    for p in map(int, sys.argv[1:]):
        cs = centers(p)
        prims = [c for c in cs if primitive(c, p)]
        print('# period %d: %d centers, %d primitive' % (p, len(cs), len(prims)), file=sys.stderr)
        for c in prims: print('%d %.15f %.15f' % (p, c.real, c.imag))
