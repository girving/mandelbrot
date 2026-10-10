"""Do a small copy's decorations come in Böttcher towers?  In M, around a primitive copy U*M (period p, centre c_U):
the renormalized parameter c' = A_c f_c^p(0) (A_c = (f_c^p)''(0)/2: the normalized critical value of the return map)
maps the copy to ≈ M; a decoration at j returns has exit point w = Φ_M(c')^(2^j), so from each j = 0 decoration
(period p + e) the theory predicts j = 1 decorations (period 2p + e) at c' = Ψ_M(±√Φ_M(c'_0)) with size ratio
|s_1|²/|s_0|² = |dc'/dc|_1² |Φ_M'(c'_0)|² / (|dc'/dc|_0² |2 Φ_M(c'_1) Φ_M'(c'_1)|²).

  python3 copy_towers.py [--cU -1.9407998065] [--p 4] [--emax 6]"""
import argparse, cmath, math
import numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('--cU', type=complex, default=complex(-1.9407998065, 0)); ap.add_argument('--p', type=int, default=4)
ap.add_argument('--emax', type=int, default=6); ap.add_argument('--rmin', type=float, default=1.3); ap.add_argument('--rmax', type=float, default=5.0)
a = ap.parse_args()
def center(c, P):
    for _ in range(80):
        z, dz = 0j, 0j
        for _ in range(P): dz = 2 * z * dz + 1; z = z * z + c
        st = z / dz; c -= st
        if abs(st) < 1e-15 * (1 + abs(c)): return c
        if not abs(st) < 1: return None
    return None
def exact_period(c, P):
    z = 0j
    for k in range(1, P):
        z = z * z + c
        if abs(z) < 1e-8 and P % k == 0: return False
    return True
def size2(c, P):   # |s|² with s = 1/(β Λ²), Λ = Π 2 z_i, β = Σ 1/Π_{j≤i} 2 z_j
    z = c; lam = 1; beta = 0
    for i in range(1, P):
        lam *= 2 * z; beta += 1 / lam; z = z * z + c
    return abs(1 / (beta * lam * lam)) ** 2
def renorm(c, p):  # c' = A_c f^p(0), and dc'/dc
    def F(c):
        z, a1, a2 = 0j, 1 + 0j, 0j   # z, dz/dz0, d²z/dz0² at z0 = 0
        for _ in range(p): a2 = 2 * (a1 * a1 + z * a2); a1 = 2 * z * a1; z = z * z + c
        return (a2 / 2) * z
    h = 1e-7 * (1 + abs(c))
    return F(c), (F(c + h) - F(c - h)) / (2 * h)
def phi_M(c, prev):
    z = c
    for n in range(1, 3000):
        z = z * z + c
        if abs(z) > 1e12: break
    else: return None
    L = cmath.log(z); N = 2 ** n
    k = round((N * cmath.phase(prev) - L.imag) / (2 * math.pi))
    return cmath.exp((L + 2j * math.pi * k) / N)
def Phi(c, steps=300):
    far = c * (100 / abs(c)); prev = far
    for t in range(1, steps + 1):
        cc = c + (far - c) * (1 - t / steps) ** 3
        prev = phi_M(cc, prev)
        if prev is None: return None
    return prev
def dPhi(c):
    h = 1e-6 * abs(c); p0 = Phi(c)
    return (Phi(c + h) - Phi(c - h)) / (2 * h) if p0 is not None else None
def Psi(zeta, steps=60):   # Φ_M(c) = ζ by Newton, continued from |ζ| = 100 along the ray
    far = zeta * (100 / abs(zeta)); c = far
    for t in range(1, steps + 1):
        zt = zeta + (far - zeta) * (1 - t / steps) ** 3
        for _ in range(30):
            f = Phi(c); 
            if f is None: return None
            d = dPhi(c); st = (f - zt) / d; c -= st
            if abs(st) < 1e-13 * abs(c): break
    return c
cU = center(a.cU, a.p); p = a.p
A_D = renorm(cU + 1e-9, p)[1]   # dc'/dc at c_U = A D
print('c_U %.12f%+.12fi, dc\'/dc = A D = %.4g%+.4gi (|AD| %.4g)' % (cU.real, cU.imag, A_D.real, A_D.imag, abs(A_D)))
# j = 0 decorations: Newton from a ring of starts c' ∈ annulus [rmin, rmax], periods p + e
found = {}
for e in range(1, a.emax + 1):
    P = p + e
    for x in np.linspace(-a.rmax, a.rmax, 41):
        for y in np.linspace(-a.rmax, a.rmax, 41):
            if not (a.rmin <= abs(complex(x, y)) <= a.rmax): continue
            c0 = cU + complex(x, y) / A_D
            c = center(c0, P)
            if c is None or not exact_period(c, P) or abs(c - c0) * abs(A_D) > 1.5: continue
            cp, _ = renorm(c, p)
            if not (a.rmin * 0.8 < abs(cp) < a.rmax * 1.5): continue
            k = (P, round(c.real, 10), round(c.imag, 10))
            found[k] = (c, P, cp)
print('%d j = 0 decorations (periods %d..%d): %s' % (len(found), p + 1, p + a.emax, sorted({k[0] for k in found})))
rows = sorted(found.values(), key=lambda t: -size2(t[0], t[1]))[:12]
for c, P, cp in rows:
    z0 = Phi(cp)
    if z0 is None: continue
    s0 = size2(c, P); d0 = renorm(c, p)[1]; dz0 = dPhi(cp)
    print('j=0  P %2d c\' %+.4f%+.4fi |s|² %.3e  ζ %.4f∠%.4f' % (P, cp.real, cp.imag, s0, abs(z0), (cmath.phase(z0) / (2 * math.pi)) % 1))
    # j returns then the atom: P_j(c') = c'_0 with P_0(c') = c', P_{i+1} = P_i² + c' (the critical orbit of z² + c')
    for j in (1, 2):
        coef = np.array([1.0 + 0j])        # polynomial P_j(c') coefficients, highest first
        Pj = np.array([1.0, 0.0], dtype=complex)  # c'
        for _ in range(j): Pj = np.polyadd(np.polymul(Pj, Pj), np.array([1.0, 0.0], dtype=complex))
        roots = np.roots(np.polysub(Pj, np.array([cp], dtype=complex)))
        dPj = np.polyder(Pj)
        for cp1 in roots:
            c1 = center(cU + cp1 / A_D, P + j * p)
            if c1 is None: print('      j=%d predicted c\' %+.4f%+.4fi: Newton failed' % (j, cp1.real, cp1.imag)); continue
            cp1f, d1 = renorm(c1, p)
            # size ∝ 1/(A_W D_W): the transversality P_j'(c') and the dynamic derivative (g^j)'(c') = Π_{i<j} 2 g^i(c')
            gd, g = 1.0 + 0j, cp1
            for _ in range(j): gd *= 2 * g; g = g * g + cp1
            pred = abs(d0) ** 2 / (abs(d1) ** 2 * abs(np.polyval(dPj, cp1)) ** 2 * abs(gd) ** 2)
            print('      j=%d predicted c\' %+.4f%+.4fi → found %+.4f%+.4fi (|Δ| %.2e, exact %s)  |s|² ratio %.4g predicted %.4g' % (
                j, cp1.real, cp1.imag, cp1f.real, cp1f.imag, abs(cp1f - cp1), exact_period(c1, P + j * p), size2(c1, P + j * p) / s0, pred))
