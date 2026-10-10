"""Exact size of a copy decoration: for the copy U (period p) and a decoration W of period P = (j+1) p + e (orbit leaves
the quadratic-like domain after j + 1 returns at x = R^(j+1)(0), R = f^p, then f^e(x) = 0):
    A_W D_W = (f^e)'(x)² · (g^j)'(c') · Q'(c),   Q(c) = A_c R_c^(j+1)(0) - A_c X(c),
with A_c = (f_c^(p-1))'(c), g the renormalized return map, c' = A_c f_c^p(0), X(c) the preimage branch of 0 under f_c^e
through x.  Checks the identity against Λ² β (the exact normal-form product) at actual centres, and splits Q' into the
quadratic-model part (dc'/dc) P_j'(c') and the rest (the atom's motion and the nonlinearity).

  python3 decor_size.py [--cU ..] [--p 4]"""
import argparse, cmath, math
import numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('--cU', type=complex, default=complex(-1.9407998065, 0)); ap.add_argument('--p', type=int, default=4)
ap.add_argument('--emax', type=int, default=12); ap.add_argument('--rmin', type=float, default=3); ap.add_argument('--rmax', type=float, default=10)
ap.add_argument('--top', type=int, default=8); a = ap.parse_args()
p = a.p
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
def orbit(c, n, z=0j):
    for _ in range(n): z = z * z + c
    return z
def dorbit(c, n, z):   # (f_c^n)'(z)
    d = 1 + 0j
    for _ in range(n): d *= 2 * z; z = z * z + c
    return d
def lam_beta(c, P):   # Λ = Π_{i=1}^{P-1} 2 z_i, β = Σ 1/Π_{j≤i} 2 z_j: A_W D_W = Λ² β
    z = c; lam = 1; beta = 0
    for i in range(1, P):
        lam *= 2 * z; beta += 1 / lam; z = z * z + c
    return lam, beta
def A_of(c): return dorbit(c, p - 1, c)
def X_of(c, x0, e):   # preimage of 0 under f_c^e through x0 (Newton)
    x = x0
    for _ in range(50):
        f, d = orbit(c, e, x), dorbit(c, e, x)
        st = f / d; x -= st
        if abs(st) < 1e-15 * (1 + abs(x)): break
    return x
cU = center(a.cU, p)
h = 1e-9
AD = (A_of(cU + h) * orbit(cU + h, p) - A_of(cU - h) * orbit(cU - h, p)) / (2 * h)
print('c_U %.12f%+.12fi, dc\'/dc = %.4g%+.4gi' % (cU.real, cU.imag, AD.real, AD.imag))
# decorations of every depth near the copy (periods p + 1 .. p + emax·?), then classify (j, e) by the orbit
found = {}
for P in range(p + 1, p + a.emax + 1):
    for x in np.linspace(-a.rmax, a.rmax, 31):
        for y in np.linspace(-a.rmax, a.rmax, 31):
            if not (1.0 <= abs(complex(x, y)) <= a.rmax): continue
            c0 = cU + complex(x, y) / AD
            c = center(c0, P)
            if c is None or not exact_period(c, P) or abs(c - c0) * abs(AD) > 1.5: continue
            found[(P, round(c.real, 11), round(c.imag, 11))] = (c, P)
print('%d decorations' % len(found))
rows = []
for c, P in found.values():
    A = A_of(c); cp = A * orbit(c, p)
    # j: the number of returns after the first while the renormalized orbit stays within |w| < 2.5
    w, j = cp, 0
    while abs(w) < 2.5 and (j + 2) * p <= P:
        j += 1; w = A * orbit(c, (j + 1) * p)
    e = P - (j + 1) * p
    if e < 1: continue
    x = orbit(c, (j + 1) * p)
    Odr = dorbit(c, e, x)                       # (f^e)'(x)
    gj = dorbit(c, j * p, orbit(c, p))          # (g^j)'(c') = (R^j)'(R(0))
    Q = lambda cc: A_of(cc) * (orbit(cc, (j + 1) * p) - X_of(cc, x, e))
    Qd = (Q(c + h) - Q(c - h)) / (2 * h)
    lam, beta = lam_beta(c, P)
    # exact normal form: A_W = Λ = (f^(P-1))'(c), D_W = ∂_c f_c^P(0) = Λ (1 + β)
    Dexact = (orbit(c + h, P) - orbit(c - h, P)) / (2 * h)
    ident = (Odr ** 2 * gj * Qd) / (lam * Dexact)
    # the quadratic model's part of Q': (dc'/dc) P_j'(c') with P_j the critical orbit polynomial of z² + c'
    Pj = np.array([1.0, 0.0], dtype=complex)
    for _ in range(j): Pj = np.polyadd(np.polymul(Pj, Pj), np.array([1.0, 0.0], dtype=complex))
    dcp = (A_of(c + h) * orbit(c + h, p) - A_of(c - h) * orbit(c - h, p)) / (2 * h)
    quad = dcp * np.polyval(np.polyder(Pj), cp)
    Wd = (A_of(c + h) * X_of(c + h, x, e) - A_of(c - h) * X_of(c - h, x, e)) / (2 * h) / dcp   # dW/dc'
    rows.append((abs(1 / (lam * Dexact)), P, j, e, cp, ident, Qd / quad, Wd, abs(Dexact / (lam * (1 + beta)))))
rows.sort(key=lambda r: -r[0])
print('1/|A_W D_W|  P  j  e   c\'                  identity O² (g^j)\' Q\'/(A_W D_W)   Q\'/quadratic part   dW/dc\'     |D/(Λ(1+β))|')
for s, P, j, e, cp, ident, qr, Wd, dchk in rows[:40]:
    print('%.3e %3d %2d %2d  %+8.3f%+8.3fi   %+.6f%+.6fi    %+.3f%+.3fi    %+.3f%+.3fi   %.6f' % (s, P, j, e, cp.real, cp.imag, ident.real, ident.imag, qr.real, qr.imag, Wd.real, Wd.imag, dchk))

med = lambda v: sorted(v)[len(v) // 2] if v else float('nan')
r0 = [r for r in rows if r[2] == 0]; r1 = [r for r in rows if r[2] >= 1]
print('p %d |dc\'/dc| %.4g: %d decorations (j = 0: %d, j ≥ 1: %d); identity max |ratio - 1| %.1e' % (p, abs(AD), len(rows), len(r0), len(r1), max(abs(r[5] - 1) for r in rows)))
print('   median |dW/dc\'| (j = 0) %.4f;  median |W\'/c\'| %.4f' % (med([abs(r[7]) for r in r0]), med([abs(r[7]) / abs(r[4]) for r in r0])))
# j ≥ 1: the quadratic model's residual beyond the atom term: Q'/((dc'/dc) P_j') - (1 - W'/P_j')
res = []
for r in r1:
    P, j, cp = r[1], r[2], r[4]
    Pj = np.array([1.0, 0.0], dtype=complex)
    for _ in range(j): Pj = np.polyadd(np.polymul(Pj, Pj), np.array([1.0, 0.0], dtype=complex))
    pd = np.polyval(np.polyder(Pj), cp)
    res.append(abs(r[6] - (1 - r[7] / pd)))
print('   j ≥ 1: median nonlinearity residual |Q\'/quad - (1 - W\'/P_j\')| %.4f' % med(res))
