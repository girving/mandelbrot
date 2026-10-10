"""Quantitative tower test with exact sizes.  For each j = 0 decoration of the copy U (outer atom W(c) = A_c X(c), e
outer steps), predict its tower members j = 1..J as the roots of the quadratic model P_j(c') = W_U + W' c' (W_U, W' the
atom and its motion at the copy's centre), polish to centres of period (j+1) p + e, and compare the exact sizes
1/|A_W D_W| with the zeroth-order prediction
    |A_W D_W|_j ≈ |(f^e)'(x)|² |dc'/dc| |(g^j)'(c')| |P_j'(c') - W'|,   (g^j)'(c') = Π_{i=1..j} 2 P_i(c')
(outer factor and dc'/dc taken at the base), i.e. the kernel K(W, W') per outer atom.

  python3 tower_size.py [--cU ..] [--p 4] [--J 3]"""
import argparse, cmath, math
import numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('--cU', type=complex, default=complex(-1.9407998065, 0)); ap.add_argument('--p', type=int, default=4)
ap.add_argument('--emax', type=int, default=8); ap.add_argument('--rmin', type=float, default=3); ap.add_argument('--rmax', type=float, default=8)
ap.add_argument('--J', type=int, default=3); ap.add_argument('--bases', type=int, default=12); a = ap.parse_args()
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
def dorbit(c, n, z):
    d = 1 + 0j
    for _ in range(n): d *= 2 * z; z = z * z + c
    return d
def AD_exact(c, P):   # A_W D_W = Λ · ∂_c f^P(0)
    z, dz, lam = 0j, 0j, 1 + 0j
    for i in range(P):
        dz = 2 * z * dz + 1; z = z * z + c
        if i < P - 1: lam *= 2 * z
    return lam * dz
def A_of(c): return dorbit(c, p - 1, c)
def cprime(c): return A_of(c) * orbit(c, p)
def X_of(c, x0, e):
    x = x0
    for _ in range(60):
        st = orbit(c, e, x) / dorbit(c, e, x); x -= st
        if abs(st) < 1e-15 * (1 + abs(x)): return x
    return x if abs(st) < 1e-11 * (1 + abs(x)) else None
def X_cont(c0, x0, c1, e, steps=40):   # continue the atom branch from c0 to c1
    x = x0
    for t in range(1, steps + 1):
        x = X_of(c0 + (c1 - c0) * t / steps, x, e)
        if x is None: return None
    return x
def Ppoly(j):
    P = np.array([1.0, 0.0], dtype=complex)
    for _ in range(j): P = np.polyadd(np.polymul(P, P), np.array([1.0, 0.0], dtype=complex))
    return P   # P_j(c'): P_0 = c', P_{i+1} = P_i² + c'
cU = center(a.cU, p); h = 1e-9
dcp = (cprime(cU + h) - cprime(cU - h)) / (2 * h)
print('c_U %.12f%+.12fi p %d dc\'/dc %.4g%+.4gi' % (cU.real, cU.imag, p, dcp.real, dcp.imag))
bases = {}
for e in range(1, a.emax + 1):
    P = p + e
    for x in np.linspace(-a.rmax, a.rmax, 41):
        for y in np.linspace(-a.rmax, a.rmax, 41):
            if not (a.rmin <= abs(complex(x, y)) <= a.rmax): continue
            c0 = cU + complex(x, y) / dcp
            c = center(c0, P)
            if c is None or not exact_period(c, P) or abs(cprime(c) - complex(x, y)) > 1.5: continue
            if not (a.rmin - 0.5 <= abs(cprime(c))): continue
            bases[(P, round(c.real, 11), round(c.imag, 11))] = (c, e)
print('%d bases' % len(bases))
bl = sorted(bases.values(), key=lambda t: -1 / abs(AD_exact(t[0], p + t[1])))[:a.bases]
allres = {j: [] for j in range(1, a.J + 1)}
for c0, e in bl:
    x0 = orbit(c0, p)
    XU = X_cont(c0, x0, cU, e)
    if XU is None: continue
    WU = A_of(cU) * XU
    Wd = (A_of(cU + h) * X_of(cU + h, XU, e) - A_of(cU - h) * X_of(cU - h, XU, e)) / (2 * h) / dcp
    O0 = dorbit(c0, e, x0)
    s0 = 1 / abs(AD_exact(c0, p + e)); s0p = 1 / (abs(O0) ** 2 * abs(dcp) * abs(1 - Wd))
    print('base P %2d e %d  c\' %+.3f%+.3fi  W_U %+.3f%+.3fi  W\' %+.3f%+.3fi  size %.3e  pred/size %.3f' % (
        p + e, e, cprime(c0).real, cprime(c0).imag, WU.real, WU.imag, Wd.real, Wd.imag, s0, s0p / s0))
    for j in range(1, a.J + 1):
        Pj = Ppoly(j); poly = Pj.copy(); poly[-2] -= Wd; poly[-1] -= WU
        roots = np.roots(poly); Pjd = np.polyder(Pj)
        P = (j + 1) * p + e
        for r in roots:
            # map c' -> c by Newton on the exact c'(c), then centre
            c = cU + r / dcp
            for _ in range(20):
                d = (cprime(c + h) - cprime(c - h)) / (2 * h); c -= (cprime(c) - r) / d
            cc = center(c, P)
            if cc is None or not exact_period(cc, P): allres[j].append(None); continue
            # membership: the orbit after (j+1) returns sits on the continued atom branch
            xj = orbit(cc, (j + 1) * p); Xj = X_cont(c0, x0, cc, e)
            ok = Xj is not None and abs(xj - Xj) < 1e-6 * (1 + abs(Xj)) and abs(cprime(cc) - r) < 0.05 * max(1, abs(r))
            if not ok: allres[j].append(None); continue
            gj = np.prod([2 * np.polyval(Ppoly(i), r) for i in range(j)])   # (g^j)'(c') = Π_{i=1..j} 2 g^i(0) = Π 2 P_{i-1}
            pred = 1 / (abs(O0) ** 2 * abs(dcp) * abs(gj) * abs(np.polyval(Pjd, r) - Wd))
            act = 1 / abs(AD_exact(cc, P))
            allres[j].append((act, pred, math.log(pred / act), abs(r)))
for j in range(1, a.J + 1):
    v = [t for t in allres[j] if t]; n = len(allres[j])
    if not v: print('j %d: none of %d found' % (j, n)); continue
    L = sorted(abs(t[2]) for t in v)
    tot_a = sum(t[0] ** 2 for t in v); tot_p = sum(t[1] ** 2 for t in v)
    print('j %d: %d/%d members found; |log pred/actual| median %.3f max %.3f; mass (Σ size²) pred/actual %.4f' % (
        j, len(v), n, L[len(L) // 2], L[-1], tot_p / tot_a))
