"""Areas of Lavaurs-model components (prototype, double precision).

The cycle through the critical point w = 1/2: R_σ(w) = F^n(g_σ(F(w))), g_σ = Ψ₊/₋(Φ_a(·) + σ).  The component is
{σ : |λ(σ)| < 1} for the multiplier λ of R_σ's fixed point near 1/2; its boundary σ(μ), |μ| = 1, comes from Newton
on (R_σ(w) - w, R_σ'(w) - μ) with continuation in μ, and area = π Σ k |a_k|² for σ(μ) - σ(0) = Σ a_k μ^k (DFT).
All maps carry 2-jets (value, ∂w, ∂σ, ∂ww, ∂wσ).  C_w = (π²/4) area (σ_model = -σ_M/2 + 3πi/8, δ = iπ/(2k + σ_M))."""
import cmath, math, sys
import numpy as np
from lavaurs import A_COEF, C_LOG, R_LOC, F, dF, _L

def phi_series3(w, attracting):
    """Φ, Φ', Φ'' of the series"""
    s = 0.25 / (w * w) + 0.25 / w + C_LOG * _L(w, attracting)
    d = -0.5 / w**3 - 0.25 / (w * w) + C_LOG / w
    dd = 1.5 / w**4 + 0.5 / w**3 - C_LOG / (w * w)
    p = 1.0; pm = 0.0
    for j, a in enumerate(A_COEF, 1):
        if j >= 2: dd += j * (j - 1) * a * pm
        d += j * a * p; pm = p; p *= w; s += a * p
    return s, d, dd

def iterate3(w, n):
    """F^n and its first two derivatives"""
    d1, d2 = 1.0 + 0j, 0j
    for _ in range(n):
        d2 = 2 * d1 * d1 + (-1 + 2 * w) * d2
        d1 = (-1 + 2 * w) * d1
        w = F(w)
    return w, d1, d2

def phi_a3(w):
    """Φ_a, Φ_a', Φ_a'', entering petal"""
    d1, d2 = 1.0 + 0j, 0j
    for n in range(100000):
        if abs(w) < R_LOC and abs(w.imag) < 0.7 * abs(w.real):
            s, d, dd = phi_series3(w, True)
            return s - n / 2, d * d1, dd * d1 * d1 + d * d2, (1 if w.real > 0 else -1) * (1 if n % 2 == 0 else -1)
        d2 = 2 * d1 * d1 + (-1 + 2 * w) * d2
        d1 = (-1 + 2 * w) * d1
        w = F(w)
        if abs(w) > 10: raise ValueError('escaped')
    raise ValueError('no petal')

def psi3(zeta, petal):
    """Ψ_petal and its first two derivatives"""
    if petal < 0:
        w, d1, d2 = psi3(zeta - 0.5, 1)
        return F(w), dF(w) * d1, 2 * d1 * d1 + dF(w) * d2
    m = max(0, math.ceil(zeta.real + max(100.0, 2 * abs(zeta.imag))))
    if m > 20000: raise ValueError('Ψ argument too far out (Re ζ = %.3g)' % zeta.real)
    zl = zeta - m
    w = 1j / (2 * cmath.sqrt(-zl))
    for _ in range(40):
        s, d, dd = phi_series3(w, False)
        st = (s - zl) / d; w -= st
        if abs(st) < 1e-17 * abs(w): break
    s, d, dd = phi_series3(w, False)
    w1 = 1 / d; w2 = -dd * w1**3
    v, e1, e2 = iterate3(w, 2 * m)
    return v, e1 * w1, e2 * w1 * w1 + e1 * w2

class Jet:
    """(value, ∂w, ∂σ, ∂ww, ∂wσ)"""
    __slots__ = ('v', 'w', 's', 'ww', 'ws')
    def __init__(s, v, w, sg, ww, ws): s.v, s.w, s.s, s.ww, s.ws = v, w, sg, ww, ws
    def apply(s, f0, f1, f2):  # f(s) for scalar f with derivatives f1, f2 at s.v
        return Jet(f0, f1 * s.w, f1 * s.s, f2 * s.w * s.w + f1 * s.ww, f2 * s.w * s.s + f1 * s.ws)

def R(w, sigma, n):
    x = Jet(w, 1, 0, 0, 0)
    x = x.apply(F(x.v), dF(x.v), 2)                               # F
    p0, p1, p2, petal = phi_a3(x.v)
    x = x.apply(p0, p1, p2)                                       # Φ_a
    x = Jet(x.v + sigma, x.w, x.s + 1, x.ww, x.ws)                # + σ
    q0, q1, q2 = psi3(x.v, petal)
    x = x.apply(q0, q1, q2)                                       # Ψ
    for _ in range(n): x = x.apply(F(x.v), dF(x.v), 2)            # F^n
    return x

def newton(w, sigma, n, mu, iters=30):
    for _ in range(iters):
        x = R(w, sigma, n)
        F1, F2 = x.v - w, x.w - mu
        a, b, d, e = x.w - 1, x.s, x.ww, x.ws
        det = a * e - b * d
        dw = (F1 * e - b * F2) / det; ds = (a * F2 - d * F1) / det
        w -= dw; sigma -= ds
        if not abs(dw) + abs(ds) < 1: raise ValueError('Newton diverged')
        if abs(dw) + abs(ds) < 1e-15: break
    return w, sigma, abs(dw) + abs(ds)

def area(sigma0, n, N=64, radial=8):
    w, s, err = newton(0.5 + 0j, sigma0, n, 0)
    center = s
    for i in range(1, radial + 1):                               # out to μ = e^{iπ/N} radially
        w, s, err = newton(w, s, n, (i / radial) * cmath.exp(1j * math.pi / N))
    pts = []
    for j in range(N):
        mu = cmath.exp(1j * math.pi * (2 * j + 1) / N)
        for t in range(1, 5):                                     # small steps along the circle
            nu = cmath.exp(1j * math.pi * (2 * j - 1 + 2 * t / 4) / N) if j else mu
            w, s, err = newton(w, s, n, nu)
        pts.append(s - center)
    pts = np.array(pts)
    mus = np.exp(1j * math.pi * (2 * np.arange(N) + 1) / N)
    a = np.array([np.mean(pts * mus**(-k)) for k in range(N)])
    A = math.pi * sum(k * abs(a[k])**2 for k in range(1, N))
    half = pts[::2]; mh = mus[::2]
    ah = np.array([np.mean(half * mh**(-k)) for k in range(N // 2)])
    A2 = math.pi * sum(k * abs(ah[k])**2 for k in range(1, N // 2))
    return A, center, abs(A - A2) / A

if __name__ == '__main__':
    # S1: σ_M = 0.6408400903+1.0869901418i  ->  σ_model = -σ_M/2 + 3πi/8 (mod 1/2), with some n
    sM = complex(sys.argv[1]) if len(sys.argv) > 1 else 0.6408400903 + 1.0869901418j
    target = float(sys.argv[2]) if len(sys.argv) > 2 else 1.724974063498964762e-3
    base = -sM / 2 + 3j * math.pi / 8
    for n in range(0, 6):
        try:
            A, c, conv = area(base - n * 0.5, n)
            print('n=%d center %.10f%+.10fi  area %.15e  (π²/4) area = %.15e  vs target: %.2e  (conv %.0e)' %
                  (n, c.real, c.imag, A, math.pi**2 / 4 * A, math.pi**2 / 4 * A / target - 1, conv))
        except (ValueError, ZeroDivisionError, OverflowError) as e:
            print('n=%d: %s' % (n, e))
