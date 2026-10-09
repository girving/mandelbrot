"""Lavaurs model at the 1/2 root (prototype, double precision).

w = z + 1/2, F(w) = -w + w² (f(z) = z² - 3/4); critical point w = 1/2, critical value v = -1/4.  Fatou coordinate Φ with
Φ(F(w)) = Φ(w) + 1/2 (fatou_series.py): attracting petals around the real axis, repelling ones around the imaginary
axis; L(w) = log w where Re w > 0 (attracting) or Im w > 0 (repelling), log(-w) on the partner petal.
  phi_a(w): F-iterate into an attracting petal, then the series, minus half the step count (and its derivative)
  psi(ζ, petal): repelling Fatou parametrizations.  F swaps the two repelling petals Q±, so there is one per petal:
          Ψ₊(ζ) = local inverse in Q+ (upper) at ζ - m (Re ζ - m ≪ 0) followed by F^{2m}, so Ψ₊(ζ + 1) = F²(Ψ₊(ζ)),
          and Ψ₋(ζ) = F(Ψ₊(ζ - 1/2)).
  transit(w, σ, cross): the Lavaurs map g_σ(w) = Ψ_s(Φ_a(w) + σ), the exit petal s = sign(Re w) of the entering
          petal P± (cross = False) or the opposite one (cross = True); both commute with F, the data decides."""
import cmath, math
from fatou_series import series

_A, _B, _C, _a = series(34)
A_COEF = [float(x) for x in _a[:30]]
C_LOG = float(_C)

def F(w): return -w + w * w
def dF(w): return -1 + 2 * w

def _L(w, attracting):
    if attracting: return cmath.log(w) if w.real > 0 else cmath.log(-w)
    return cmath.log(w) if w.imag > 0 else cmath.log(-w)

def phi_series(w, attracting):
    s = 0.25 / (w * w) + 0.25 / w + C_LOG * _L(w, attracting)
    d = -0.5 / w**3 - 0.25 / (w * w) + C_LOG / w
    p = 1.0
    for j, a in enumerate(A_COEF, 1):
        d += j * a * p; p *= w; s += a * p
    return s, d

R_LOC = 0.05

def phi_a(w, max_steps=100000):
    """Attracting Fatou coordinate and its derivative; None if the orbit does not enter a petal"""
    dw = 1.0 + 0j
    for n in range(max_steps):
        if abs(w) < R_LOC and abs(w.imag) < 0.7 * abs(w.real):
            s, d = phi_series(w, True)
            return s - n / 2, d * dw
        dw *= dF(w); w = F(w)
        if abs(w) > 10: return None
    return None

def psi_plus(zeta):
    """Ψ₊(ζ) and Ψ₊'(ζ)"""
    # local region: Re ζ' ≤ -max(400, 2|Im ζ|)  (|w| ≲ 0.025)
    m = max(0, math.ceil(zeta.real + max(400.0, 2 * abs(zeta.imag))))
    zl = zeta - m
    w = 1j / (2 * cmath.sqrt(-zl))  # 1/(4w²) = ζ', upper repelling petal
    for _ in range(30):
        s, d = phi_series(w, False)
        step = (s - zl) / d; w -= step
        if abs(step) < 1e-17 * abs(w): break
    s, d = phi_series(w, False)
    dw = 1 / d
    for _ in range(2 * m):
        dw *= dF(w); w = F(w)
    return w, dw

def psi(zeta, petal):
    """Ψ_petal(ζ) and its derivative (petal = +1 upper, -1 lower)"""
    if petal > 0: return psi_plus(zeta)
    w, dw = psi_plus(zeta - 0.5)
    return F(w), dF(w) * dw

def phi_a_entry(w):
    """(Φ_a, Φ_a', entering petal sign): the petal is the sign of Re of the point where the orbit enters"""
    dw = 1.0 + 0j
    for n in range(100000):
        if abs(w) < R_LOC and abs(w.imag) < 0.7 * abs(w.real):
            s, d = phi_series(w, True)
            return s - n / 2, d * dw, (1 if w.real > 0 else -1) * (1 if n % 2 == 0 else -1)
        dw *= dF(w); w = F(w)
        if abs(w) > 10: return None
    return None

def transit(w, sigma, cross=False):
    """g_σ(w), ∂/∂w and ∂/∂σ; the starting point's petal (sign of Re after an even number of steps) decides the exit"""
    r = phi_a_entry(w)
    if r is None: return None
    p, dp, petal = r
    q, dq = psi(p + sigma, -petal if cross else petal)
    return q, dq * dp, dq
