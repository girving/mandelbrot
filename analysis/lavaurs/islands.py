"""Two-transit components as island preimages: Θ(σ) = H(ζ0 + σ) + σ - ζ0 with H = Φ_a ∘ Ψ the horn map; at a
single-transit center σ_u (excursion n_u) H has a critical point and Θ(σ_u) = σ_u - (n_u + 1)/2 (check); the two-
transit components with final excursion n_c are the σ with Θ(σ) = σ_c (+ j/2 shifts) for single-transit centers σ_c.
Prototype (double)."""
import cmath, math, sys
import model_area as M
from lavaurs import phi_a_entry

zeta0, _, petal0 = phi_a_entry(-0.25 + 0j)

def Theta(sigma):
    """Θ(σ), Θ'(σ)"""
    q0, q1, q2 = M.psi3(zeta0 + sigma, petal0)
    p0, p1, p2, petal = M.phi_a3(q0)
    return p0 + sigma - zeta0, p1 * q1 + 1, (p2 * q1 * q1 + p1 * q2)

def solve(target, s, iters=60):
    for _ in range(iters):
        t, d, dd = Theta(s)
        st = (t - target) / d; s -= st
        if abs(st) < 1e-15: return s
    return None

if __name__ == '__main__':
    S2, S1 = -0.28523807327370199 + 0.95845771894111875j, 0.17957995250940523 + 0.63460217927426621j
    t, d, dd = Theta(S2)
    print('Θ(S2) = %s,  S2 - (0+1)/2 = %s,  Θ\'(S2) = %s, Θ\'\' = %s' % (t, S2 - 0.5, d, dd))
    # preimages of σ_c = S1 (+ j/2) in S2's island: quadratic seeds σ ≈ S2 ± sqrt(2(target - Θ(S2))/Θ'')
    for j in range(-3, 4):
        target = S1 + j / 2
        for sgn in (1, -1):
            seed = S2 + sgn * cmath.sqrt(2 * (target - t) / dd)
            s = solve(target, seed)
            if s is None: continue
            for n in range(0, 4):
                M.TRANSITS = 2
                try:
                    A, c, conv = M.area(s, n)
                except Exception:
                    continue
                if abs(c - s) < 1e-8:
                    print('j=%+d branch %+d: σ = %.10f%+.10fi  n=%d  (π²/4) area = %.10e  ratio to C(S1) %.5f  |Θ\'|^-2 = %.5f' %
                          (j, sgn, s.real, s.imag, n, math.pi**2 / 4 * A, math.pi**2 / 4 * A / 1.7249740635e-3, abs(Theta(s)[1])**-2))

def track(center, target, sgn, steps=400):
    """Follow Θ(σ) = Θ(center) + t (target - Θ(center)), t: 0 → 1, on branch sgn (Θ ≈ Θ_c + δ + Θ''δ²/2 at the center)"""
    t0, d, dd = Theta(center)
    s = None
    for i in range(1, steps + 1):
        tgt = t0 + (i / steps) * (target - t0)
        if s is None:  # quadratic seed with the linear term
            b, c = d / dd * 2, -2 * (tgt - t0) / dd
            s = center + (-b + sgn * cmath.sqrt(b * b - 4 * c)) / 2
        for _ in range(40):
            v, dv, _ = Theta(s)
            st = (v - tgt) / dv; s -= st
            if abs(v - tgt) < 1e-11 * max(1, abs(tgt)) or abs(st) < 1e-15: break
        else:
            return None
    return s
