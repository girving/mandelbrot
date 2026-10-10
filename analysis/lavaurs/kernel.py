"""The universal island kernel: children of a small Lavaurs-model component.

For a small single-transit component U (center σ_U), the critical orbit's first excursion passes the critical point at
u = D_U (σ - σ_U) (D_U = ∂(F^{n_U} ∘ Ψ)/∂σ, |D_U| ~ C_U^{-1/4}), so F of it is v + u², and the next exit coordinate is
p_2 = Φ_a(v + u²) + σ.  A child over the single-transit target t at shift j solves Q(u) = Δ, Δ = σ_t + j/2 - σ_U, with
Q(u) = Φ_a(v + u²) - ζ0.  From the quadratic normal form of the return map (area_σ = A_card / |A D|² for a small
component, A = R''/2, D = ∂R/∂σ):
    C_child = C_U C_t |Φ_a'(v)|² / (K |Q'(u)|⁴),   K = (π²/4) A_card,
with relative corrections O(C_U^{1/4}).  Both preimages ±u count (the two branches).

  kernel(Δ) -> list of (u, weight = C_child / (C_U C_t))"""
import cmath, math
from lavaurs import phi_a, phi_a_entry

V = -0.25 + 0j
_Z0 = None
A_CARD = 3 * math.pi / 8
K = math.pi ** 2 / 4 * A_CARD

def zeta0():
    global _Z0
    if _Z0 is None:
        e = phi_a_entry(V); _Z0 = (e[0], phi_a(V)[1])
    return _Z0

def Q(u):
    r = phi_a(V + u * u)
    if r is None: return None
    s, d = r
    return s - zeta0()[0], 2 * u * d

def solve(delta, u0, its=60):
    u = u0
    for _ in range(its):
        r = Q(u)
        if r is None or r[1] == 0: return None
        f, df = r
        step = (f - delta) / df
        u -= step
        if abs(step) < 1e-13 * (1 + abs(u)): return u
        if abs(u) > 3: return None
    return None

def kernel(delta, grid=12, R=1.5):
    """All solutions u of Q(u) = Δ found from a grid of starts in |u| < R (u and -u both), with weights"""
    cv = zeta0()[1]
    sols = []
    for a in range(grid):
        for b in range(grid):
            u0 = complex(-R + 2 * R * (a + 0.5) / grid, -R + 2 * R * (b + 0.5) / grid)
            u = solve(delta, u0)
            if u is None: continue
            if not any(abs(u - s) < 1e-8 or abs(u + s) < 1e-8 for s in sols): sols.append(u)
    out = []
    for u in sols:
        dq = Q(u)[1]
        w = abs(cv) ** 2 / (K * abs(dq) ** 4)
        out += [(u, w), (-u, w)]
    return out

if __name__ == '__main__':
    # D2: U = S1, t = bulb, j = 0 (S1 is not small: expect O(1) corrections)
    s1 = complex(0.17957995250940523, 0.63460217927426621); bulb = complex(-1.0074583370365449, 0.16135210336429348)
    for j in (0, -1, -2):
        ks = sorted(kernel(bulb + j / 2 - s1), key=lambda x: -x[1])
        print('S1 over bulb, j=%d: Δ=%s; predicted C/(C_U C_t) for the solutions: %s' % (j, bulb + j / 2 - s1, ['%.3f' % w for u, w in ks[:6]]))
    print('actual D2: C/(C_U C_t) = %.3f; the R-branch partner: %.4f' % (1.8577745536e-4 / (1.724974063498964762e-3 * 0.20812104826488634), 4.9974309665762845e-06 / (1.724974063498964762e-3 * 0.20812104826488634)))
