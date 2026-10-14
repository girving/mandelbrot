"""The straightened parameter of a component in a source U's copy: c(σ) = A (R_U(1/2) - 1/2), R_U(w) = F^{n_U}(g_σ(F(w)))
the source's own return map (quadratic-like near the critical point w = 1/2) and A = R_U''(1/2)/2 =
(F^{n_U} ∘ Ψ)'(ζ0 + σ) Φ_a'(v).  At the center of U*X, c(σ) ≈ c_X (exactly in the small-copy limit): a test for tuning
that needs no guess of where U*X lies."""
from lavaurs import F, dF, psi, phi_a_entry, phi_a

V = -0.25 + 0j

def straighten(sigma, n_u):
    e = phi_a_entry(V)
    if e is None: return None
    p, dp, petal = e
    x, dx = psi(p + sigma, petal)
    d = dx * dp                       # d/dv of Ψ(Φ_a(v) + σ)
    for _ in range(n_u): d *= dF(x); x = F(x)
    return d * (x - 0.5)

if __name__ == '__main__':
    tests = [('bulb*H', complex(-1.0407424741922733, 0.36350852754039714)), ('bulb*A', complex(-1.0833588520516131, 0.48701507889481238)),
             ('bulb*Q', complex(-1.1127436313454888, 0.54137059618035221)), ('3.5e-5 r=4', complex(-0.80458, 0.32210)),
             ('bulb 1/4 sat', complex(-0.82360, 0.17627)), ('bulb 3/4 sat', complex(-1.17167, 0.14431))]
    for nm, s in tests: print('%-12s c = %s' % (nm, straighten(s, 1)))

def passages(sigma, n_u, r):
    """For each transit i < r: |F^{n_u}(x_i) - 1/2| (distance of U's critical passage) and the steps from x_i to the
    next gate entry.  A U-tuned component (U single-transit, excursion n_u) passes near the critical point after
    exactly n_u steps at every transit."""
    from lavaurs import R_LOC
    e = phi_a_entry(V); p, _, petal = e; p = p + sigma
    out = []
    for i in range(1, r):
        x, _ = psi(p, petal)
        w = x
        for _ in range(n_u): w = F(w)
        dist = abs(w - 0.5)
        w2, steps = x, 0
        while not (abs(w2) < R_LOC and abs(w2.imag) < 0.7 * abs(w2.real)) and steps < 100000: w2 = F(w2); steps += 1
        out.append((dist, steps))
        e = phi_a_entry(x); p = e[0] + sigma; petal = e[2]
    return out

def tuned_by(sigma, r, n_u):
    """|R_U^r(1/2) - 1/2| with R_U(w) = F^{n_u}(g_σ(F(w))) the return map of the single-transit source U (excursion
    n_u) at the phase σ.  NOT a tuning test: F commutes with g_σ, so R_U^p = F^{p n_u + p - 1} g^p F is the return map
    of any p-transit component with that excursion count; this only checks n ≡ -1 (mod p) in the given
    representative (necessary for tuning, not sufficient)."""
    from lavaurs import transit
    w = 0.5 + 0j
    for _ in range(r):
        t = transit(F(w), sigma)
        if t is None: return float('inf')
        w = t[0]
        for _ in range(n_u): w = F(w)
        if abs(w) > 10: return float('inf')
    return abs(w - 0.5)
