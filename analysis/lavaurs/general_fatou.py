"""Formal Fatou coordinate at the root of the p/q bulb (prototype, double precision via numpy).

At c_0 = λ/2 - λ²/4 (λ = e^{2πi p/q}) the fixed point α_0 = λ/2 has multiplier λ; in w = z - α_0 the map is
f(w) = λ w + w², and f^q(w) = w + A w^{q+1} + … has q attracting and q repelling petals permuted by f.  The Fatou
coordinate satisfies Φ(f(w)) = Φ(w) + 1/q.  Ansatz Φ(w) = Σ_{j=-q}^{N} a_j w^j + β log w: the coefficient of w^m in
Φ(f(w)) - Φ(w) involves a_m (λ^m - 1) plus lower-order a's, so a_m is determined for m ≢ 0 (mod q), while for m ≡ 0
the equation at the next order fixes the remaining freedom (a_m for m ≡ 0 is a gauge; we set a_0 = 0 and solve the
resonant orders for the coefficients they determine).  Here: solve the truncated linear system by least squares and
check the functional equation numerically.

  python3 general_fatou.py p q"""
import sys, cmath, math
import numpy as np

def series_matrix(p, q, N, Nrows=None):
    """Unknowns: a_{-q..N} (excluding a_0) and β.  Equations: coefficients of w^m, m = -q+1 .. N, of
    Φ(f(w)) - Φ(w) - 1/q, where f(w) = λ w + w² = λ w (1 + w/λ)."""
    lam = cmath.exp(2j * math.pi * p / q)
    M = N + 2 * q + 4          # series length for expansions
    idx = [j for j in range(-q, N + 1) if j != 0]
    ncol = len(idx) + 1        # + β
    rows = list(range(-q + 1, (Nrows if Nrows is not None else N) + 1))
    A = np.zeros((len(rows), ncol), dtype=complex)
    rhs = np.zeros(len(rows), dtype=complex)
    # (λ w (1 + w/λ))^j = λ^j w^j Σ_k binom(j, k) (w/λ)^k  (generalized binomial)
    for col, j in enumerate(idx):
        coef = 1.0 + 0j
        for k in range(0, M):
            m = j + k
            if k > 0: coef *= (j - k + 1) / k / lam
            term = lam ** j * coef
            if m in rows:
                A[rows.index(m), col] += term
        if j in rows: A[rows.index(j), col] -= 1
    # β L(w) with L = log(w^q)/q on each petal (a branch per petal, invariant under w → λ w): its change under f is
    # log((1 + w/λ)^q)/q = log(1 + w/λ)
    col = len(idx)
    for k in range(1, M):
        m = k
        if m in rows: A[rows.index(m), col] += (-1) ** (k + 1) / k / lam ** k
    if 0 in rows: rhs[rows.index(0)] = 1.0 / q
    return A, rhs, idx, lam

def solve(p, q, N=30):
    """Only the constant a_0 is a gauge (set to 0).  The resonant a_m (m ≡ 0 mod q) are fixed by the order above them,
    so unknowns run to N + 2q and equations to N + q; coefficients up to N are trusted."""
    A, rhs, idx, lam = series_matrix(p, q, N + 2 * q, N + q)
    keep = list(range(len(idx) + 1))
    sol, res, rank, sv = np.linalg.lstsq(A[:, keep], rhs, rcond=None)
    a = {idx[i]: sol[n] for n, i in enumerate(keep[:-1]) if idx[i] <= N}
    beta = sol[-1]
    return a, beta, lam, np.abs(A[:, keep] @ sol - rhs).max()

def phi(w, a, beta, branch=0):
    """Φ on a petal: L(w) = log(w^q)/q on the branch log(w^q) + 2πi branch"""
    q = -min(a)
    return sum(c * w ** j for j, c in a.items()) + beta * (cmath.log(w ** q) + 2j * math.pi * branch) / q

if __name__ == '__main__':
    p, q = int(sys.argv[1]), int(sys.argv[2])
    a, beta, lam, res = solve(p, q)
    print('p/q = %d/%d: λ = %s; linear-system residual %.1e' % (p, q, lam, res))
    print('  leading a_{-q..-1}:', ' '.join('%.6g%+.6gi' % (a[j].real, a[j].imag) for j in range(-q, 0)), ' β =', beta)
    # check the functional equation at sample points in an attracting petal direction (small |w|)
    for r in (0.02, 0.01):
        errs = []
        for k in range(8):
            w = r * cmath.exp(2j * math.pi * (k + 0.37) / 8)
            fw = lam * w + w * w
            # compare Φ(f(w)) - Φ(w) with 1/q, the log handled by the continuous change log(fw) - log(w) = log λ + log(1 + w/λ)
            d = sum(c * (fw ** j - w ** j) for j, c in a.items()) + beta * cmath.log(1 + w / lam)
            errs.append(abs(d - 1 / q))
        print('  |Φ(f(w)) - Φ(w) - 1/q| at |w| = %.2f: max %.2e' % (r, max(errs)))
