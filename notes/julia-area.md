# Areas of filled Julia sets by the area transfer operator

Goal: compute μ(K(c)) for hyperbolic c to many digits, as a first step toward faster estimates of μ(M)
(see the discussion of 2026-10-06: dynamical self-similarity makes K(c) easy where M's parabolic points make
it hard).

## Method (`julia.h`, `julia.cc`, `julia_area`)

L g(z) = Σ_{f(w)=z} g(w)/|f'(w)|² pulls area back under f(z) = z² + c.  On an annulus A = {r1 < |z| < r2} with
r1 − |c| > r1² and r2² − r2 > |c| (so |c| < 1/4), every point of A either escapes through
X = {z ∈ A : |f(z)| ≥ r2} or lies in K, and

  area K(c) = π r2² − ∫_X h,   h = (1 − L)⁻¹ 1.

f⁻¹(A) ⊂⊂ A makes L compact on analytic functions and h analytic, so collocation converges exponentially.
Chebyshev points in log r (the nearest singularity of h is at |z| = |c|) times odd equispaced angles;
barycentric interpolation at the preimages ±√(z − c), which share a radius, so L is stored as N × nr and
N × nt factors and applied in O(N²) without forming it.  Geometry comes from arb.  The solve is iterative
refinement: residuals in S ∈ {double, Expansion<2>, Expansion<3>}, corrections by GMRES in double, ~15 digits
per round.  The area integral is trapezoid in θ times Gauss–Legendre in r over [ρ(θ), r2].

## Results (Expansion<3>, 18-core M5 Pro)

| c | area K(c) | agreement | time (100 × 301) |
|---|---|---|---|
| 0 | π | 1.4e-48 (exact) | — |
| −0.2 | 3.06738856513448320124239642202731435742915152392 | 1.5e-47 over three annuli | 4.6 s |
| 0.15 + 0.15i | 3.02818653139968324843367793715458647564610444333 | 1e-47 (100×301 vs 120×361) | 5 s |
| −0.24 | 3.03633886599404550603058017155174706553220476153 | 1e-47 | 3.8 s |
| 0.24i | 3.01781547412201129153904150682454212709438613264 | 1e-47 | 4.9 s |
| 0.24 | 2.957848470000562745298079018889773… | 3.5e-36 (140 × 421) | 25 s |

(c is the nearest double to the decimal shown.)  Near the parabolic cusp c = 1/4 convergence slows, as
expected: the annulus barely contracts and L's spectral radius approaches 1 (GMRES iterations 150 vs 85).
Time splits evenly between the Expansion<3> residual matvecs and the double GMRES; a fused Expansion dot
product would speed the former.

Reproduce: `julia_area --c -0.2 0 --nr 100 --nt 301 --prec td [--r1 0.33 --r2 1.3] [--verbose]`.
