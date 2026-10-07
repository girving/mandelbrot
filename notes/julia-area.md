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

## The parabolic cusp c → 1/4 (2026-10-06/07)

c = 1/4 − ε, ε = 2^-k, k = 3 … 23, in double (`scratch/julia/cusp_scan.py`, fits `scratch/julia/cusp_fit.py`).

**Area.**  area K(1/4 − ε) = A₀ + 3.67368043 ε − 10.567656 ε^{3/2} + … in half-integer powers from ε¹; integer powers
alone fit badly and a √ε term comes out zero.  Fits over k ≥ 9 give the cauliflower

  **area K(1/4) = 2.93000204710204**, stable to ~1e-13 across power sets,

against Monte Carlo 2.93000 ± 0.00008 (4.3e9 samples; the disk |z| ≤ 1/2 maps into itself, so orbits entering it
are certified interior).

**Spectral gap.**  L's leading eigenvalue is |f'(q)|^{-2} + O(ε²) (ρ − μ⁻² ≈ 21 ε²) at the repelling fixed point q,
μ = 1 + 2√ε: the variational bound from the delta measure at q is nearly sharp, and 1 − ρ ≈ 4√ε.  The other
eigenvalues near 1 form the cluster μ^{−2−a−b} of dilations in q's Koenigs coordinate.

**Escape tail at c = 1/4** (Monte Carlo, `scratch/julia/cusp_tail.cc`): area escaping after n steps T(n) ~ n^{-3.2}.
The exterior near the parabolic point is a thin cusp (width ~d² at distance d) and escape from distance d takes
~1/d steps.  M's tail is T(k) ≈ 1.27/k: a single dynamical parabolic point is 2–3 powers lighter, so M's heavy
tail must come from its dense set of roots and the different time scaling there (π/√δ in parameter space).

**Solver work this needed.**  h ~ 1/|z − 1/2| down to the scale √ε, so: Möbius angular grading (`--grade`);
the parabolic preconditioner (1 − L₊)⁻¹ = Π (1 + L₊^{2^j}) over the branch fixing q, applied through a 2×
oversampled grid plus local stencils (~⅓ matvec; local stencils on the collocation grid alone fail because L₊
carries even Nyquist-scale modes with weight ≈ 1 near q); a vectorized matvec; and `--fast`, an NUFFT-style
L (exponential of semicircle kernel) for the double corrections, O(N (nr + nt)), 2e-12 accurate.  A k = 20
solve went from ~30 min to ~1.5 min (before `--fast`).  A two-grid preconditioner did not pay.

**Direct solve at c = 1/4** (prototype `scratch/julia/parabolic_proto.py`, dense).  With r1 = 1/2 the formula still
holds, but h ≈ c₀₀/(3|z − 1/2|) at the parabolic point (Fatou coordinate: backward branch is ζ ↦ ζ − 1,
h ≈ Σ |ζ|⁴/|ζ − n|⁴).  Plain collocation is garbage near 1/2.  Enrichment h = h_reg + Σ c_ij χ H_ij with
H_ij = (1 − L₊)⁻¹ m_ij (orbit sums), a cutoff χ (H alone has a seam on the principal branch's cut), and c_ij
pinned by Taylor constraints on s = 1 + L₋h (c₀₀ = 1 + h(−1/2) ≈ 3.4637) captures the singular profile and gives
2.9300020(5) ± 2e-7 at 60 × 121, an independent check of the extrapolation without assuming its power series.
Convergence is only algebraic: h_reg keeps weak singularities at every scale (poles of the complexified Fatou
sum at z = 1/2 − 1/n), so full precision would need a Fatou-coordinate patch.
