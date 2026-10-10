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

## Satellite roots: centering, and dynamical escape tails (2026-10-07)

**Centering at α.**  `julia_area` now centers the annulus at the attracting fixed point α when |c| ≥ 1/4 (a disk
|z − α| ≤ r1 maps into itself for r1 < 1 − |λ|).  But the interior disk must contain the critical value: L's
weight 1/(4|z − c|) makes h singular like 1/|z − p| at every postcritical point p in A.  Around α that needs
|α||1 − α| < 1 − |λ|, so on the real axis only c ∈ (−0.394, 0.236).  Near the satellite roots the critical orbit
creeps into α through a pinch (the 2-cycle at α ± i√ε near −3/4), and any forward-invariant region containing it
is far from round: reaching them needs a boundary-fitted interior region, e.g. a level curve of α's Koenigs
coordinate, or patches.

**Escape tails at parabolic points** (Monte Carlo, `scratch/julia/satellite_tail.cc`, 2.7e8 samples; interior
certified in the attracting petals Re u > 50 of the Fatou coordinate u = −1/(q a w^q), f^q = w + a w^{q+1} + …):
the area escaping after n steps behaves like

  T(n) ~ n^{−(1 + 2/q)}  for a parabolic point with q petals:

| c | q | predicted D_j/D_{j−1} | measured (octaves 7–10) |
|---|---|---|---|
| 1/4 (cusp) | 1 | 0.125 | n^{−3.2} overall |
| −3/4 | 2 | 0.250 | 0.205, 0.228, 0.250, 0.225 |
| root of the 1/3 bulb | 3 | 0.315 | 0.250, 0.275, 0.298, 0.285 |
| root of the 1/4 bulb | 4 | 0.354 | 0.264, 0.308, 0.340, 0.324 |

Reason: exterior points near the parabolic point lie in a band of bounded width in the Fatou coordinate u ∝ w^{−q},
where f^q is translation, escape time is ~Re u, and dA_w ∝ |u|^{−2−2/q} dA_u; integrating the band beyond time n
gives n^{−1−2/q}.

**For parameter space:** the exponent tends to 1 as q → ∞, matching M's tail T(k) ≈ 1.27/k.  A hypothesis worth
testing: M's 1/k law is a superposition over satellite-type parabolic points of all periods, each with a tail
k^{−1−2/q} (with parameter-space time scaling, which differs from the dynamical one), so that the sum over q
(roots weighted by their sizes) gives 1/k up to slowly varying corrections.

**Band constants** (`analysis/parabolic_band.cc`, local; `analysis/satellite_tail.cc`, global; 2026-10-07).  In the
repelling Fatou coordinate u ≈ −1/(q a w^q) exterior points form a translation-invariant band, so near α
T_local(n) ≈ B_q n^{−(1+2/q)}, B_q = |E₀| |a|^{−2/q} q/(q + 2), with |E₀| the exterior area per unit length of the
repelling Écalle cylinder.  Every preimage p of α hosts a scaled copy weighted by |(f^m)'(p)|^{−2}, so the whole
tail is W B_q n^{−(1+2/q)} with W = 1 + |f'(−α)|^{−2} h(−α) = 1 + h(−α) (|f'(−α)| = |λ| = 1).

| root | q | \|a\| | local B_q | \|E₀\| | global B_total | W |
|---|---|---|---|---|---|---|
| cusp 1/4 | 1 | 1 | 1.56 ± 0.04 | 4.7 | ≈ 5.6 (noisy) | ≈ 3.5 (c₀₀ = 1 + h(−½) = 3.464 directly) |
| −3/4 | 2 | 2 | 0.640 ± 0.01 | 2.56 | 2.13 | 3.3 |
| 1/3 bulb root | 3 | 4.58 | 0.405 ± 0.01 | 1.86 | 1.37 | 3.4 |
| 1/4 bulb root | 4 | 11.7 | 0.294 ± 0.006 | 1.51 | 1.02 | 3.5 |
| 1/5 bulb root | 5 | 32.7 | 0.220 ± 0.01 | 1.24 | | |
| 2/5 bulb root | 5 | 28.1 | 0.224 ± 0.01 | 1.19 | | |

Observations: |E₀| depends on q but barely on p (1/5 vs 2/5 agree to 4%), falling roughly like 0.46 + 4.2/q;
B_q ≈ 1.1–1.6/q; and the total preimage weight W ≈ 3.3–3.5 at every root measured, i.e. h(−α) ≈ 2.4 (h is
2.3–2.5 at generic points of A at the cusp too).  The cusp's W agrees with the enrichment prototype's c₀₀.

## Parameter space: the tail's shape is universal (2026-10-07)

Boxes of half-width R (the bulb radius ~ sin(πp/q)/q²) centered on ten roots (cardioid 1/2 … 1/7, 2/5, 3/7; −5/4;
the primitive −7/4), `escape_tree --box … --base 32 --depth 8`, escape-time octaves 2^6 … 2^24, on a cluster CPU
node (`analysis/results/roots-boxes.log`, ~5 min total).

* Every box decays like the whole of M: escaping area per octave ∝ 2^{−1.03 … −1.04 j}, i.e. T(k) ~ 1/k with the
  same slow drift.  A box around a root does not isolate that root: it contains component-boundary arcs lined with
  satellite roots of every p/q, so 1/k is a property of the boundary everywhere (single roots are much steeper:
  ~k^{-3} satellite, ~k^{-5} primitive, matching the Böttcher octave energies j^{-4}, j^{-6} near roots).
* The ratio of the global tail (production run) to each box's tail is constant to within statistics over
  j = 14 … 23, while the global k·D(k) itself drifts 20% (1.698 → 1.366):

| box | 1/2 | 1/3 | 1/4 | 1/5 | 2/5 | 1/6 | 1/7 | 3/7 | −5/4 | −7/4 |
|---|---|---|---|---|---|---|---|---|---|---|
| global / box | 3.48 | 12.7 | 32.3 | 98 | 31.8 | 202 | 386 | 59.7 | 27.0 | 575 |

  So the tail separates, T_region(k) = w(region) Φ(k), with one universal time profile Φ for all of ∂M (satellite
  and primitive regions alike) and only the amplitude depending on location.  The amplitudes do not collapse under
  simple size scaling (w ∝ R^{1.4} across q at p = 1, but ∝ R^{2.3} across p at fixed q).

**Use for μ(M).**  The deep tail's shape can be measured in a small box instead of over all of M, and the
amplitude at moderate k.  E.g. stop the global run at k0 = 2^26 and get A(2^26) − A(2^32) (≈ 3.8e-8) as
w Σ Φ, which needs w to ~3e-4 relative for 1e-11, cheap at moderate k; Φ beyond 2^32 from a deep run in one box
replaces the extrapolation.  This reduces the tail and deep-orbit cost, not the statistical error of A(k0) itself.

**More roots, and near-parabolic crossovers** (`analysis/results/parabolic-univ.log`, cluster CPU, 2026-10-07):

| root | cycle | global B_total | local B (summed over cycle) | W |
|---|---|---|---|---|
| 1/6 | P 1, q 6 | 0.68 ± 0.02 | (run invalid: r0 inside r_min, \|a\| = 101) | |
| 2/7 | P 1, q 7 | 0.62 ± 0.02 | (invalid, \|a\| = 207) | |
| 3/7 | P 1, q 7 | 0.61 ± 0.02 | (invalid, \|a\| = 205) | |
| −5/4 | P 2, q 2 | 2.2 ± 0.2 | 0.465 ± 0.01 | 4.7 ± 0.5 |
| −7/4 | P 3, q 1 (primitive) | too sparse | 0.65 ± 0.03 | |

Exponents n^{−(1+2/q)} hold at q = 6, 7 and for the 2-cycle at −5/4 and the primitive 3-cycle at −7/4.  The
global constants at 2/7 and 3/7 agree to 2%, like 1/5 vs 2/5 locally: the constants depend on q, not p.  But W is
not universal across root types: 4.7 at −5/4 (a parabolic 2-cycle off the cardioid) against 3.3–3.5 at the
cardioid's roots.

Near-parabolic (c inside the cardioid near −3/4, 1 − |λ| = δ): the n^{-2} tail holds up to a crossover n* and then
collapses; n* grows like roughly δ^{-0.7} between δ = 2^-9 and 2^-12 (between the 1/δ and 1/√δ guesses), and the
curves do not yet collapse cleanly.  Near the cusp the tail is too thin for uniform sampling; both need
importance sampling near α.

**GPU runs** (`julia_tail.h/cc`, H200, `analysis/results/parabolic-gpu.log`; CPU and GPU histograms agree
exactly).  2^34 global samples at −3/4 take 27 s.

* Higher q on the cardioid: local B₆ ≈ 0.19 (1/6), B₇ ≈ 0.16 (1/7, 2/7, 3/7 agree); global 0.68 (1/6),
  0.60–0.62 (1/7, 2/7, 3/7).  W ≈ 3.6–3.8: the cardioid's roots keep W ≈ 3.4–3.8 for q = 1 … 7.
* −7/4 (primitive 3-cycle), 2^38 global samples: B_total still falling at octave 10 (4.30, 3.58, 2.86 ± 0.45), so
  W ≲ 4.4 there; not yet asymptotic.
* **Near-parabolic crossover collapses**: with δ = 1 − |λ| (λ the multiplier of α), near −3/4 the ratio of the
  tail to the root's, B_j(δ)/B_j(root), depends only on j + log₂ δ (≈ 0.71 at −2, 0.45 at −1, 0.15 at 0, for
  δ = 2^-8 … 2^-12), i.e. T(n; δ) ≈ T_root(n) F(n δ) with crossover n* ∝ 1/δ.  Near the cusp (δ = 1 − λ = 2√ε) the
  cutoff likewise moves one octave of n per octave of δ.  So in both cases the crossover time is 1/(1 − |λ|), the
  multiplier's distance from the unit circle.

For parameter space: the natural clock near a root is the multiplier, not the parameter.  Near a satellite root
1 − |λ| ∝ |c − c₀|; near a cusp 1 − λ ∝ √|c − c₀|; the parameter-space gate passage times have the same two
scalings, which suggests the dynamical crossover and the parameter-space escape time are one function of λ.

**The parameter-space clock near a root is the multiplier** (`analysis/root_escape.cc`, cluster CPU,
`analysis/results/gate-times.log`).  For exterior parameters c(λ) on circles |λ − λ0| = ρ around the cardioid's
p/q roots (2^22 angles, ρ = 2^-4 … 2^-14), the shortest escape times (10% quantile) satisfy

  n |λ − λ0| → 2π/q:  6.2831 (cusp), 3.1415 (1/2), 2.0943 (1/3), 1.5707 (1/4), 1.2566 (2/5), 0.8976 (1/7)

against 2π/q = 6.28319, 3.14159, 2.09440, 1.57080, 1.25664, 0.89760: the parabolic gate time, depending on q and
not p.  Higher quantiles sit at integer multiples of 2π/q (multiple passes through the gate).  The escaping
fraction of each circle is ∝ ρ (the thin exterior cusps between tangent components), with a coefficient
growing with q (0.3, 0.8, 1.5, 2.5, 4.7, 7.6 for q = 1, 2, 3, 4, 5, 7).  The dynamical crossover near a root
(n* ∝ 1/(1 − |λ|)) runs on the same clock.

**−7/4 at 2^40 samples** (107 s on one H200): global B_total plateaus at 3.45 ± 0.2 (octaves 9–10); with local
0.65, W ≈ 5.3 ± 0.4.  So W ≈ 3.4–3.8 on the cardioid's roots, 4.7 at −5/4 (2-cycle), 5.3 at −7/4 (primitive
3-cycle): constant within a family, not across families.

## The family expansion (2026-10-07)

Sum M's escape-time tail over types of root families rather than regions:

  T(k) = Σ_types ∫ G_type(k/τ) dN_type(τ),

with G_type the universal local profile of a family type (p/q satellite, primitive cusp, tuned copy), measured
once, and N_type the distribution of that type's occurrences over all of M by area and time scale, whose Mellin
transform ζ_type(s) = Σ_occurrences (scale)^s is the family's zeta function (computable from the component tree:
satellite scalings ~ sin(πp/q)/q², tuned copies by r_W as in D(s)).  One term per family then accounts for the
family's contribution wherever it occurs.  Premise under test: every occurrence of a root type looks alike in its
parent's multiplier coordinate, up to area scale |dc/dλ|² and a time shift log2 P (`analysis/multiplier_profile`).

**Premise test** (`analysis/results/family-profiles.log`, cluster CPU, 2^28 samples each): escape-octave profiles in
λ-area units near seven root occurrences, compared after shifting time by log2 P.  Ratio to the cardioid's
profile, flat to a few % over shifted octaves 6 … 16:

| occurrence | parent | q | amplitude vs cardioid |
|---|---|---|---|
| −5/4 | period-2 disk | 2 | 1.10 |
| 1/2 root of the 1/3 bulb | period-3 satellite | 2 | 1.17 |
| 1/2 root of the airplane | period-3 primitive (a cardioid) | 2 | 1.00 |
| 1/2 root of the 1/4 bulb | period-4 satellite | 2 | 1.25 |
| 1/3 root of the period-2 disk | period-2 disk | 3 | 1.10 |

So every occurrence of a type has the same time profile G_q in its parent's multiplier coordinate; amplitudes are
1 for cardioid-shaped (primitive) parents and 1.10–1.25 for disk-shaped (satellite) parents, perhaps from the
multiplier map's distortion over the sampled disk.  The family expansion's premise holds: T(k) = Σ_q G_q ⊛ N_q
with N_q the occurrences weighted by |dc/dλ|² (and a parent-shape factor near 1) and shifted in time by log2 P.

**Ownership: the families partition the exterior** (`analysis/family_owner`, `analysis/family_shares.py`,
`analysis/results/family-owner.log`).  λ-disks overlap, so the expansion needs each escaping parameter assigned to
one family.  Dynamics decides: near a p/q root the critical orbit passes the root's gate, lingering by a
near-parabolic p-cycle with m − 1 ≈ −q²(λ/λ0 − 1).  The owner is the outermost gate passed: the smallest Q ≤ 64
whose orbit returns |z_{i+Q} − z_i| < 0.1 and whose Newton-refined Q-cycle has |m − 1| < μ0, then the catalog root
of period Q nearest c.  Roots come from hyperbolic's component centers by continuing the own multiplier to 1.  Three
earlier rules failed instructively:
- "most frequent return period": a child's gate also returns at its parent's period, so the parent steals it;
- "smallest mean log return": orbits rotating near a bulb's cycle at high q′ (gate period 2q′ > scan) fall back
  to the bulb's period;
- "gate explains the escape time": on a circle |λ − λ0| = ρ, nρ comes in multiples of 2π/q, one per pass through
  the gate (root_escape: 50% one pass, 25% two, 10% ≳ 8).  Multi-pass escapes are seahorse-valley copies inside the
  root's disk, so the root's family has to include them.

Nested ownership does that.  Per root, the owned tail is ∝ 1/n (owned area × 2^j flat over octaves 8 … 16), so each
family term is A_r P_r g_q / n.  Within the cardioid's family, g_q collapses across occurrences (μ0 = 1/2: q = 5
0.064, 0.065; q = 7 0.034, 0.038, 0.038).  With g_q calibrated on the cardioid's roots alone,
Σ_r A_r P_r g_q predicts what each period's satellites own:

| period | μ0 = 1/2 | 1 | 2 | 4 |
|---|---|---|---|---|
| 5, 7 (cardioid only) | 1.01, 0.97 | 0.96, 1.02 | 0.99, 1.02 | 1.01, 1.01 |
| 4 (incl. −5/4) | 1.23 | 1.56 | 0.88 | 0.73 |
| 6 (incl. 1/2 of 1/3 bulb, 1/3 of 2-disk) | 1.23 | 1.75 | 0.42 | 0.25 |
| 8 | 1.18 | 1.73 | 0.67 | 0.45 |

Small disks show the parent-shape factor (1.1–1.5 for disk parents; the excess sits at the disk's rim, where
|m − 1| < μ0 bulges differently per parent, while inner λ-shells match the premise test's ≈1.0–1.1).  Large disks
nest bulb children into their ancestors' families, so those ratios fall below 1.

**Coverage** (share of the long-time tail, octaves 10 … 13):

| μ0 | gates ≤ 8 | cardioid's q ≤ 8 families | gates 9 … 64 | no gate ≤ 64 |
|---|---|---|---|---|
| 1/2 | 15% | 13% | 19% | 66% |
| 1 | 38% | 29% | 30% | 31% |
| 2 | 77% | 68% | 19% | 3% |
| 4 | 88% | 81% | 12% | 0% |

The shares are flat across octaves 8 … 19, so the families' sum is octave-independent as universality predicts.
With μ0 = 4, the cardioid's p/q families for q ≤ 8 own 81% of M's tail: q = 2 28%, 3 21%, 5 13%, 4 8%, 7 7.5%,
8 2.7%, 6 1.6%.  (Even q falls below the odd q around it: at μ0 = 4 each q gate disk reaches its neighbors, so
nesting is order-dependent.)  A_r = sin²(πp/q) for the cardioid, so the cardioid's family term is
Σ_q (Σ_p sin²(πp/q)) g_q(μ0).  Its sum over all q is the leading term of a sublinear accounting of M's tail: a few
calibrated profiles g_q times an explicit sum over roots.  Open: g_q(μ0) for large q (≈ q^-2.4 at μ0 = 4 over
q = 5 … 8), the residual 12–20% in gates beyond period 8 (children of bulbs and copies, each a family with its own
A_r), and turning an accurate tail model into μ(M) (μ(M) = area(n) − T(n), so a 1% model of T(n) only buys a
factor 100 over raw escape counting at the same n).

## Renormalization cascade over the component tree

Every flat truncation of μ(M) converges algebraically: escape time 1/n, Böttcher terms 1/log N, components by
period P^-2 at cost 2^P (each family term of the family expansion is A P g_q / n).  The exponential convergence for
K(c) came from an operator with a spectral gap; in parameter space the operators that contract are renormalization
operators.  Assume, as the CI already does, that μ(M) = Σ_W area(W) (density of hyperbolicity and μ(∂M) = 0), and
organize the components as a tree.  The question: is each child's area (parent's multiplier data) × (a universal
function of the combinatorics with fading memory)?  Tools: `analysis/bulb_areas` (exact satellite bulb areas of any
parent: Newton to the child's center, then the boundary by continuation, then Green's theorem relative to the
center), `analysis/satellite_tree`, `analysis/bulb_cascade.py`; results in `analysis/results/bulb-cascade.log`.

**The satellite tree is nearly everything.**  Of the hyperbolic area through period 16 (1.4987444), the components
reached from the cardioid by satellite bifurcations alone hold 99.92%, and primitive copies hold 0.079% (2% of the
period-16 area, growing slowly).

**Bulb areas have fading memory.**  Normalize the cardioid's p/q bulb as F = area q^4 / (π sin^2(πp/q))
(|c'(λ0)|^2 = sin^2(πp/q)).  F ∈ [1.04, 1.33] is not a smooth function of p/q (63/127: 1.073; 63/128: 1.241); it is a
function of the continued fraction, dominated by the last digits.  Fix the last L digits and vary the first one over
1, 2, 3, 5, 10.  The spread falls geometrically, about 0.53 per added digit 1 and 0.30 per added digit 2, i.e. like
q_suffix^-1.3 (−1.27, −1.31, −1.35 for three digit patterns).  The prefix dependence factorizes,
δF ≈ c(suffix) φ(a1), with φ's shape converging in L.  Large final digits have clean expansions,
F([2, n]) ≈ 1.060 + 0.84/n.  All this matches near-parabolic renormalization (Inou–Shishikura): the operator acts on
rotation numbers by the Gauss map and contracts, so outer levels are forgotten.  Equivalently, F(p/q) ≈ Ψ(y) for
y = (p^-1 mod q)/q = [0; a_k, …, a_1], the reversed continued fraction, with Ψ continuous.  At q ≥ 100 the spread of
F within bins of y is 0.011, against 0.17 in bins of p/q.  Prepending a digit moves y by O(q^-2), and Ψ is Hölder
with exponent about 0.65, which explains the q^-1.3 fading.

**Context is carried by the combinatorics.**  Bulbs of other parents, normalized by the parent's |c_W'(λ0)|^2 and
divided by the cardioid's F: away from the parent's root (p/q ∈ [1/4, 3/4], q ≥ 20) R = 1 within 0.06–0.25%
(median).  Near the parent's root (child 1/n) R dips to 0.90–0.96, with a profile set by the parent's own rotation
number.  Three parents attached at 1/2 (the period-2 disk, its 1/2 child at −1.31, the 1/3 bulb's 1/2 child) give
R(1/n) = 0.9466/0.9500/0.9505, 0.9348/0.9359/0.9356, 0.9473/0.9474/0.9472, 0.9623/0.9622/0.9622 at n = 5, 10, 20, 40.
Parents at 1/3, 1/4 and 2/5 each have their own profile.  The airplane copy's bulbs match the cardioid's to 2e-4
(primitive parents straighten almost conformally).

**But one level does not sum exponentially.**  The cardioid's bulbs total 0.291067, 0.292753, 0.293159, 0.293266
through q = 16, 32, 64, 128 (truncation error ~Q^-2).  The sum weighs F against sin^2(πp/q), and p^-1 is
decorrelated from p (Kloosterman cancellation), so only the per-q average F̄(q) matters.  F̄(q) drifts smoothly
(≈ 1.276 − 0.32 q^-0.54) but fluctuates arithmetically by about ±0.007 even over primes.  A drift model fitted on q ≤ 64
predicts the exact Σ_{64<q≤128} to 0.5%, cutting the truncation error 200× (1.1e-4 → 5.6e-7), but the rate stays
Q^-2.  An exponential method would need the sum over continued fraction words as a transfer-operator problem (Mayer's
operator for q^-4 weights is nuclear, so it converges super-exponentially for analytic weights).  But Ψ is only
Hölder, so that would also be algebraic unless the fading has its own spectral expansion (corrections q^-1.3, q^-2.6,
…, each with a factorized eigenfunction).  The factorized shape of δF hints at that, but it is untested.

Parameter-space reading: near-parabolic renormalization is a dynamical-space operator with a hyperbolic fixed point,
and parameter-space bulb areas inherit its contraction (fading memory in the Gauss map) and its context dependence
(the parent's word).  That is the dynamical → parameter transfer working.  Whether it yields an exponentially
convergent μ(M) depends on the spectral structure of the fading.

**The fading has a spectral expansion** (`bulb_cascade.py modes`).  Take periodic suffix families s_L = block^L +
[last] (blocks [1], [2], [3], [1,2], [4]), 12 prefixes u, and D(L, u) = F(u s_L) − F(u_0 s_L).  If the fading is a
sum of modes, D = Σ_k a_k λ_k^L φ_k(u), and a matrix pencil over the prefixes recovers the λ_k.  872 bulbs up to
q = 40000 (accurate to 1e-9: mirror pairs agree to 1.1e-9):
- D is numerically low-rank: singular values 0.94, 0.04, 0.01, 3e-4, 6e-5, 5e-6 for block [1]; 0.77, 0.04, 6e-4,
  5e-5, 1e-5, 3e-7 for [2].  Each added mode cuts the fit residual 10–30×; for [1,2] three modes leave 4e-9.
- The eigenvalues are set by the Gauss map's contraction along the periodic orbit [block]^∞: λ = (Π y_i^2)^γ with
  universal exponents.

| block | Π y_i^2 | γ (leading ± pair) | γ (next) |
|---|---|---|---|
| [1] | 0.382 | 0.678, 0.718 | 0.943 |
| [2] | 0.172 | 0.693, 0.698 | 1.045 |
| [3] | 0.0917 | 0.686, 0.697 | 1.044 |
| [1,2] | 0.0718 | 0.707, 0.707 | 1.033 |
| [4] | 0.0557 | 0.663, 0.775 | — (too few L) |

So F has the structure of a weighted Gauss transfer operator with weight |T'|^-γ: a leading pair λ = ±(Π y^2)^0.70
(the sign from the orientation reversal of each inverse branch), then γ ≈ 1.04, then more.  Equivalently,
Ψ(y) is Hölder with exponent ≈ 0.70 on continued fraction cylinders, with a discrete spectrum of corrections.  It is
a Brjuno-type function (the Brjuno function has γ = 1/2).  If F satisfies such a functional equation with analytic
data, the sum over all bulbs is the resolvent of a two-sided Mayer operator.  Mayer's operator is nuclear on
holomorphic functions, so that sum would converge super-exponentially: the exponential route.  Next test: fit the
functional equation directly and check that it reproduces the bulb sum to many digits.

**An exact telescoping identity puts the number theory in a Mayer operator** (`bulb_cascade.py telescope`).  With T
the Gauss map and Δ(x) = F(x) − F(Tx), F(x) = Σ_j Δ(T^j x).  Prepending a digit maps x ↦ 1/(a+x) and multiplies q by
a+x, so summing over prefixes gives, exactly,

  Σ_x π sin²(πx) F(x) / q^4 = Σ_x Δ(x) G(x) / q^4,   G = (I − L)^-1 π sin²(π·),   (Lf)(x) = Σ_a (a+x)^-4 f(1/(a+x)).

G is analytic, and Chebyshev collocation gets it to 1e-16 with 36 modes in 0.4 s (check: Σ_{a≥2} a^-4 G(1/a) =
(π/2)(ζ(3) − 1)/ζ(4)).  Because Δ fades, the right side converges faster.  Its increments per doubling of Q are 2.2e-4,
3.0e-5, 3.9e-6 (Q^-2.9), against 4e-4, 1.1e-4 (Q^-2) for the direct sum.  The cardioid's bulbs total
0.2933010 (Q = 128 telescoped, ± a few 1e-7).  The rate is limited by words with a large first digit, where Δ does not
fade.

**Large digits** (`bulb_cascade.py large`, to n = 16384).  This needed cusp coordinates in `bulb_areas`: in c, a bulb
of radius 1e-9 near 0.25 loses most of its digits.
- A large last digit is analytic: F([2, n]) − F_∞ halves per doubling of n (ratios 2.00–2.03), F_∞ = 1.06004.
- A large first digit (deep in the reversed continued fraction) is not.  F([n, 2]) → 1.2393643 and F(1/n) → 1.0808167
  with differences shrinking 3.3 → 2.2 per doubling: a fractional power (about n^-2γ ≈ n^-1.4 to n^-1.7) giving way to
  a small 1/n term.  So the parabolic-implosion limit is approached with the renormalization exponents.

**Functional equations tried.**
- Ψ(y) = α(y) + β(y) Ψ(Ty), with α, β polynomials: rms 2.9e-3 (from 9.3e-3 for α alone).  Not exact.
- Δ(x) p^{2γ} = h(x) (x = p/q, since the denominator of Tx is p): not smooth.
The structure the data point to is a linear response through the renormalization tower.  A front digit perturbs the
starting state, and the perturbation propagates through later levels by matrices M(b, x) depending on the forward
variable x (each level's renormalization depends on its own rotation number).  Then Σ_x Δ(x) G(x)/q^4 is the
resolvent of a matrix-valued Mayer operator (𝓛V)(x) = Σ_b (b+x)^-4 M(b, x) V(1/(b+x)).  That operator is still
nuclear when M is analytic, hence super-exponential.  The ± pair and γ ≈ 0.70, 1.04 are M's spectrum along periodic
orbits.  Next: learn M (K = 2–4 modes, analytic in x) from bulbs on many short words, and test whether the resolvent
reproduces the telescoped totals to many digits.

**Bulb areas are a rational series** (`analysis/wfa/`).  The response Δ(a·w) = F(1/(a + x_w)) − F(x_w) over 12,231
short words w (q ≤ 200) and front digits a = 1…10 is low-rank.  Its singular values are 109, 18, 8.2, 1.7, 0.066,
7.5e-3, 6.8e-4, 5.9e-5, 2e-6, about ×10 per mode after the fourth (111k bulbs, 14 min, all converged).  A linear
recursion c(b·w) = M(b, x_w)ᵀ c(w) on those coefficients fits only to 3e-3, uniformly in x_w and q_w, because the
single-digit prefixes do not span a space closed under adding a digit.

The right object is the Hankel matrix H[u, s] = F(u·s) over prefixes u and suffixes s.  If F is a weighted finite
automaton over the digit alphabet, F(w) = αᵀ N(a1)⋯N(a_{k−1}) β(a_k), then H has finite rank (Fliess).
- Singular values: 128 prefixes (q_u ≤ 14) × 1101 suffixes (q_s ≤ 60) give 470, 1.8, 0.49, 0.34, 0.18, 0.042, 0.028,
  9e-3, 5.7e-3, 1.9e-3, 8.6e-4, …, 5.8e-5.  256 prefixes (q_u ≤ 20) resolve 28 of them, down to 2.7e-6: about ×0.5
  per mode early, ×0.7 by mode 28.
- Learning: spectral learning (N(b) = P[U_b]⁺ P[U_b·b]) gives an 8-state automaton that predicts F on held-out
  long words up to q = 40000 with median error 2–3e-4.  With q_u ≤ 20 the best is 7.3e-5 at K = 10–12.  K is
  capped by how few prefix rows stay in the set after appending a digit.
- Large digits: first digits b ≤ 400 and last digits c ≤ 400 measured directly (7k bulbs; interior digits in the
  hundreds need BULB_TOL = 1e-9).  This cut the error on F([n, 2]) from 8e-3 to 3e-5.

**The automaton sums exactly.**  Build words by prepending digits.  Then

  V(x) = π sin²(πx) α + Σ_b (b+x)^-4 N(b)ᵀ V(1/(b+x)),   S_model = Σ_c c^-4 β(c)ᵀ V(1/c),

a matrix-valued Mayer operator (nuclear, solved by Chebyshev collocation).  As a control variate on the tail,
S ≈ S_exact(Q) + S_model − S_model(Q).  From bulbs with q ≤ 200 alone (K = 8–10), the estimates for Q = 64…200 are
0.29330126, 0.29330142, 0.29330146, 0.29330148, 0.29330149.  So the cardioid's bulbs total 0.2933015 ± ~1e-8.
That is about 1000× better than direct truncation at q = 200 (error 1.5e-5) and 50× better than plain telescoping.
(An earlier version that extrapolated the large digits looked flat at 0.29330159, but it was wrong on large first
digits by 8e-3; its flatness was luck.)

**Status.**  This is a dynamical → parameter result: near-parabolic renormalization makes bulb areas an
approximately rational series on continued fraction digits, and Mayer's operator sums it.  The accuracy is set by
the automaton's.  That accuracy improves with the Hankel block (2.1e-4 → 7.3e-5 when the prefix rows double), and the
singular values decay geometrically, so the method could converge exponentially.  But that is not shown yet, and
three things stand between this and many digits:
1. learning that isn't capped by prefix closure: shifted Hankel blocks H_b[u, s] = F(u·b·s) computed directly;
2. bulb data beyond double (Hankel singular values below ~1e-8 sit under the noise of 1e-11 bulb areas; hyperbolic's
   Expansion<2> polish would do);
3. the rest of the satellite tree: the parent-word dependence suggests an automaton over words of words, untested.

**The cardioid's bulbs to 14 digits** (`analysis/wfa/wfa_hp_sum.py`, `analysis/results/bulb-total.log`).  Two
pieces of engineering:
1. Double-double bulbs (`bulb_areas` with BULB_EXP=1).  Continuation in double, then Newton polish of the center and
   each boundary point in Expansion<2>, with Green's sum in Expansion<2>.  The 1/2 bulb comes out π/16 to 3e-31;
   N = 64 boundary points converge to ~1e-20 (N/2 subrule).  It is about 2× faster than the double mode, which
   computed N and 2N separately.
2. Learning without the prefix-closure cap.  N(b) for b ≤ 3 from suffix closure (H[u, b·s] = P N(b) Q[:, s]); for
   4 ≤ b ≤ 256 from shifted blocks F(u·b·s) on a 56 × 56 basis of short words chosen by pivoted QR (smallest
   singular values 4e-2, 3e-3 in orthonormal coordinates).  Cubic interpolation in 1/b between measured digits.
With clean data the Hankel singular values (256 × 1101) fall smoothly to 8e-10 by mode 68; held-out error on long
words (q > 2000) is 1.7e-8 at K = 40 and 4.5e-10 at K = 76 (median).

Exact head (every bulb with q ≤ 1000, double-double, fsum) plus the automaton tail, K = 36, 40, 44:

| | K = 36 | K = 40 | K = 44 |
|---|---|---|---|
| model total S_model | 0.2933015125432736 | 0.2933015125432689 | 0.2933015125432746 |
| estimate, Q = 200 | …432776 | …432779 | …432790 |
| estimate, Q = 1000 | …432848 | …432868 | …432870 |

The model reproduces the exact head through q = 1000 to ~1e-14 in total.  Chebyshev degree 40 → 56 changes
nothing; the digit cutoff 2000 → 4000 changes 2.5e-15.  The coarse 16 × 16 blocks used before for digits 65–400
had O(1) errors on words with such interior digits (−6e-10 in total); proper blocks bring that class to −1e-17.
**The cardioid's bulbs total 0.293301512543286 ± 1e-14**, from about 300k bulbs (an hour of laptop time).  Direct
summation of the same bulbs is good to 1.5e-6 (q ≤ 1000).

**The tree is a words-of-words automaton** (first checks).  Each satellite component is a sequence of rotation
numbers (r1, …, r_m).  Joining their continued fractions with a separator keeps the indexing a regular language, so
if each level's factors are automaton-representable, the tree sum is still one linear resolvent.
- Depth 2 is low-rank across the separator.  Bulbs of all 40 cardioid children with q1 ≤ 16, over child words
  q2 ≤ 30: the matrix F(r1; r2) has singular values 120, 0.58, 0.38, 0.10, 0.015, 0.011, 6.5e-3, 2.3e-3, 1.6e-3, ….
- The parent's derivative profile |c_W'(e^{2πiθ})|² is real-analytic around the circle (log's Fourier coefficients
  fall ~10× per harmonic; exactly constant for the 2-disk).  The geometric transport of the parent's curvature
  κ = λ c_W''/c_W' explains little of the child-to-child variation (residual 2.2e-3 → 1.7e-3 at first order), so the
  context is mainly combinatorial.
- Children-to-parent area ratios ρ(r1) (children with q2 ≤ 30) are 0.10–0.14, with the same continued fraction
  structure as F: ρ(1/n) → 0.1154, ρ([n, 2]) → 0.135, ρ([2, n]) → 0.124, deeper words 0.13–0.14.  Consistent with
  (μ(M) − 3π/8) / (cardioid bulbs) = 0.3285 / 0.2933 ≈ 1.12.
Plan: the subtree totals satisfy a linear fixed point τ(node) = 1 + Σ_children (A_child/A_node) τ(child).  With
the per-level factors as automata (digits, separator) times analytic profiles, that is a resolvent of one
matrix-Mayer operator for the whole satellite tree; the primitive copies (0.08%) enter through tuning.
- Context fades across levels about 20× per level.  Depth-3 bulbs (grandparents 1/2, 1/3, 2/5, 1/4, 3/7; parents
  1/2, 1/3, 2/5 of those; children q3 ≤ 20): at fixed parent, the grandparent changes F(r1, r2; r3) by median 5e-4
  (max 1e-2, at the smallest child words); at fixed grandparent, the parent changes it by 1.1e-2 (max 0.15).  Child
  profiles depend on the grandparent at 2e-4 of their size.
- Context is more than a new start vector.  Fitting F(r1; r2) = α(r1)ᵀ N(r2) β with the cardioid level's N, β and
  a free α per parent leaves 1e-3 in-sample and 3–5e-2 held out: near the parent's root the grandparent's parabolic
  point adds directions the root-level automaton never needed.  So the automaton has to be learned jointly over
  digits and the separator (`analysis/tree/`: gen.py, run.sh, learn.py).
- The full per-level cascade factor is low-rank across the separator.  λ(r1; r2) q2^4 = A(r1, r2) q2^4 / A(r1), which
  includes the parent's derivative profile, over 40 parents × 276 children has singular values (relative) 1, 0.1,
  0.04, 0.01, 5e-3, 1e-3, 7e-4, 3e-4, …, 6e-7 at mode 20: about ×0.5 per mode, like the within-level Hankel.

**Design for the tree sum (linear despite the multiplicative cascade).**  Per-level normalized factors (F and the
profile Ω) multiply along a path, so a model that makes each of them a linear functional of an automaton state would
need tensor-power states.  The way out is to make the cascade factor itself the linear object.  A node α carries a
K-dimensional context vector c(α), scale included, with each child's area linear in it:
A(α·r) q_r^4 = c(α)ᵀ f(r), and c(α·r) = c(α)ᵀ E(r).  The rank above says K ≈ 40 for 1e-12.  E(r) and f(r) are
functions of the child word r: analytic in x_r (profile) times rational in r's digits (F).  Read right to left as in
Mayer's recursion, the analytic factor is applied when a level's word is complete (x = x_r) and the digits through
matrices N(b).  So
  Σ_tree A = (3π/8) + Σ_{r1} q1^-4 c(r1)ᵀ Z,   Z = (I − Ē)^-1 f̄,   Ē = Σ_r q_r^-4 E(r),   f̄ = Σ_r q_r^-4 f(r),
with Ē and f̄ matrix-valued Mayer sums over all words.  Physically c(α) is a finite-dimensional proxy for the
derivative kernel c_α'(μ) conj(c_α'(ν)): children's areas are linear functionals of it, and children's kernels are
linear transforms of it (by composition with the child's multiplier map).

**Learning the tree's automaton (first attempts, `analysis/tree/learn2.py`–`learn4.py`).**
- A small mixed Hankel block over the extended alphabet (632 candidate prefixes, 361 suffixes, 80 × 80 basis) gives
  a poor automaton: median errors of 1e-2 even at the root level.  Its spectrum decays far more slowly (6e-5 at mode
  76) than the single-level block's, because a split inside a child word mixes context and within-level modes.
- Digits learned jointly from root rows and context rows (r1#u, 11 contexts r1 with q1 ≤ 8) over 489 digit suffixes
  predict the children of those contexts geometrically in K: depth-2 children are good to median 2e-6, 6.7e-8,
  6.4e-9, 4.9e-10, 3.8e-11 at K = 20, 30, 40, 50, 60.  So within a level the context enters as a start vector α(ctx),
  provided the digit matrices are learned with context rows present.
- The map from context to start vector is the hard part.  Regressing α(r1) on the root state after reading r1's
  digits (127 contexts, 70% training) leaves 1e-2 even in training, and analytic dependence on x_{r1} overfits.  The
  root-level state carries what predicts the cardioid's continuations, not what predicts r1's children; the joint
  state has to be learned from a block with many contexts, with the separator from row closure
  P[(r1, ·)] = P[root r1] N(#).  A context-rich block (128 contexts × short prefixes as rows, digit suffixes and
  separator futures as columns, ~226k bulbs) is computing.
- Context-rich block (`learn5.py`; 1457 rows: root prefixes and 111 training contexts r1 (q1 ≤ 20) × prefixes
  q_u ≤ 4; 617 columns: digit suffixes and separator futures; ~226k bulbs, ~55 min).  Its spectrum decays smoothly
  (2.3e-6 at mode 95).  The separator from row closure fits its training contexts well (residual 3.7e-4 at K = 20,
  4e-9 at K = 100).  But digit matrices from row closure are unstable: they are fixed only on the span of training
  row states, so long words amplify spurious eigenvalues and predictions blow up.  Digits want suffix closure (as in
  learn3, which is stable); the separator wants many contexts.
- Contexts span many directions (`ctx_pca.py`).  F(r1; s) − F_card(s) over 127 contexts × 489 suffixes has singular
  values 10, 0.8, 0.6, 0.1, 0.04, …, 4e-6 at 30: ×0.6–0.7 per direction.  Weighted by each suffix's share of the
  sum, 2e-8 at 30.  So ~40 coefficients per context, readable from ~40 children (not ~500), but the map from r1 to
  those coefficients still has to be learned.

**Where the tree stands.**  Proven so far: within a level, context enters as a start vector of a jointly learned
digit automaton, geometrically in K (3.8e-11 at K = 60); context fades ~20× per level; every quantity measured
(normalized areas, the full cascade factor, profiles) is low-rank across separators with geometric spectra.  Open:
(1) a stable learner for the separator, i.e. the map from a context's combinatorics to its start vector (suffix
closure for digits, separator from many contexts); (2) a linear form for the multiplicative cascade, since products
of per-level automaton outputs need tensor-power states unless the cumulative product itself has low rank, which
fast fading suggests (to be measured at depth ≤ 3); (3) the primitive copies (0.08% of the area, via tuning, with
copy bulbs universal to 2e-4).
- `learn6.py` (stable suffix-closure digits on a digit-only block with 111 contexts, separator by row closure):
  the separator residual is 0.13.  A context's start vector is not a linear function of the root-level state of r1.
  That state predicts the cardioid's bulbs near r1, not r1's own children.  A joint state needs separator-visible
  directions, and stable digit learning in those directions needs columns w#s' (digits, then a separator), which by
  brute force is ~2M bulbs.  Next: a better representation rather than more data.  Candidates: per-context
  coefficients on ~40 context directions, read from ~40 children each, with their own fading-memory model in r1; or
  a physically motivated context state (the parent's multiplier-map jet and the near-parabolic germ at its root).

**The cascade is effectively linear** (`analysis/tree/cascade_rank.py`, `context_rank.py`, `depth_saturation.py`;
about 130k small bulbs, minutes on the laptop).  Products of per-level automaton outputs would need tensor-power
states unless the cumulative area has low rank.  It does, nearly.
- The context map's rank (`context_rank.py`, rows: 91 prefixes u of a cardioid-level parent word; columns (s, r2)).
  For a single child r2, the child factor λ(u·s; r2) q2^4 has rank 34–39 at 1e-10, like the parent's own F (31).
  For all 21 children together it is 77: about 2.5× the scalar rank, not 21×.
- Depth (`depth_saturation.py`, 155 rows q_u ≤ 16): the cumulative area without q's over (s, r2) has rank 15 / 36 / 61
  / 86 at 1e-4 / 1e-6 / 1e-8 / 1e-10; over (s, r2, r3), 17 / 38 / 74 / 111.  A third level adds +2, +13, +25 rather
  than multiplying, and context fades ~20× per level, so the rank should saturate.
- At depth 2 → 3 with context rows (`cascade_rank.py`): ranks 37 (F), 39 (G2), ≥ 92 (G3, capped by 121 rows), the growth
  coming from the next level's dependence on the parent word.
So one automaton with K ≈ 100–150 should carry the whole satellite tree's cascade (scale and context together) to
~1e-10.  The tree sum is then one matrix-Mayer resolvent, as designed above.  What remains is learning it stably:
digit closure must see separator-visible directions, which by brute force means ~1e6 bulbs — the point where a GPU
bulb_areas pays off.

**GPU bulb areas** (`bulb.h`, `bulb.cc`, `analysis/bulb_batch.cc`, `bulb_test.cc`).  The bulb_areas algorithm as one
`__host__ __device__` routine, one thread per bulb, jobs sorted by period; CPU threads or a CUDA kernel, bit-for-bit
identical (|z| by sqrt rather than hypot, twiddles from the host, -ffp-contract=off).  Continuation is
predictor-corrector (Euler along (z', c') = (−b, a)/det from the Newton Jacobian, then Newton, splitting the step
on failure), and the Expansion<2> polish is adaptive (a second step only when the first still moves).  That is
2.3× faster on the CPU than fixed 64 radial steps and 4 substeps (24k bulbs with q ≤ 400 in 3.2 s on the laptop).
- Validation: π/16 for the 1/2 bulb to 1e-29.  An independent 50-digit Python computation of the 5/21 bulb
  (8.96139057855149651874608981050830e-06) agrees with both this code and the old bulb_areas to ~1e-30.  Across all
  24k bulbs with q ≤ 400 the two codes agree to median 6e-26, max 1.3e-23 (near-cusp bulbs, the N = 64 limit).
- Correction: the earlier observation that the old tool's 3-step polish was "short by ~1e-16" was an artifact.  The
  comparison built Decimals from printed %.17g strings of the high parts, which differ from the doubles by up to half
  an ulp; exact conversion (Decimal(float(s))) removes it.  One polish step alone really does leave ~1e-16 on
  near-parabolic bulbs.
- The area is π Σ k |a_k|² from a DFT of the polished boundary values, so the polish needs c(μ_j) only: residual in
  Expansion<2>, Jacobian from double (simplified Newton).  It agrees with the Green's-sum version to median 6e-26
  and with the 50-digit 5/21 bulb to 6e-31.
- H200, bit-identical to the CPU:

| batch | first port | predictor-corrector | + Fourier area, simplified polish |
|---|---|---|---|
| cardioid bulbs q ≤ 400 (24k) | 1.3 s | 0.53 s | 0.44 s |
| cardioid bulbs q ≤ 1000 (152k) | 9.0 s | 1.9 s | 1.3 s (117k bulbs/s) |
| cardioid bulbs q ≤ 2000 (608k) | 55 s | 10.7 s | 5.9 s (100k bulbs/s) |
| 1/3 bulb's children q ≤ 1000 (304k) | 44 s | 8.8 s | 4.8 s (63k bulbs/s) |

  The pod's 22 CPU threads take 3.1 s for the first row; the laptop's 10, 2.1 s.

## The copy layer (primitive components)

μ(M) by hyperbolic areas needs more than the satellite tree: primitive copies hold ~0.1–0.2% of the area.  Every
component factors uniquely under tuning into non-renormalizable letters: satellite letters (the cardioid's p/q
bulbs) and non-renormalizable primitive components (NRPs).  So μ(M) = S_tree + Σ_X μ(M_X) over maximal primitive
copies X = σ ⋆ P (σ a satellite path, P an NRP).  First measurements (`analysis/copies/`, from tuning's new
TUNING_DUMP: 65,242 roots through period 16 with angle words, tuning, centers, areas):
- NRP area per period decays like p^-3 (p³ × area ≈ 0.01–0.026, prime periods higher; total 8.9e-4 through period
  16, ~3e-5 beyond).  Summing all periods needs structure, not enumeration.
- Copies are far from area-conformal: κ = area(W0⋆V) A_card / (area(W0) area(V)) for primitive W0 has median 0.998
  but ranges 0.79–1.61 (the airplane copy: −14% to +61%).  Most of it is the copy cardioid's profile
  |c_X'(e^{2πiθ})|² differing from sin²(πθ); the copy's own bulbs, normalized by its own profile, are universal to
  2e-4.  So a copy's satellite-tree total is a linear functional Σ_k ĝ_k(X) T_k of its profile's Fourier
  coefficients, with universal T_k: per copy, the multiplier map only, no bulbs.
- Satellite tuning of primitive letters is wildly non-uniform (κ from 0.2 to 235), so primitive components must be
  handled where they sit, in the decorations of satellite nodes.
- External-angle words are the wrong alphabet.  g(w) = area of the component whose root receives angle w̄ is
  defined for every binary word, but its Hankel block (255 × 510, |u| ≤ 7, |s| ≤ 8) has rank 205 at 1e-8 with
  singular values decaying very slowly, for any geometric weight.  Lavaurs' pairing is non-local, and satellite
  structure (Sturmian in binary) is irregular there.
- Misiurewicz families are exactly summable.  The real primitives accumulating at −2 (words 01^{p−1}, periods 3–16)
  have area ratios → 1/256 = 4^-4 (diameter ~16^-n: both the distance to −2 and the component's scale shrink by
  the multiplier 4), with deviations 1.1e-3, 4.4e-4, 1.6e-4, …, 2.5e-9, shrinking ×0.3–0.4 per step, toward 1/4.
  That is Tan Lei's asymptotic similarity with geometric corrections.
Proposed structure: near a Misiurewicz point a, the copies at depth n come from backward orbits of the critical
point under f_a landing near a, with areas ~|Df^N|^-4.  Their total should be a Ruelle transfer operator on J(a)
with weight |f'|^-4 (pressure P(4) < 0, so geometrically convergent): a dynamical-space operator.  The
Misiurewicz points are then summed over the satellite tree (limb branch points and tips), whose automaton can
carry them.  Open: assigning each NRP to its family canonically (nested scales; on the real line, kneading).

**Ruelle operator at a Misiurewicz point: works locally, but the hard tail is parabolic** (`analysis/copies/`:
near_minus2.py, families.py, t_analytic.py, levels.py, tail_where.py, cusp_family.py; 4223 real NRPs in (−2, −1.7)).
- The standard size estimate s = 1/(βΛ²) (Λ = Π 2z_i along the center's critical orbit, β = Σ 1/Π_{j≤i} 2z_j) gives
  area(X) = (A_card/4)|s_X|² (1 + O(4^-n)): the ratio → 0.25000 as the depth n (steps the orbit lingers at β = 2)
  grows, with spread 0.03, 0.009, 0.0025, 6e-4, 1.4e-4, … from n = 2, down to 1e-6 by n = 10.
- Copies near −2 are labeled (n, w0): w0 the exit point from β, converging to a precritical point w0* of f_a = z² − 2
  (±√2 for m = 1, then 1.6629, 1.1111, …).  area(n, w0) = 256^-n |(f_a^m)'(w0*)|^-4 Φ(4^-n, w0*), and
  H(w0*) = lim Φ is smooth: a quartic in w0* fits log H to 7e-5 over the 24 converged families.  Summed over w0*,
  that is the Ruelle operator L g(z) = Σ_{f(w)=z} |f'(w)|^-4 g(w) on J(−2) = [−2, 2].
- But only deep in.  At fixed n, Φ_n is smooth in w0* only for n ≳ 7 (rms of log Φ_n about smooth fits: 1.9 at
  n = 0, 0.25 at n = 4, 1.5e-3 at n = 7, 1e-4 at n = 8).  Farther from −2, f_c differs too much from f_a.
  Extrapolating families to n = 0 in t = 4^-n is too unstable to test analyticity there.
- And near −2 is not where the copy area is.  At every period ~99% of the area in (−2, −1.7) has n = 0.  At large
  period the area instead sits in long lingering near near-parabolic cycles: 88% of the period-16 area spends ≥ 12
  steps nearly repeating a cycle of period ≤ 4.  It is one copy (c = −1.7414, area 2.0e-12), whose orbit makes 5
  passes through the period-3 gate of the airplane's cusp at −1.75: the intermittency family approaching the cusp
  from c > −1.75 (largest copy per period: d = c + 1.75 = 0.018 at p = 11, 0.010 at p = 14, 0.0086 at p = 16).
So Misiurewicz families are geometric (256^-n) and comparatively harmless.  The slowly converging (p^-3) part of the
copy area is parabolic: gate families at the cusps of primitive components (intermittency) and, by the same
mechanism, at satellite roots (seahorse-valley copies).  That is the physics of the escape-time tail and of the
satellite automaton's large-digit limits.  The copy layer should therefore attach to every component's root a
gate-family sum with a Lavaurs-universal profile in the number of passes k, carried by the automaton like F and ρ.
Next test: the real intermittency family at −1.75 out to k ~ 100 passes (periods ~300; real centers by
bracketing, areas by multiplier-map continuation), checking for a power law in k with an expansion in 1/k.

**Intermittency families at the airplane cusp: a(k) = C_w (k + σ_w)^-6 (1 + O(k^-2)), copies converge like k^-2**
(`analysis/copies/intermittency.py`, `intermittency_fit.py`; bulb.cc now takes P = 0 jobs: a component from its own
center).  Family A = kneading L(RLL)^k C (period 3k+2), B = L(RLL)^k RLC (3k+4), k = 2..120 (periods to 364,
d = c + 1.75 down to 1.4e-5), 238 components in 0.1 s on 2 CPU threads, conv ≤ 6e-23.
- d k² → π²/49 for both (universal for this cusp), and the phase σ_w = lim π/(7√d) − k is family-specific:
  σ_A = 1.04165822, σ_B = 1.42608326.  This is the Lavaurs phase of the member in the implosion coordinate.
- With that σ (from the centers alone), a(k)(k+σ)^6 = C (1 + 0/(k+σ) + c2/(k+σ)² + …): the first-order coefficient
  vanishes to the accuracy of σ (1e-4 relative), and c2/C = 0.157.  No log terms are visible down to 1e-11.
  C_A = 6.218507090694e-9, C_B = 9.787611545731e-10 (≈ 11 digits by Neville from k ≤ 120).
- Family sums: Σ_A = 3.593729090006e-8, Σ_B = 1.018759365844e-9, with tails beyond k = 120 of 5e-20 and 7e-21.
  A whole family costs a few dozen members at any target accuracy.
- The copies converge too: area(p/q child)/area(parent) → R_w(p/q) with error O(k^-2), giving 13 digits from
  k ≤ 120 (A: 1/2 0.1666965728127, 1/3 0.0238324911885, 1/4 0.0051626980869, 2/5 0.0040144451880; B: 1/2
  0.1666696296695, 1/3 0.0238275284340, …).  They are near-conformal (M's own 1/2 ratio is 1/6 exactly: deviation
  1.8e-4 for A, 1.8e-5 for B) but family-dependent: the limit copy is the Lavaurs-limit copy at phase σ_w.
- Engineering: for ~1e-10-sized parents near −1.75 the normalization w = |c_W'(λ0)|² (radial continuation in
  double) is noisy at 1e-9, so use area ratios, not F; 431 of 2618 children fail at the area stage for k ≳ 6
  (scattered).  Both need parent-local coordinates (c − c_parent with c_parent in Expansion<2>).
Parameter-space picture: near a primitive cusp c0, d = c − c0 ≈ π²/(a(k + σ)²) maps the implosion phase σ (periodic,
one unit per gate pass) onto the cusp neighbourhood with Jacobian |dd/dσ|² ∝ |k + σ|^-6.  So the cusp's copies are
the pullback of one σ-periodic Lavaurs parameter set M_L, and their area is ≈ ∫_{M_L ∩ strip} Σ_k |dd/dσ|²(k + σ) dA(σ):
a known kernel integrated over a fixed set, with corrections in 1/k.  The family index w (exit itinerary) labels
the components of M_L, and C_w ∝ their σ-areas.  This turns the p^-3 copy tail into (sum over Lavaurs components)
× (rapidly extrapolable k-sum), and each copy's tree is its limit tree plus O(k^-2).  Open: the sum over w (M_L is
itself a Mandelbrot-like set, presumably with its own tree/automaton), multi-pass words (several gate visits),
complex families (seahorse valley at satellite roots), and cusps of every primitive (the automaton must carry
σ-profiles per node).

**Seahorse-valley families at −3/4: a(k) = C_w k^-4 (1 + clean 1/k series), copies converge too**
(`analysis/copies/seahorse.py`, `seahorse_fit.py`).  In the cardioid limb k/(2k+1) (wake words (01)^{k-1}001 /
(01)^{k-1}010) the period ≤ 16 catalog shows primitive families continuing in k: S1 = (01)^{k-1}00101/00110 (period
2k+3, the limb's largest primitive), S2 = (01)^{k-1}0011 / (01)^{k-1}0100 (2k+2), S3 = (01)^k 0001/0010 (2k+4).
Centers by complex Newton from an extrapolation of 1/δ (δ = c + 3/4), continuity-checked; k = 2..120 (periods to
244, |δ| down to 0.013), 357 components in 0.1 s, plus 1785 children (q ≤ 4) in 2.7 s, no failures.  The
continuation reproduces the catalog at k = 5, 6 (not seeds).
- δ ≈ iπ/(2k + σ): the phase σ(k) = iπ/δ − 2k converges like b/k (Δσ ≈ b/k², b = 1.5, 4.3, 4.0), with no
  log k drift (a log term would make k Δσ tend to a constant; it is < 1e-4).  σ_S1 = 0.640840091 + 1.086990138i,
  σ_S2 = 1.570476147 + 0.439279045i, σ_S3 = 1.732522820 + 0.377716461i.
- area ∝ |dδ/dk|² ∝ k^-4: C = lim a k^4 = 1.72497406351e-3, 2.70222138515e-5, 1.31521485164e-5 (Neville, 11–12
  digits).  Unlike the cusp, the 1/k term of a|2k+σ|^4 does not vanish (relative +1.8, −3.5, +4.3), but the 1/k
  series is clean.  Family sums (k ≥ 2): 1.835358844693e-5, 6.712399003378e-6, 5.638923883294e-6, with tails beyond
  k = 120 of 3.3e-10, 5.0e-12, 2.5e-12 taken from the expansion.
- Copies: child/parent area ratios converge to 10–12 digits (S1: 1/2 0.1611239094, 1/3 0.0385247154, 2/3
  0.0150451474, 1/4 0.0086550192, 3/4 0.0029873262).  The limit copies are strongly non-conformal (S1's 1/3 and
  2/3 children differ by a factor 2.6; M's are equal), so each family carries its own limit tree.
Parameter-space picture, now for both kinds of parabolic point: near a root (multiplier e^{2πip/q}, here −1) the
limbs k/(2k+1) and their decorations are the image of a σ-plane Lavaurs set under δ = iπ/(2k + σ) (Jacobian ∝ k^-4);
near a primitive cusp, of a σ-periodic set under d = π²/(a(k+σ)²) (Jacobian ∝ k^-6).  Either way each Lavaurs
component w gives a family whose areas and copy trees have clean 1/k expansions with limits to ~11 digits from
k ≤ 120.  For the satellite automaton, k is the large CF digit ([0; 2, k] here), so per-node copy totals should get
the same large-digit treatment as F.  Open: the sum over w (its decay with the family's extra period j, which
needs limb catalogs past period 16, e.g. lavaurs() to 24), and joining the families to the automaton.

**The sum over seahorse families: 2^{j-1} canonical families per extra period, algebraic decay led by two-pass words**
(`angles.h`: lavaurs(P, wake) restricted to one wake, to period 29; `analysis/limb_families.cc` stats/keys/jobs/custom;
`analysis/copies/limb_constants.py`, `two_pass.py`).
- Families are canonical: stripping (01)^{k-1} from the angle words of the non-renormalizable primitives of limb
  k/(2k+1) gives keys independent of k for extra period j ≤ 2k + 1 (limb 2 differs from 3 from j = 6, 3 from 4 from
  j = 8, 4 from 5 from j = 10), and the limit count is 2^{j-1} families at extra period j.
- Centers at any k from the keys: both parameter rays traced to near the root, Newton from each, must agree (1e-9|δ|).
  Reproduces the earlier continuation (S2 = j1_0, S1 = j2_0, S3 = j3_3).  All 2046 families with j ≤ 11 (from limb 5)
  at 11 k in 12..64: 22,506 components, 0 failures, 37 s on 2 threads.  C_w by Neville in 1/k: median error 3e-9.
- Σ_w C_w per j = 2.7e-5, 1.7e-3, 2.1e-5, 1.0e-4, 5.9e-6, 1.0e-5, 1.5e-6, 1.6e-6, 4.2e-7, 3.4e-7, 3.7e-8 (j = 1..11):
  decaying but alternating, and carried by a few families (top 10% hold > 97% from j = 7).  Total through j = 11:
  1.893e-3 (phase-plane area 16ΣC/π² = 3.07e-3).
- The dominant families are two-pass words: a second run of the gate word inside the key, e.g. E0 = 00(10)^m 1 /
  00(10)^{m-1}110, O1 = 0011(01)^m / 0011(01)^{m-1}10 (the orbit turns k times at the fixed point, leaves, returns and
  turns m more times).  Followed to m = 24 (j ≈ 50, periods to ~250; 672 components): C(m) ∝ m^-6.2 (local exponents
  6.1–6.3), σ(m) drifting with Δσ ≈ 0.08–0.1/m, and pairs of sequences merge at large m (E2 and O1 to 2.41e-12 and
  2.42e-12 at m = 24, E0 and OL to 7.41e-12 and 7.40e-12): what follows the second run fades, as in the bulb tree.
- But the limits do not commute: the k → ∞ extrapolation at fixed m degrades as m grows (error 1e-5 at m = 1 to 50% at
  m ≈ 10 for E0 and OL, with k from max(12, j) to 4×).  A second, partial pass is at most a full transit (~k turns), so
  a(k, m) depends on m/k: the two-pass layer needs a two-dimensional scaling a(k, m) ≈ k^-4 m^-6 Φ(m/k), not iterated
  one-dimensional limits.
Implications: the family layer is a word problem with the same shape as the satellite one.  Keys are binary words,
gate passes are runs of 01 (the analogue of large CF digits), memory fades past a run, and runs contribute algebraic
power laws in their lengths.  So C_w should be learnable as a weighted automaton over run-length-encoded keys, with
large-run asymptotics (m^-6 per extra pass) and the joint scaling in m/k at the outer pass.  Next: a(k, m) on a 2D grid
(m up to k ~ 100) for the scaling function Φ, and a run-length Hankel test on the census constants.

**Cluster census (j ≤ 14) and two-pass grid; family constants are a size-estimate sum** (`cluster_families.sh` on the
cpu queue: 147,447 census components and 24,448 grid components in 61 s on 30 cores, 1 failure; `family_size.py`,
`family_wfa.py`, `two_pass_grid.py`).
- Kneading is the right label: every family's kneading sequence is 1^{2k} 0 t, and the 2^{j-1} families of extra
  period j are exactly the 2^{j-1} binary tails t of length j - 1 (for |t| < 2k; t = 1^{2k} is the bulb's period
  doubling).  So family constants form a function G(t) on the free binary monoid, with no invalid words.
- Σ_w C_w by j (limb 7, all 16,383 families with j ≤ 14): … 3.4e-7, 1.3e-7, 9.7e-8, 4.7e-8, 3.3e-8 (j = 10..14), i.e.
  ∝ j^-7 (local exponent 7.0 over 10→12 and 12→14), all of it in the top 10% of families.  Brute enumeration
  by j cannot reach 1e-12 (tail ~ J^-6); the long-run families must be summed structurally.
- Spectral WFA on raw G over kneading bits (block |u| ≤ 5, |s| ≤ 6, tested on all 8192 tails of length 13): the
  Hankel singular values fall geometrically (1e-14 by rank 40) but predictions are poor (total off by 1–6%, single
  tails by orders of magnitude).  99.5% of the length-13 mass is in the 144 tails with a run ≥ 9, which is beyond
  what the block saw (power laws in run length have no finite rank), and the 1e-3..1e-30 dynamic range swamps the
  rest.
- But the constants are almost entirely the standard size estimate: F = area/(A_card |s|²) with s = 1/(βΛ²) from the
  center's critical orbit tends to 1 in k: F_∞ median 0.99999967, 80% within 1.6e-4 of 1, 98% within 1.2e-3; the
  outliers are the few large families (S1: 1.077; the leading E0, O1, E2 members: 0.97–1.06).  So Σ_w C_w ≈ A_card Σ_w
  lim k^4 |β_w Λ_w²|^-2, a weighted sum over critical orbits of the limiting (Lavaurs) dynamics: transfer-operator
  form with weight |Λ|^-4, plus small corrections F - 1 concentrated in a handful of large families.
- Two-pass grid a(k, m), m < k ≤ 128, m ≤ 64: g = a k^4 m^6 collapses in t = m/k only roughly (E2 at t = 1/4:
  5.54e-4, 5.41e-4, 5.24e-4 for m = 8, 16, 32; other columns drift 10–40% per doubling).  Consistent with
  a ≈ k^-4 m^-6 Φ(m/k) (a partial pass of m turns starts ~m^{-1/2} from the degenerate fixed point, derivative
  ∝ m^{3/2}; the throat width set by δ ∝ 1/k enters through m/k) with sizable 1/m corrections; m ≤ 64 is too
  small to pin Φ.
Next: compute families' constants from the dynamics (the size estimate in the Lavaurs limit) rather than from
components, which makes long runs cheap, and learn the small correction F - 1 (O(1e-4), smooth) by WFA; get Φ(t)
from the model map z + a z³ + perturbation instead of from M.

**Family areas split into an exact orbit factor and a learnable O(1) correction; Φ(m/k) from the gate model**
(`limb_families size`, `family_size.py`, `family_wfa_F.py`, `two_pass_grid.py`).
- Size-estimate engine (`limb_families size`): centers by rays only at a family's first three k, then Newton
  continuation in k from a quadratic extrapolation of 1/δ (continuity-checked, rays as fallback), and |s|², |Λ|, |β|
  from the orbit; no boundary tracing.  On the two-pass grid it reproduces the centers to 1e-14 |δ| with 24 ray
  pairs for 864 points.  Limit: centers are doubles, so components below ~1e-13 across (two-pass m ≳ 100–160)
  need Newton on δ = c + 3/4 or Expansion<2>.
- WFA on h = F_t(k) − 1 at fixed k (exact, no k-extrapolation noise; |h| median 1.8e-3, range −0.035..0.080):
  trained on tails ≤ 12, tested on all 8192 tails of length 13, the error in F is median 3.9e-8, 99% 1.5e-6, max
  3.4e-5 at rank 32 (k = 64; k = 32 alike), against 1–6% total error for the raw constants.  So
  area = A_card |s|² F with |s|² exact from the critical orbit and F learnable to ~1e-6: the remaining hard part
  is the sum of |s|² over families, which is where the dynamics (transfer operator) has to come in.
- Gate model for Φ: near α = −1/2 with δ = iπ/(2k), f² is a rotation by θ ≈ π/k plus a cubic term, linear in
  u = 1/w²: u ↦ e^{2iθ}(u + 2).  A partial pass of m turns ending at the exit starts at
  |u_0| = 2|sin mθ/sin θ| ≈ (2k/π) sin πt and contributes |u_0|^{3/2} to Λ (β is unaffected), so
  area ∝ m^-6 (πt/sin πt)^6.  Checks: the landing depth max|z+1/2|^-2 of the second pass follows (sin πt)
  (ratios 0.724, 0.930 between t = 1/2, 1/3, 1/4 vs 0.707, 0.918); with that factor g = a k^4 m^6 is flat in t for
  E0 and OL (m = 32: 7.82e-4, 8.40e-4, 8.41e-4 at t = 1/2, 1/3, 1/4).  E2 and O1 land deeper than the sine law at
  finite t (their Λ^-4 relative to E0 runs 0.05 → 0.34 from t = 0.8 to 0.25, while β agrees to 3–8%), so the
  landing law depends on the family; and at fixed t all sequences drift by a further ~m^-0.6 over m = 4..32 that
  the model lacks (log-drifting phase σ(m)?).  Larger m needs the δ-coordinate centers above.

**Deep areas and the family pipeline** (`bulb.cc` local jobs, `family_jobs.py`, `family_constants.py`, `skeletons.py`,
`family_sum.py`, `cluster_pipe.sh`).
- `bulb_job_local`: a P = 0 job from a double-double center, tracked in c = c0 + L Δ with L = |1/(βΛ²)| and the whole
  continuation in Expansion<2>.  Double tracking cannot work there: a double periodic point carries ~1e-16 n |Λ|²
  relative error (the ordinary path survives small components only because its tolerance is absolute in c).
  Ordinary jobs bitwise unchanged; local agrees with ordinary to 3e-24 where both run; 428 two-pass components to
  m = 1024, k = 16781 (areas to 4e-39, period 35,611) all succeed in 27 s on 2 threads.
- The size estimate in the k → ∞ limit is excellent for long runs: F_∞ − 1 = +2.7e-4, +3.1e-5, +3.8e-6 (E0) and
  −2.4e-4, −3.0e-5, −3.7e-6 (E2) at m = 64, 256, 1024, i.e. ∝ m^{-3/2}.
- Pipeline: digit word → keys_for → size engine (centers) → bulb_batch --local (areas) → C = lim k^4 area (polynomial
  in t = n/2k with the gate shape divided out).  Pilot class ('1', 3, *) even: matches the census to 4e-9 relative for
  its heaviest word (3e-6 at n = 10); C n^6 still falls at n = 40 (local exponent −0.46), so classes need large n.
- Every word belongs to one class (skeleton = first maximal digit as '*', plus its parity); the 60 heaviest classes
  hold 99.76% of the census mass at digit sum ≥ 9.  Campaign: 2546 words (n ≤ 64 dense, sparse to 512), 163,564
  points, on the cluster's cpu queue.
- Campaign 1 (60 classes, 2546 words, 163,564 points, 8 min on 30 cores, 0 failures) agrees with the census on 206
  overlapping words (median 1e-7, limited by the census's k ≤ 64).  The 100 heaviest families recomputed with
  k ≤ 1024 (Neville in 1/k, subsets agree to ~1e-20): S1 = 1.724974063498964762e-3 (the k ≤ 120 fit: 1.72497406351e-3);
  census Σ weighted error 2.8e-10 → 2.4e-13.
- Σ_w C_w ≈ 1.893452650832e-3: census 1.893391954437e-3, beyond it 6.0548e-8 from the 60 classes (computed to
  n = 64, then fitted tails ~1e-11), other classes ≈ 1.5e-10 (crude: their census share 2.45e-3).  Phase-plane area
  16Σ/π² = 3.06955e-3.  Campaign 2 (classes 61–800 to n = 48, the four dominant classes refined) under way.
- Long runs decouple words: the classes' large-n limits come in groups of four (the entry contexts '0 1', '1 1'
  even/odd, '1' are asymptotically equivalent) and factorize: the four leading classes → 64·3.44e-4 = 0.0220 (the
  universal second pass), a digit 2 after the run multiplies by ≈ 5.0e-3, a 4 by ≈ 1.75e-4, prefixes ('1 3',
  '1 2 1', '1 1 2') by ≈ 3.1e-3 and ('0 2', '0 1 1', '1 1 1') by ≈ 1.75e-3.  So asymptotically C(u, n, v) ≈
  L(u) 0.022 n^-6 R(v): the large-digit matrix of a digit automaton is rank one, and sums over words with a long
  digit factor into prefix sum × n-sum × suffix sum.
- All of M in run-length digits is not geometric (area share by number of digits 0.67, 0.30, 3.9e-3, 1.7e-2, 9e-4,
  6.5e-3, … through period 16; period-doubling-type satellites need many short runs): kneading runs are the right
  coordinates at a parabolic point, the satellite tree keeps its continued-fraction / tuning coordinates.
- Campaign 2 (classes 61–800 to n = 48 plus n ≈ 64, and the four dominant classes with k ≤ 1024: 17,072 words,
  764,080 points, 32 min on 30 cores; 1.67M ray-pair fallbacks dominate the cost, so the small-k predictor needs
  work; 247 families stopped early, 16 words skipped for too few large-k points, 695 without a census base for keys).
  **Σ_w C_w = 1.8934526215e-3 ± 4e-12** (phase-plane area 16Σ/π² = 3.06954977e-3): census 1.893391954437e-3 (±2.4e-13),
  758 classes beyond it 6.066641e-8 (extrapolation ±2.9e-12, fitted tails ±1.9e-13), classes past 800 ≈ 6.6e-13
  (census-share estimate; the 60-class version of it, 1.5e-10, was 25% high), words with only small digits beyond the
  census ~1e-13 (geometric, ×0.15–0.2 per unit digit sum), skipped words ≲ 1e-12.  So an exponentially large family
  layer sums to ~2e-9 relative from ~0.9M components, through its structure: kneading digits, classes, and n^{-1/2}
  expansions with factorized limits.
- Caveat: Σ C_w is the k → ∞ (phase-plane) constant, not an area of M.  The families' actual area is
  Σ_w Σ_{k ≥ k_w} a(k, w) with a(k, w) far below C_w k^-4 at small k (S1 at k = 2: 7.75e-6 vs C/16 = 1.1e-4), families
  exist only from k_w ≈ (j−1)/2, and small limbs' longer tails belong to other parabolic points.  Assembling μ(M)
  means attaching such family layers to the satellite tree's nodes (roots and cusps), with the small-k region in the
  tree itself.

**The 1/2-root layer is a series over full gate transits** (`limb_families limb`, `diag_census.py`, `diagonal.py`,
`diag_negative.py`, `diag_sum.py`; lavaurs/maximal_tuning now to period 62).
- Assembly target: μ(M) = A_card + Σ_{W ∈ NR} a(W) R_W (NR = cardioid bulbs and all NRPs, R_W the copy factor); the
  1/2 root's new layer is the NRPs of the cardioid limbs k/(2k+1) and their mirrors: T(K) = 2 Σ_{k≥K} A(k).
- Beyond the stable range (tail > 2k) limbs are not just "more of the same": the unstable band has a k-independent
  structure in the offset j - 2k (missing tails and duplicated kneading sequences, identical for k = 2..5) and real
  mass: by number of tail runs of transit length (≥ 2k-2), limb 3: 75.7 / 23.6 / 0.6%, limb 4: 80.7 / 18.9 / 0.4%,
  limb 5: 84.4 / 15.6 / ~0% (truncated at j ≤ 2k+8).  These are orbits with a second full transit, e.g. tail
  1^(2k+2) (limb 4: area 2.3e-7, the third largest NRP).  A pass can be longer than the first.
- Their keys are k-independent up to one 2-bit insertion per limb (14,123 labels (u, form, d, v), tail u R v with R a
  run of 2k + d; limb 6 predicted exactly), so the constants D = lim a(k, label) k^4 follow from limbs 4 and 5.
  D(-, 1^(2k+d), -): d = 2: 1.858e-4, -2: 1.366e-5, 4: 1.317e-5, 6: 1.9e-6, odd d ~1e-7..1e-6; at negative d the words
  are stable (keys_for) and D ~ |d|^-6 (D|d|^6 ≈ 0.012–0.016 for d = -10..-14).
- Two-transit sector Σ D ≈ 2.2729e-4 (14,121 limb labels 2.254362e-4, further negative offsets 1.85e-6, |d|^-6 tail
  4e-9): +12% on top of the one-transit S = 1.8934526e-3.  It is truncated in |u| + |v| (offsets d ≥ 2 missing past
  |u| + |v| = 4; masses 2.18e-4, 6.7e-6, ≥4.4e-7, ≥1.4e-7 at |u| + |v| = 0, 2, 4, 6), and unlike long partial passes a
  full transit does not decouple the word (D(u, d, -)/D(-, d, -) varies ×5 over d).
- So the k → ∞ constant of A(k) k^4 is S_0 + S_1 + S_2 + … over the number of full transits (S_1/S_0 ≈ 0.12, S_2 a
  few % of S_1): a renewal over iterations of the Lavaurs map through the gate, with excursion words between
  transits.  The natural exact object is the Lavaurs-phase transfer operator; computing via M at large k works but
  needs each sector's own word census.

**Lavaurs model at the 1/2 root: phase-plane areas computed directly** (`analysis/lavaurs/`: fatou_series.py,
lavaurs.py, centers.py, calibrate.py, model_area.py; Python prototype in double precision).
- w = z + 1/2, F(w) = -w + w² (multiplier -1), F_δ = F + δ.  Formal Fatou coordinate with Φ(F(w)) = Φ(w) + 1/2:
  Φ = 1/(4w²) + 1/(4w) + (11/8) L(w) + Σ a_j w^j (a_1 = -5/16, a_2 = 75/64, …; exact rationals, Borel-type growth
  |a_j|^{1/j} ≈ 1.7 → 2.9 for j = 10 → 60; L = log on one petal of each pair, log(-·) on its partner).  The odd 1/(4w)
  term is the source of the m^{-1/2} expansions seen in M.  Φ_a by iteration into the attracting petal (|w| < 0.05,
  30 terms); F swaps the repelling petals, so two parametrizations: Ψ₊ (local inverse in the upper petal, then F^{2m})
  and Ψ₋(ζ) = F(Ψ₊(ζ - 1/2)).  Lavaurs map g_σ = Ψ_±(Φ_a + σ), the exit petal equal to the entering one.
- Calibration against M (point-cloud match of single-transit centers with the census phases, 34 of the top 40 at 5e-9):
  σ_model = -σ_M/2 + 3πi/8 + 1/2 with excursion length n = j - 1 (δ = iπ/(2k + σ_M)), hence C_w = (π²/4) area_model.
- Component areas by multiplier-map boundary tracing in σ (2-jets through Φ_a and Ψ, N = 64, DFT) reproduce the M
  constants: S1 1.724974063499687e-3 (M: …498965e-3, +4e-13), S2 +5e-11, S3 -3e-11, E0m2 +1e-11, E2m1 -9e-12; the
  limb's bulb itself appears as the component 0.2081.  So every sector's constants are areas of components of one
  fixed σ-plane set, computable without limbs, k-extrapolation or ray tracing.
- Since g_σ commutes with F, a cycle F^n g F^{n'} g F … depends only on the total n and the number of transits r: an
  r-transit component is σ with g_σ^{r-1}(v) in the basin and F^n(g_σ^r(v)) = 1/2.  Excursions between transits enter
  only through which basin component g_σ(·) lands in.
- C++ model (`lavaurs.h/.cc`, `lavaurs_test`, `analysis/lavaurs/lavaurs_area.cc`): double-double, ~0.3 s per component,
  C(S1) = 1.7249740634989648e-3 (M: 1.724974063498964762e-3).  All 16,383 census families on the cluster (6 min, 0
  failures, all centers distinct): median agreement 1.6e-8 with the 1/k-extrapolated census (its own error), Σ
  differing by 2.5e-14.  Two-transit labels: rule n = |u| + |v| + d - 1 at the standard seed (r transits: σ - 1/2 ↔
  n + r); 14,129 of 14,505 match the M-side D (median 3.4e-7); the rest are poor seeds.
- M-side seeds fail for long digits (σ_M extrapolation in 1/k has t = n/2k corrections), so `lavaurs_area --walk`
  walks a class in its long digit, predicting centers by Lagrange extrapolation in n^{-1/2} (a long digit is a pass
  landing ~n^{-1/2} from the parabolic point) and accepting within 0.3 radii; seeded with the reliable model centers
  (census + campaign words whose model C matches), all top classes run to n ≈ 130.
- **S_0 = 1.8934526216147e-3** on model constants (census 1.893391954412168e-3, beyond it 6.0666442e-8 with tail
  spread 2.4e-14, classes past 800 ≈ 6.6e-13, small-digit words ≈ 1e-13): ± ~5e-13, agreeing with the M-based value.
- Island structure for the transit series: at a single-transit center σ_u the horn map H has a critical point, and
  Θ(σ) = H(ζ0 + σ) + σ - ζ0 satisfies Θ(σ_u) = σ_u - (n_u + 1)/2, Θ'(σ_u) = 1; two-transit centers are exactly
  Θ^{-1}(single-transit centers), two per island.  The largest, D2 (-|1|2|-), lies in S1's island over the limb's
  bulb: Θ(D2) = the bulb's center.  Targets include tuned components (bulb, satellites), and island u over target u
  is u's period doubling, so a model-native census of the NRP sector needs a tuning filter.

**Rigor plan (discussed 2026-10-09)**: layers and what each can be.
- L0: Σ_hyperbolic ≤ μ(M) is a theorem; equality needs no queer components (MLC) and μ(∂M) = 0 (open), so the
  rigorous end products are a certified lower bound and a conditional "μ(M) ∈ [L, L + T]".
- L1: the countable index set (tuning tree, NRPs per limb, kneading digits, transit counts, island × target pairs) is
  rigorous combinatorics in principle; our empirical rules (stable-range kneading bijection, 2-bit insertions, the
  universal unstable band, island preimages) need proofs, and the combinatorics-to-component bridge (ray landing,
  tuning, Lavaurs' theorem) would enter Lean as axioms.
- L2: each component's area can be certified (Krawczyk on |μ| = r > 1, Cauchy tail bound; model components need
  certified Fatou coordinates).
- L3: tails and accelerations (1/k, n^{-1/2} expansions, Neville, automata, class and transit truncations) are
  heuristic; rigorous tail bounds via Koebe distortion and fading-memory estimates are plausible at modest precision.
- L4–L5: ball arithmetic and Rust + Aeneas → Lean for the kernels are feasible; CUDA stays a bitwise-checked mirror.
- Practice from now on: every computed component carries a combinatorial label (M kneading tails / digit words for
  one transit; for two, the constituent labels with branch side, shift j and the horn map's critical point), so a
  later certified pass reuses the enumeration unchanged.

**Two-transit census in the model (islands)**: two-transit centers solve Θ(σ) = σ_c + j/2 for single-transit
centers σ_c; tracked from a source center by Euler-Newton path lifting with adaptive steps.  Island S1 over the
limb's bulb gives D2 (L branch) and 10|1|0|- (R branch) exactly.  Over the top 30 × 30 pairs: 264 distinct
components; 28 are doublings (source u over target u: u's period-doubling satellite, tuned, absent from the NRP
labels, mass 4.5e-5); of the 236 others 99 are M labels and 137 (7.0e-7) were missing from the label census (most
likely the 166 ambiguous labels skipped there).  The immediate basin's Ψ-preimage is one region holding a whole
sequence of single-transit centers (S1, j4_2, j6_10, j8_42, …), where H has many critical points; so components are
deduplicated by center and labeled by the critical point they reach, not by the source they were tracked from.
- Doublings and shifts: a family u's period-doubling satellite is the two-transit component targeting u itself at
  shift j = -(n_u + 1), attached to u (1.6 radii from its center, area ratio ≈ 0.16 like any 1/2 satellite); it is
  found from many source regions, so it is identified geometrically (within 2.5 radii of its target's center).  The
  target ranges over all representatives σ_c - m (period 1 in σ, n + 2m): shifts j ≤ -2 are second excursions with a
  long first run, decaying like a long digit (bulb target: 1.3e-5, 1.9e-6, 4.0e-7 at j = -2, -4, -6).
- S_1 (union of the validated label components and the island censuses (top 30 × 30 with j ≥ -14, 38,883 pairs
  with j ≥ -1), deduplicated, 12 doublings removed): **S_1 = 2.6737e-4** (was 2.273e-4 from labels alone: the label
  census had dropped "ambiguous" labels, e.g. -|1|0|01 worth 3.9e-5).  Island-found mass converges fast in source and
  target rank (rank < 10 already within 5e-8 at j ≥ -1).
- r transits: Θ_r(σ) = p_r - ζ0 with p_{i+1} = H(p_i) + σ; (r-1)-transit centers are the tracking sources (Θ_r' = 1
  there).  D2's region over the bulb (three transits): C = 3.336e-5 = 0.18 D2.  The transit series decays slower
  than first estimated (S_1/S_0 ≈ 0.14; chain S1 → D2 → D3 ratios 0.11, 0.18).

**Tunings and the cusp classifier (2026-10-09 overnight).**  The multiplier map's cusp measure |σ'(1)|/|a_1| (~1e-17 at
a primitive cusp, O(1) at a satellite root) separates primitive components from satellites (doublings, the bulb's
children, the cardioid's bulbs), but a cusp-primitive r-transit component can still be tuned: U*X with U an r'-transit
component and X primitive of period r/r' in M.  Its cycle is (F^n g^r F)^p: r p transits, excursion p n + p - 1.
`lavaurs_area` with LAVAURS_TUNE computes U*X from U's multiplier map (σ ≈ center + 2 a_1 c_X, refined by
interpolating the copy map through the tunings already found, linear fallback for distorted satellite copies).
- Top 1000 single-transit families + the bulb: the airplane and the period-4 primitive (c = -1.9408) tunings all found,
  all cusp-primitive, no duplicates; ratios area(U*A)/area(U): median 3.450e-4 (S1: 3.334e-4; bulb: 7.9e-4), for
  -1.9408: 9.74e-7.  Σ C(U*A) over NRP families ≈ 6.3e-7, which must leave the r = 3 sector (it is in R_U).
- Copies are far from area-conformal in σ: doubling ratios 0.161 (S1), 0.184 (j4_0), 0.153 (j4_2) vs 1/6 in M.  So
  R_W at the 1e-6 level needs per-copy evaluation (profile functional plus primitive sub-copies), not μ(M)/A_card.
- The D-chain over the bulb (r = 2..7: 1.858e-4, 3.336e-5, 2.834e-6, 2.42e-8, 2.20e-9, 2.94e-11) converges to
  σ ≈ -0.033 + 0.601i with irregular ratios; a single chain says little about S_{r+1}/S_r (branching dominates).

**Finite k.**  Per limb, the census (16,383 one-transit families, complete at k = 16..64) and the two-transit labels
(13,844 present at all k = 16..96) give G_r(k) = k^4 Σ area, fitted as model constant + Σ b_j k^-j:
- one transit: b1..3 = 6.957e-4, -2.846e-2, -0.159 (b1/S = 0.37), Σ_{k≥16} G/k^4 = 1.6458005e-7, Σ_{k≥32} =
  2.0152846e-8, Σ_{k≥64} = 2.46939912e-9 (fit orders 6-8 agree to ~1e-15).
- two-transit labels (model D = 2.2475848e-4): b1 = 3.528e-4 (b1/D = 1.57), b2 = -5.99e-3; Σ_{k≥16} = 2.0230606e-8,
  Σ_{k≥64} = 2.9678829e-10.
- The Jacobian of δ = c + 3/4 = iπ/(2k + σ_M) explains only part of the 1/k term; the rest is family-specific, and
  the phases drift as (σ_M(k) - σ_M∞) k = -1.50-0.16i, -1.98+0.34i, -1.91-0.64i at k = 64 (j2_0, j4_0, j4_2), still
  moving (log k / k terms from the Fatou coordinate's log).  So model centers do not seed M at small k for small
  components; island-only components have no M data yet.
- Error budget for 1e-12 in μ(M): at K = 64 the tail needs S_tot to ~3e-7 and b1 to ~10%; at K = 16 S_tot to ~5e-9,
  b1 to ~1e-7 absolute, b2 to ~2e-6.  K = 64 is the natural cut for this layer; limbs 16..63 need M-side data for
  every component above ~1e-9 (or a next-order Lavaurs model).

**Sectors after cusp classification, and the Farey limbs in the strip.**  All 23,372 distinct two-transit centers
(labels + island censuses) and 9,887 three-transit centers (island3: top 300 two-transit sources × bulb + top 30
families, shifts -4..1) rerun with the cusp measure (cluster cpu queue, 27 min); the cusp histogram has an empty gap
from 1e-12 to 1e-1, so the classification is unambiguous.
- r = 2: **S_1 = 2.93675e-4** over 23,359 primitive components (labels 2.2842e-4, island-only 6.53e-5); 12
  satellites (the cardioid bulb B_2 = 1.5065e-2, S1's doubling 2.779e-4, j4_2's, j1_0's, a chain of doublings of
  the j(2m+1)_(4^m-1) families).  Island-found mass is complete in the shift j (geometric decay past j = -6).
- r = 3 (partial census): 5.0805e-5 primitive and not tuned (93% bulb-targeted; D3 = 3.336e-5 dominates), 3
  satellites (S1's 1/3 satellite 2.6e-5 among them).  Ratios S_1/S_0 = 0.155, S_2/S_1 ≥ 0.17.
- The σ strip holds more than the limb k/(2k+1).  Rotation numbers between θ_k = k/(2k+1) and θ_{k+1} belong to
  Farey descendants B_a (the mediant (2k+1)/(4k+4) is the cardioid bulb B_2 at -0.5009 + 0.0418i, area ratio to the
  bulb 0.072 ≈ (q_1/q_2)^4), each with its own limb of NRPs.  The bulb-source island census shows the mediant
  limb's one-transit layer: bulb~t|R|-6 for single-transit families t, with C ≈ (0.005–0.013) C_t (j2_0: 0.0126,
  j4_0: 0.0083, j4_2: 0.0093, j1_0: 0.0047), a non-conformal scaled copy of the census at r = 2.
- So the transit count mixes limbs: a Farey limb at level a contributes at ~a times the transits, and S_r gets a
  power-law part (~r^-3 from φ(m) m^-4 sizes) on top of the bulb limb's own series.  Bookkeeping: the strip
  computes all limbs with rotation numbers in [θ_k, θ_{k+1}), which is the consistent definition of this layer
  (the satellite/tree side then handles rotation numbers < θ_K); the sum should be organized per limb (the
  bulb's layer plus Σ_F ρ_F × copy, ρ_mediant ≈ 0.0115 for the one-transit part), not by raw transit count.
- Open: a wake classifier in the model (which limb a component's α-rotation puts it in); the M-side labels are
  bulb-limb by construction, and limb-5 σ_M matching is too noisy (finite-k drift O(1) at k = 5).
- Confirmed in M (`limb_families wake J p q`: every NRP of the p/q limb, full words): the mediant limb 9/20 at k = 4
  has one dominant NRP, period 22 = q + 2, σ_M = 1.713 + 1.972i (model bulb~j2_0|R|-6 in the representative with
  n = 3: 1.883 + 1.941i; drift as expected at k = 4), a k^4 = 6.9e-6 (model C = 2.18e-5), the rest ≤ 1.9e-7.
- Structure: the Farey fraction between θ_k and θ_{k+1} with weights (a, b) has rotation number exactly θ_x =
  x/(2x+1) at x = k + b/(a+b).  So the strip is a continuum of "fractional limbs" t = x - k ∈ [0, 1), with limbs at
  rationals t = b/m, bulbs scaled ~m^-4 (B_2 at t = 1/2: 0.072 ≈ 2^-4 · 1.15), and a gate passage costing m
  transits.  If each limb's layer is κ m^-4 times the bulb limb's, the Farey part is κ Σ_{m≥2} φ(m) m^-4 =
  κ (ζ(3)/ζ(4) - 1) = 0.11 κ; the mediant's one-transit copy suggests κ ≈ 0.18, i.e. ~2% of S_tot (~4e-5),
  converging like M^-2 in the denominator cutoff.
- Bulb-source census (island5: bulb × all 16,383 census families, shifts -10..1, both branches): shifts j and j - 4
  on the two branches are the same component one period apart, leaving three full copies of the one-transit census
  (mod 1): A (R -6 ≡ L -2) Σ = 2.3008e-5 = 0.01215 S_0 at Re σ ≈ 0.56 next to B_2: the mediant limb's one-transit
  layer (per-family ratio median 0.0036, max 0.0126: far from uniform); B (R -8 ≡ L -4) Σ = 4.10e-5 = 0.0217 S_0 at
  Re σ ≈ 0.27: the bulb limb's own two-transit family "-|1|d|tail" (its j2_0 member is the validated -|1|0|01);
  C (L 0) Σ = 3.3e-6, bulb-limb labels.  New mass only 4.41e-7 (R -10).  **S_1 = 2.94116e-4** in the strip
  (bulb limb 2.711e-4, mediant limb's one-transit copy 2.301e-5).
- Limb leading NRPs at k = 8 (a k^4): 8/17: 1.317e-3, mediant 17/36: 1.254e-5 (ratio 9.5e-3; model 0.0126), m = 3
  limbs 25/53: 9.8e-7, 26/55: 4.8e-7 (ratios to the mediant's 0.078, 0.038 vs (2/3)^4 = 0.20): the Farey part falls
  faster than m^-4 at this k, roughly 1.15 × the mediant's, ~1.3% of S_tot.
- Tail now (tail_sum.py, R left out, both mirror sides): 2Σ_{k≥16} = 4.0046e-7 ± 5.8e-9, 2Σ_{k≥64} = 5.958e-9 ± 6e-11,
  the error from unseen transit mass (≥ 3 transits, estimated 3e-5 ± 2e-5); census and label pieces are exact to
  1e-15 from M.  For 1e-12 at K = 64 the unseen mass must be known to ~7e-7.
- Satellites must be island sources too.  B_2 (the mediant cardioid bulb, r = 2) as a three-transit source gives the
  m = 3 Farey bulbs (C = 3.06e-3 at Re σ = -1/3, -2/3: t = 1/3, 2/3; ratio to the bulb 0.0147 ≈ 3^-4 · 1.19), the
  m = 3 limbs' leading NRPs (B2~j2_0|L|-2: 1.52e-6 at -0.302 + 0.102i, |L|-4: 9.6e-7 at -0.645 + 0.101i, matching
  the M limbs 25/53 and 26/55 at k = 8), and a large primitive B2~bulb|R|-8 = 6.91e-5 in the mediant region: S_2 is
  well above the 5.08e-5 found from primitive sources.  Doubling regions (S1*H) give only small NRPs (≤ 1.3e-7).
- Each Farey limb has a heavy NRP in its bulb's region over the limb's bulb B_1, plus copies of the one-transit
  census over all targets (left branch, even shifts).  Mediant (m = 2, B_2 source, r = 3): B2~bulb|R|-8 = 6.91e-5
  (R shifts -10, -12: 6.0e-7, 3.8e-8; in M: limb 9/20's largest NRP, P = 29, a k^4 = 1.05e-5 at k = 4, larger than
  its period-22 copy-A member), copies L -6, -2, -4, 0: 3.2e-6, 1.6e-6, 1.0e-6, 1.6e-7 (partial over targets).
  m = 3 (B_3 sources, r = 4): 1.58e-5 + 7.1e-6 + 4.5e-6 = 2.7e-5 = 0.39 × the mediant's (2 (2/3)^4 = 0.40), and the
  m = 4 Farey bulbs at Re σ = -1/4, -3/4 (C = 9.6e-4 = 4^-4 · 1.19 × bulb, the same 1.19 as m = 2, 3).
- So per Farey denominator m the heavy pieces scale ~φ(m) m^-4 (heavy part alone ≈ 6.9e-5 · 16 · (ζ(3)/ζ(4) - 1) ≈
  1.2e-4), and the strip's transit series is not geometric (S_2 ≈ 1.47e-4 > S_1/2): r-transit mass is dominated by
  Farey limbs of denominator ~r - 1.  The bulb limb's own series (S_0, 2.711e-4, …) looks geometric (~0.15).
  Organize: strip = bulb limb (transit series) + Σ_F limb F, each limb F ≈ a distorted copy, with a scaling law in m
  (bulbs: 1.19 m^-4 exactly in the pattern so far) and a tail ~M^-2 in the denominator cutoff.
- Farey edge t = 1/m (bulb B_m as an r = m source over the limb's bulb at r = m + 1 gives B_{m+1} and limb m's heavy
  NRPs): bulbs C = 1.5065e-2, 3.061e-3, 9.645e-4, 3.896e-4, 1.849e-4, 9.83e-5, 5.68e-5 for m = 2..8, i.e. m^4 C =
  0.241, 0.248, 0.247, 0.244, 0.240, 0.236, 0.233 (the bulb itself: 0.208); the m = 4 step also gives the t = 2/5
  bulb (4.09e-4 at -0.400).  Heaviest NRP from B_m's region: 6.91e-5, 1.58e-5, 1.19e-6, 8.06e-7, 1.85e-7, 6.3e-8 for
  m = 2..7, ratios to the bulb 4.6e-3, 5.2e-3, 1.2e-3, 2.1e-3, 1.0e-3, 6.4e-4: falling faster than the bulbs (~m^-6),
  so the Farey part converges faster than φ(m) m^-4 would suggest, but these are single components, not limb
  totals.
- Three-transit sector (island3 + island6's top 10,000 of 91,700 importance-ordered pairs + island7: 12 two-transit
  satellites as sources, B_2 × all 16,384 targets, shifts -12..1): **S_2 ≥ 1.4226e-4** (5.08e-5 + 1.90e-5 + 7.25e-5;
  airplane tunings 4.8e-9 removed), satellites 8.17e-3 (the m = 3 Farey bulbs 2 × 3.06e-3, the bulb's 1/3 child).
- Tail (tail_sum.py with island-only r = 2 mass 6.57e-5, r = 3 1.4226e-4, unseen ≥ 4 transits 5e-5 ± 3e-5):
  S_tot ≈ 2.3905e-3; 2Σ_{k≥16} = 4.2196e-7 ± 9.1e-9, 2Σ_{k≥32} = 5.1218e-8 ± 8.6e-10, 2Σ_{k≥64} = 6.2545e-9 ±
  9.1e-11 (before the copy factor R; × 1.279 ≈ 5.40e-7, 6.55e-8, 8.00e-9).  Census and label pieces are exact to
  1e-15 from M; the error is the unseen transit mass and the 1/k coefficient of the model-only mass.
- Four transits (island8, r = 4: the two m = 3 Farey bulbs × top 2000 targets, other r = 3 satellites and top 300
  r = 3 primitives × bulb + top 30; only the top ~2000 of 13,519 pairs survive, the Job hit its activeDeadline):
  S_3 ≥ 6.08e-5 (the m = 3 limbs' heavy pieces 1.58e-5, 7.1e-6, 4.5e-6; a new 1.40e-5 at -1.214 + 0.246i).  Sector
  sequence in the strip: S_0..S_3 = 1.8935e-3, 2.941e-4, ≥1.423e-4, ≥6.08e-5 (ratios 0.155, 0.48, ≥0.43): beyond two
  transits the Farey limbs dominate and the raw transit series decays slowly; the rerun (island9) covers the rest.
- Tail with these (model-only mass 2.688e-4, unseen 6e-5 ± 4e-5): S_tot ≈ 2.461e-3; 2Σ_{k≥16} = 4.356e-7 ± 1.2e-8,
  2Σ_{k≥64} = 6.442e-9 ± 1.2e-10 before R.  1e-12 at K = 64 needs the unseen mass to ~7e-7: the Farey limbs must be
  summed by structure (per-limb layers with a scaling law in the denominator m), not by raw transit censuses.

**Summing by Farey limb (2026-10-09, day).**
- Limb classifier (`limb_class.py`, `kneading_side.py`): the first internal-address step q = position of the first
  kneading 0; in the representative convention x_i sits at (2k+1) i, so the first 0 at offset 2b - 1 from x_m means the
  limb t = b/m (q = m(2k+1) + 2b).  The kneading partition near 1/2 is R_{1/6} ∪ A ∪ R_{2/3} where the arc A through the
  critical basin component is not the real segment: at limb parameters the partition is the preimage of the critical
  value's ray, which near v is the shifted parameter ray rising from the component's root, so A is the lift under Φ_a
  (2:1 at the critical point) of the vertical half-line above Φ_a's critical value: it leaves 0 at 45° and 225° and
  meets ±1/2 from above/below (|Im| ≤ 0.12).  In the repelling petals R_{2/3} sits at Im ζ = 2.1598 (± 1e-4, periodic),
  above every component, so deep gate steps are α side.  With this partition all 17 components of known limb are
  right (labels, bulbs B_1..B_4, mediant components confirmed in M incl. 9/20's P = 29 NRP, the m = 3 limbs); on the
  census, labels all land in the bulb limb and invalid limbs (b not prime to m) carry 5.6e-6 of 5.0e-4.
- Per limb (one strip side, r = 2..4 censuses): bulb limb 2.706e-4 + 5.34e-5 + 2.92e-5 (r = 2, 3, 4) beyond S_0;
  mediant 1.070e-4; t = 1/3: 1.33e-5; 2/3: 1.76e-5; 1/4, 3/4: ≥ 2.7e-7, 1.2e-7 (only r ≤ 4 seen).
- Tails: a limb-F component has kneading 1^(q_F - 1) 0 τ and pairs with the bulb-limb component of the same tail τ
  (for one-transit families τ is the census word; tokens '0' and runs 1^(2k g + c), normalized to the limb index).
  The copy ratio ρ_F(τ) = C_F(τ)/C_bulb(τ) depends mostly on the number g of gate passages in τ: mediant g = 0, 1, 2:
  0.0121, 0.325, 0.070; t = 1/3: 0.00084, 0.018; t = 2/3: 0.00053, 0.071; t = 1/4, 3/4 (g = 0): 0.00015, 0.00006.  So
  L_F ≈ Σ_g ρ̄_F(g) S_bulb(g) with S_bulb(g) = 1.893e-3, 2.706e-4, 5.34e-5, ≥ 2.92e-5: mediant 1.147e-4, 2/3: 2.02e-5,
  1/3: 6.5e-6 (half its census mass has no bulb partner yet).  g = 0 ratios fall ~m^-6.5; g = 1 ones are 20-130×
  larger and dominate.
- Structure: a tail run spanning g gates is itself what a Farey limb's first run is, so the whole strip is a renewal
  over kneading runs 1^(a_1) 0 1^(a_2) 0 … with a_i = 2k g_i + c_i: the Farey limbs are "first run spans m ≥ 2 gates",
  not a separate sum.  Next: a transfer operator over runs (weights per run type (g, c), learned from the census and
  checked against the per-limb data above), which would sum all limbs and transit counts at once.

**Toward a transfer operator over runs (2026-10-09, afternoon).**
- Low rank at a gate: C over (source, branch) × (target, shift) for the two-transit census has singular values 1,
  6.4e-3, 2.2e-3, 9e-4, 8e-5, 2e-5, 1.6e-6 (mass-weighted rank-3 error 1.2e-6), but entrywise and held-out errors
  are O(1): low rank compresses the heavy entries, it does not generalize.  A plain matrix-rank automaton is not it.
- Smoothness in the target: for one source (the bulb, all 16,383 targets) C/C_t is a smooth function of the target
  position y = σ_t + j/2: nearest neighbours within 1e-3 agree to 0.1% (90th percentile), within 1e-2 to 0.5–1.7%.
- **Area from the center**: from the return map's normal form (A = R''/2, D = ∂R/∂σ, area_σ = A_card/|AD|² for a
  small component), an r-transit component over target t has C ≈ C_t / |Θ_r'(σ) Π_{i<r} H'(p_i)|² (H = horn map;
  for r = 2, C_t/|H'(H'+1)|²).  Checked: r = 2, 479 components, median log10 ratio 0.000, IQR ±0.015 (targets below
  1e-4: median 0.0002); r = 3, 3000 components, median 0.0000, IQR ±0.0011 (±0.25%); large targets (bulb, j2_0) are
  1.5× off (the Jacobian varies over them: exact areas there).  `lavaurs_theta` returns Π H'; `--locate` prints it.
- A universal kernel (U small: C_child = C_U C_t |Φ_a'(v)|²/(K |Q'(u)|⁴), Q(u) = Φ_a(v + u²) - ζ0, with the island
  offset Θ(σ_U) = σ_U - (n_U+1)/2) gets the shape of the shift dependence but is off by 2–20× per target: the island
  map σ ↦ F^{n_U} Ψ(ζ0 + σ) carries O(1) Koebe distortion over |u| ~ 1, so there is no universal kernel.
- Children scale with the parent: over 77 single-transit sources from C = 1.7e-3 to 1.2e-12, their children's mass
  over a fixed target set (bulb + top 10, shifts around -(n_U+1)) is W_U ∝ C_U^(1 - 0.018), W/C median 0.017 (IQR
  0.009–0.031), scattering with the local distortion rather than following a smooth function of σ_U.  So beyond the
  top ~100 sources a level's sum hardly depends on the small ones (~5e-10 per level); completeness over big
  components, chains and satellite sources (the Farey bulbs) is what matters.
- ζ formulation: Θ(σ) = E(ζ0 + σ) - 2ζ0 with E(ζ) = H(ζ) + ζ, and H(ζ+1) = H(ζ) + 1, so all two-transit children of a
  target (every shift j) are the weighted preimages of one value under the cylinder map exp(4πi E): a genuine 1D
  transfer operator, smooth in the target; its critical points (E' = 0, H' = -1) sit next to the single-transit
  centers, which is why sources enumerate it.  Half of the Θ-preimages (the wrong petal parity for odd shifts of the
  bulb, and the targets themselves, where H' = 0) are not centers.
- Fast census (`lavaurs_area --children`, `fast_census.py`): children by Newton on Θ_r from the source's local
  quadratic, verified as centers (`lavaurs_center`), weighted by the area formula, deduplicated mod 1; heavy children
  exact (area, cusp).  300 sources × 101 targets reproduce 2.7532e-4 of S_1 = 2.94116e-4 in 2.6 min (every found
  component's area exact to the formula's accuracy); the missing 1.8e-5 are children of large sources whose islands
  the local start does not reach (H232 1.37e-5, H9386 3.0e-6, N-series): more starts per source, or path tracking for
  the top sources.
- Fast census levels 3–6 (sources: two-transit components > 1e-10 and satellites, heavy sources with 26 Newton
  starts, small ones with 2 and the top 30 targets; heavy children exact with cusp; satellites kept as sources):
  ~1.5 h on 4 cores.  It recovers the exact three-transit set (1.42239e-4 of 1.42260e-4) and finds what the
  path-tracked censuses missed (e.g. the children of H232, 1.37e-5, absent from island3's sources).
- Tunings are everywhere past three transits: U*X for every primitive X of period r in M (periods 3–6: 1, 3, 11, 20
  primitives, `primitives.py`; period 4 has three: -1.9408 and -0.1565 ± 1.0322i) and V*airplane for two-transit
  sources (B_2*A = 2.115e-5 at r = 6).  Tuned with LAVAURS_TUNE over the bulb + 200 families and 40 two-transit
  sources; removed by center (`assemble.py`).
- **Strip sectors (one side, all censuses merged, tunings and satellites removed)**: r = 1..6:
  1.8935e-3, 2.9412e-4, 1.5130e-4, 1.0378e-4, 7.683e-5, ≥2.265e-5 (level 6 least complete).  Ratios 0.155, 0.51,
  0.69, 0.74: beyond two transits the decay is slow (Farey limbs at r ≈ m + 1); the tail past r = 5 is 1–2e-4 if
  it continues, the dominant uncertainty in S_tot ≈ 2.65–2.75e-3.
- Levels by limb (fast census, components > 1e-9, tunings by single-transit U and two-transit V removed): r = 3:
  bulb limb 6.05e-5, mediant 8.37e-5, m = 3 2.6e-6; r = 4: bulb limb 6.88e-5, m = 3 2.8e-5, mediant 4.0e-6;
  r = 5: bulb limb 2.41e-5, mediant 1.60e-5, m = 3 1.23e-5, m = 4 9.8e-6, unclassified 1.44e-5; r = 6: bulb limb
  7.7e-6, m = 5 4.9e-6, m = 4 1.2e-6, unclassified 6.7e-6.  Each Farey denominator enters near r = m + 1.
- **Open: satellite tunings.**  Tunings of the satellites (the limb's bulb B_1, the Farey bulbs, doublings) by M's
  primitives are cusp-primitive components of the census and are not yet removed except where the copy-map guess
  found them (bulb*A 1.65e-4, bulb*Q 1.16e-5, B_2*A 2.1e-5).  The bulb's copy is strongly distorted: guesses for most
  of its period 4–6 tunings fail, and a 3.51e-5 four-transit primitive at -0.80458 + 0.32210i (n = 7 = the bulb-tuning
  excursion 4·1 + 3, found from the predicted bulb*P4b location) is likely bulb*P4b.  Estimate: the bulb's tuned
  primitives total ~0.208 × (M's primitive-copy fraction 1-2e-3) ≈ few e-4, comparable to the sectors, so r ≥ 4 above
  are upper bounds contaminated at the 1e-5 level.  Neither the limb classifier (tuned orbits pass the critical point
  closely, where the fixed arc is not the partition), nor a quadratic normal-form straightening (bulb*A → c = -2.99,
  not -1.75), nor passage distances discriminate.  Next: copy-map interpolation from many anchors (the bulb's own p/q
  satellites, H, A, Q) to place every bulb*X, B_m*X, and remove them; or a renormalization (little Julia set) test.
- **Satellite tunings removed via anchors** (`anchors.py`, `lavaurs_area --satellites`, `lavaurs_multiplier_point`):
  the bulb's satellite tree is placed exactly (the p/q satellite of a component is rooted at its multiplier map's
  point e^{2πip/q}, r q transits, excursion q n + q - 1; the M side by the same recursion from the cardioid), giving
  55 anchors (c, σ) for the copy map; a local affine fit in (c, c̄) predicts bulb*X to 0.002 (airplane), 0.004 (the
  period-4 primitive -0.1565 + 1.0322i: the 3.51e-5 four-transit component is bulb*P4), 0.03 (-1.9408, the stretched
  real direction).  All 35 primitives of periods 3–6 tuned (32 distinct components, 3 pairs landing on one σ):
  total 2.6618e-4 of bulb tunings (airplane 1.646e-4, period 4: 3.5e-5 + 1.4e-5 + 1.2e-5, period 5 ~3e-5, ...).
- Strip sectors after removing them: r = 1..6: 1.8935e-3, 2.9412e-4, 1.5130e-4, 5.461e-5, 4.675e-5, ≥1.331e-5.
- By limb (r ≤ 6): bulb limb 1.8935e-3, 2.711e-4, 6.05e-5, 1.97e-5, 7.6e-6, 4.2e-6 (Σ 2.2566e-3); mediant 2.30e-5,
  8.37e-5, 4.0e-6, 1.60e-5, 1.3e-6 (Σ 1.28e-4: transits in pairs, its own blocks decaying ~0.15–0.2 like the bulb
  limb's); t = 1/3, 2/3: 2.2e-5 each; 1/4, 3/4: ≥ 5.4e-6, 6.0e-6; fifths ≥ 5.1e-6 in all (entry level only).
- **Each limb's NRP layer is ∝ its bulb's area**: L_F/area(B_F) = 0.0108 (bulb limb), 0.0088 (mediant), 0.0072
  (m = 3), ≥ 0.006 (m = 4), slowly varying.  So the Farey part ≈ κ̄ Σ_F area(B_F) ≈ 0.007 · 0.24 (ζ(3)/ζ(4) - 1) ≈
  1.9e-4 (an m^-4.5 fit of the per-limb ratios gives 2.0e-4): S_strip ≈ 2.26e-3 + 2.0e-4 ≈ 2.46e-3 ± 5e-5.
- Picture: the NRP layer of a limb p/q is κ(p/q) × area(bulb p/q) with κ a slowly varying function of the
  continued fraction — the same fading-memory structure as the bulb areas' F(p/q).  A κ automaton over CF digits,
  computed per limb from Lavaurs-type models at each root, would sum all limbs, not only those near 1/2.
- Not a tuning test: |R_U^p(1/2) - 1/2| = 0 at σ_Y (U's return map iterated) only checks Y's excursion count
  n ≡ -1 (mod p) in the given representative, since F commutes with g_σ and R_U^p = F^{p n_U + p - 1} g^p F is the
  return map of every such component (necessary for tuning, not sufficient).  Dropping satellite sources does not
  avoid tunings either: Newton from primitive sources (H232's descendants) finds bulb*X too.  So deep levels need
  tunings placed explicitly (anchors) for every primitive X of each period — M has ~54, ~100, ~250 primitives of
  periods 7, 8, 9 — or the per-limb structure below instead of raw levels.
- Best current estimate (one side, k → ∞): bulb limb 2.2566e-3 (+ tail 4e-6), mediant 1.280e-4 (+4e-6), m = 3
  4.40e-5 (+5e-6), m = 4 1.10e-5 (+2e-6), m = 5 6.5e-6 (entry level × 1.3), m ≥ 6 1.1e-5 (κ = 0.005 × their bulb
  areas): **S_strip ≈ 2.4725e-3, uncertain at ~1e-5** (Farey extrapolation, per-limb tails).  At K = 64 that is
  ~3e-11 in μ(M) (before R); 1e-12 needs S_strip to ~3e-7.

**Toward a κ(p/q) automaton: per-limb NRP layers and a Lavaurs model at every root (2026-10-09, evening).**
- Direct limb censuses in M (`limb_families wake J p q` + bulb_batch: 485k NRPs for the 29 limbs p/q ≤ 1/2 with q ≤
  13 at extra period J ≤ 14 in 70 s): κ(p/q) = Σ NRP area / bulb area is 0.0023 (1/2), 0.0039 (1/3), 0.0058, 0.0077,
  0.0095, 0.0110, 0.0122, 0.0133, 0.0141, 0.0148, 0.0153, 0.0156 (1/4 .. 1/13), and along [2, k] (k/(2k+1)) 0.0045,
  0.0058, 0.0070, 0.0079, 0.0079 → 0.0108 (the σ-strip limit); [a, 2], [a, 3] smooth in a.  So κ depends on both ends of
  the CF, smoothly.  But each pass through a limb's own gate costs ~q in period (1/13 jumps 0.0142 → 0.0156 at
  j = q + 1, the second-pass family): direct enumeration (2^J) cannot reach J ≫ q, so large-q limbs need the pass
  (transit) structure at their own root.
- General-root Lavaurs model (`general_fatou.py`, `general_lavaurs.py`, prototype): at the p/q root f(w) = λw + w²,
  λ = e^{2πip/q}; the formal Fatou series Φ(f(w)) = Φ(w) + 1/q with the resonant log term β log(w^q)/q (per-petal
  branch), solved as a linear system (only a_0 is a gauge; the resonant a_m are fixed one order up); reproduces q = 2
  exactly (1/4, 1/4, 11/8) and gives for 1/3: a_{-3} = -0.071429 + 0.013746i, a_{-1} = 0.047619 - 0.137464i, β =
  1.773243 + 0.031420i (double precision is ill-conditioned past N ~ 12; C++ needs acb).  Φ_a by iteration into the q
  attracting petals, Ψ_k for the q repelling ones; single-transit centers f^n(Ψ_s(ζ0 + σ)) = -λ/2.  q = 2 reproduces
  S1 and the bulb (shifted by the branch constant (11/8)π); q = 3 gives a bulb-like component (area_σ ≈ 0.125) and
  NRPs 5.33e-4, 5.24e-4, 1.76e-5, 1.63e-5, … — top/second ratio 30 vs M's 29 in limb 10/31.  Next: calibrate σ
  against M near 1/3 (centers' 1/ε + q²k phases at k = 10..18), then port the pipeline.
