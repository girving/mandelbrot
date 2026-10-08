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
