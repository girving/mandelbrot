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
