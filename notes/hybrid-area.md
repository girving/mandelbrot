# Hybrid estimates of the Mandelbrot set's area

Notes from September 2026, branch `hybrid-area`.  Everything below is reproducible with the tools in this
repo plus the coefficient file `f-k27.npy` (2^27 Böttcher coefficients, 2 GiB) from
`https://storage.googleapis.com/mandelbrot/numpy/`.

## Summary

- **Estimate:** μ ≈ 1.507 ± 0.003, from the empirical law U_k ≈ μ + a/k for the Böttcher upper bounds
  (k = log₂ of the number of terms).  It is consistent with pixel counting (1.5066) and with an
  independent lower bound from hyperbolic components.
- **Lower bound:** the areas of all hyperbolic components of period ≤ 16 sum to **1.4987444**, and every
  component lies inside M.  The component areas are computed to ~1e-23 in `Expansion<2>`, but the sum is
  not interval-certified.
- **Structure:** Böttcher octave energy decays locally like a power of log n at parabolic points
  (j^-4 at satellite roots, about j^-6 at primitive cusps, with j = log₂ n), and like a power of n at
  Misiurewicz points.  Over the whole circle it decays like j^-2, which is a 1/log N tail.
- **Rigorous piece (elementary proof, §4):** the tail Σ_{n>N} n|b_n|² is at least c (log N)^-5, by
  shadowing real orbits near the cusp 1/4.  So the Böttcher bounds can never converge like a power of N,
  and β_M(2) = 1 for the integral means spectrum of the exterior map.

## Notation

ψ(w) = w + Σ b_n w^-n maps the outside of the unit disk to the outside of M.  The area is
μ = π(1 - Σ n|b_n|²), and U_k = π(1 - Σ_{n<2^k} n|b_n|²) is the upper bound after 2^k terms.
The tail is T(N) = Σ_{n>N} n|b_n|², so U_k - μ = π T(2^k).  Octave j is S_j = Σ_{2^j ≤ n < 2^(j+1)} n|b_n|².

## 1. Hyperbolic component areas (`hyperbolic`)

For each period p the tool:

- finds all 2^(p-1) roots of f_c^p(0) by Newton from rings of 4d points just outside |c| ≤ 2, keeps those
  of exact period p, and asserts the count against OEIS A000740;
- traces each component's boundary by continuing (z, c) along f^p(z) = z, (f^p)'(z) = λ from λ = 0 to
  |λ| = 1 in double, polishing in `Expansion<2>`;
- integrates Green's theorem, area = ½∫Re(c̄ λ c'(λ)) dθ, with the trapezoid rule.

| p | components | a_p | Σ_{q≤p} a_q | p³ a_p |
|---|---|---|---|---|
| 1 | 1 | 1.1780972 | 1.1780972 | 1.18 |
| 2 | 1 | 0.1963495 | 1.3744468 | 1.57 |
| 4 | 6 | 0.0232751 | 1.4542643 | 1.49 |
| 8 | 120 | 0.0044680 | 1.4858401 | 2.29 |
| 12 | 2010 | 0.0021014 | 1.4949217 | 3.63 |
| 13 | 4095 | 0.0008511 | 1.4957728 | 1.87 |
| 14 | 8127 | 0.0011339 | 1.4969067 | 3.11 |
| 15 | 16365 | 0.0009449 | 1.4978516 | 3.19 |
| 16 | 32640 | 0.0008928 | 1.4987444 | 3.66 |

Checks:
- The cardioid and the period-2 disk match 3π/8 and π/16 to within 2e-32.
- Doubling the number of boundary points changes each period's area by at most 7e-30.
- `Expansion<2>` Newton residuals are ≤ 8e-24 at p = 16.

Periods 1–16 take about 2 minutes on 18 cores.  Prime periods follow the cardioid-bulb law
a_q ≈ π(φ(q) - μ(q))/(2q⁴), where μ(q) is the Möbius function; composite periods add bulbs on bulbs.
Fitting a_p ~ p^-κ on p = 9..16 gives κ ≈ 2.2 and an extrapolated total of about 1.509.  Fixing κ = 3
gives about 1.504.

## 2. Where the Böttcher energy lives (`octaves`, `octscan`, `census`)

Octave j is Fourier-transformed at full resolution: n|b_n|² with n ∈ [2^j, 2^(j+1)) becomes an energy
density over external angle θ at resolution 2^-(j+1).

**Local laws.**  Energy within a few bins of a point, fitted over j = 14..26:

| point class | predicted | measured γ (S ~ j^-γ) | measured α (S ~ 2^(-αj)) |
|---|---|---|---|
| satellite roots (−3/4, −5/4, 1/5, 2/5 bulbs, …) | γ = 4 | 3.5–4.2, stable | drifts 0.30 → 0.23 |
| primitive cusps (−1.75, period-4 roots) | γ = 6 | 5.7–7.5 | drifts |
| main cusp 1/4 | γ = 6, contaminated by nearby 1/q bulb roots (γ = 4) | 4.8–5.0 | drifts |
| Misiurewicz (−2, i, angle 1/4) | power of n | γ drifts 22 → 30 | α ≈ 1.9–2.1, stable |
| whole circle | γ = 2 (§3) | 1.82 → 1.86 | 0.155 → 0.116 |

The predictions come from shell geometry.  Near a point c₀ where escape takes about τ(ε) iterations at
distance ε, the set {g_M < δ} reaches distance ε(δ) with τ(ε) ≈ log₂(1/δ).  The local shell area is the
area of the exterior of M within that distance:

- **Satellite root** (e.g. −3/4): escape takes about π/ε iterations (the "π in the Mandelbrot set"
  phenomenon).  The exterior near the root is the thin gap between two tangent circles, with area ~ε³.
  The shell area is ~L^-3 with L = log(1/δ), so T ~ L^-3 and S_j ~ j^-4.
- **Primitive cusp** (e.g. 1/4): escape takes about π/√ε iterations.  The exterior is a horn
  |Im c| ≲ 2(Re c - 1/4)^{3/2}, with area ~ε^{5/2}.  This gives ~L^-5 and S_j ~ j^-6.
- **Misiurewicz points:** ψ is Hölder there, so the local energy decays like a power of n.

**Stationarity.**  At a fixed angular resolution, the distribution of energy over θ converges as j grows.
Each hyperbolic component's gap set (the angles inside its wake but outside its sub-wakes of period ≤ 20)
keeps a nearly constant share of S_j for j ≥ 20.  Energy concentrates heavily near roots: at j = 26, 82%
of it lies within ±4 bins of the roots of periods ≤ 20, which cover 13.5% of the angles.

Superseded along the way:
- a "moving front" of primitive roots, which was an artifact of windows saturating the circle;
- an apparent match π·(peak energy near period-p primitives) ≈ a_p, which depended on the window width;
- multiplicative tuning scaling S_j ~ S_{j/2}, which is not supported by the data.

## 3. The 1/log N law and the estimate

**Heuristic.**  Every point of a hyperbolic component's boundary ∂W is neutral: the cycle's multiplier has
modulus 1.  A parameter just outside ∂W, where the multiplier has modulus 1+η, lingers about C/η cycle
steps before escaping.  So g_M ≈ 2^(-p C/η) there, and the shell {g_M < δ} has normal width
ε ≈ p C ln 2 / (L |λ'(c)|).  Integrating along ∂W,

  A(δ) ≈ (C ln 2 / L) · Q,  with Q = Σ_W p_W ∫|c'_W(e^{iθ})|² dθ,

so T(N) ≈ K / ln N and U_k ≈ μ + a/k, where a = πK/ln 2.

**Test on the actual bounds.**  With μ = 1.5066, the quantity (U_k - μ)·k ln 2/π should be constant:

| k | 8 | 10 | 12 | 14 | 16 | 18 | 20 | 22 | 24 | 26 |
|---|---|---|---|---|---|---|---|---|---|---|
| K | 0.835 | 0.859 | 0.869 | 0.864 | 0.871 | 0.872 | 0.868 | 0.870 | 0.869 | 0.867 |

Fitting U_k = μ + a/k with μ free:

| fit window | μ | a | rms residual | error in predicting U_27 |
|---|---|---|---|---|
| k = 12..20 | 1.50803 | 3.899 | 1.3e-3 | +8.7e-4 |
| k = 14..22 | 1.50909 | 3.882 | 9.8e-4 | +1.3e-3 |
| k = 16..24 | 1.50560 | 3.949 | 9.2e-4 | +2.6e-4 |
| k = 18..27 | 1.50469 | 3.970 | 6.1e-4 | (in sample) |
| k = 12..27 | 1.50723 | 3.911 | 1.0e-3 | (in sample) |

For comparison, fixed power laws in n predicted U_27 with errors of 5e-3 to 2.6e-2, and extrapolated to
1.58–1.66.  Adding a b/k² term makes μ unstable (1.504–1.524), so only the leading 1/k term is
determined.  Hence **μ ≈ 1.507 ± 0.003**.

Structural consistency: the weights Q_p = p Σ_{W of period p} ∫|c'_W|² dθ from `hyperbolic` are

| p | 1 | 2 | 3 | 4 | 6 | 8 | 10 | 12 | 14 | 16 |
|---|---|---|---|---|---|---|---|---|---|---|
| Q_p | 3.142 | 0.785 | 0.340 | 0.187 | 0.108 | 0.072 | 0.053 | 0.051 | 0.032 | 0.029 |

- The sum through p = 16 is 5.14, and the tail falls off roughly like p^-1.3, so Q ≈ 5–7.
- K = 0.867 then implies C = πK/(ln 2 · Q) ≈ 0.6–0.8.  That is O(1), as the heuristic needs.  It is not a
  sharp test, because C isn't computed independently.
- The heuristic is weakest at boundary points with irrational internal angle, and where bulbs crowd the
  boundary.

## 4. A rigorous lower bound on the tail

**Theorem.**  There are explicit constants c, N₀ > 0 such that for all N ≥ N₀,

  T(N) = Σ_{n>N} n|b_n|² ≥ c (log N)^-5.

Hence U_N − μ ≥ πc (log N)^-5.  In particular T(N) is not O(N^-a) for any a > 0, which means
β_M(2) = 1 for the integral means spectrum of ψ.  No extrapolation of the Böttcher bounds as a power
of N can be correct.

The proof is elementary.  Its only inputs are the Böttcher map Φ : C∖M → {|w| > 1}, Grönwall's area
theorem, and estimates on real orbits near the cusp c = 1/4.  All three exist in or near `ray`.

**Step 1 (shell identity).**  Let g_M = log|Φ| be the Green's function of M, and let
A(δ) = area{c ∉ M : g_M(c) < δ} = area ψ({1 < |w| < e^δ}).  Applying the area theorem to ψ on
|w| > R (after rescaling, ψ_R(w) = ψ(Rw)/R) gives the missed area πR² − πΣ n|b_n|² R^{-2n}.
Subtracting the R = 1 case,

  A(δ) = π(e^{2δ} − 1) + π Σ_n n|b_n|² (1 − e^{−2nδ}).

**Step 2 (the shell area is controlled by the tail).**  Use 1 − e^{−2nδ} ≤ 2nδ for n ≤ M, and ≤ 1 for
n > M.  Also Σ_{n≤M} n·n|b_n|² ≤ M Σ n|b_n|² ≤ M, since the missed area is nonnegative.  So

  A(δ) ≤ π(e^{2δ} − 1) + 2πδM + πT(M).

Take M = ⌊δ^{-1/2}⌋.  This gives πT(M) ≥ A(M^{-2}) − O(M^{-1}).

**Step 3 (a horn of slowly escaping parameters at the cusp).**  Fix 0 < ε ≤ 1/16 and the real
parameter c = 1/4 + ε.  Let z_0 = 0 and z_{n+1} = z_n² + c.  Set u_n = z_n − 1/2 and
φ(u) = u² + ε, so that u_{n+1} = u_n + φ(u_n).  The sequence u_n increases from −1/2 to +∞.

- *(3a) Slow passage.*  The map u ↦ u + u² + ε is increasing for u > −1/2, and sends −√ε to −√ε + 2ε.
  So the first u_n ≥ −√ε is at most −√ε + 2ε.  While |u_n| ≤ √ε, each step is at most 2ε.  Hence
  z_n ≤ 1/2 + √ε ≤ 1 for all n ≤ n₁ := ⌈1/√ε⌉ − 2.
- *(3b) A telescoping identity.*  φ(u_{n+1}) = φ(u_n)·(1 + 2u_n + φ(u_n)).  So
  2z_n = 1 + 2u_n ≤ φ(u_{n+1})/φ(u_n).
- *(3c) The key sum.*  Let n₂ be the first n with u_n ≥ 1.4, i.e. z_n ≥ 1.9.  Then
  Σ_{k=1}^{n₂} 1/φ(u_k) ≤ 10 ε^{-3/2}.

  *Proof, by telescoping.*
  - *Terms with u_k < −√ε.*  Put v = −u.  Then v_{k+1} ≤ v_k(1 − v_k), so
    1/v_{k+1}³ − 1/v_k³ ≥ 3/v_k².  Since 1/φ ≤ 1/v², these terms sum to at most ε^{-3/2}/3.
  - *Terms with |u_k| ≤ √ε.*  Each step is at least ε, so there are at most 2/√ε + 1 of them, each at
    most 1/ε.  They sum to at most 2.25 ε^{-3/2}.
  - *Terms with √ε < u_k < 1.4.*  Write u_{k+1} = u_k(1 + ρ) with u_k ≤ ρ ≤ 1.65.  Then
    1/u_k³ − 1/u_{k+1}³ ≥ (3/2.65³)/u_k², so these terms sum to at most 6.2 ε^{-3/2}.
  - *The last term,* at k = n₂, is at most 1/1.96.

  Numerically the three parts tend to 0.14, 1.28 and 0.14 times ε^{-3/2}.  The whole sum tends to
  (π/2) ε^{-3/2}.
- *(3d) Shadowing.*  Let c' = c + h with |h| ≤ ε^{3/2}/1000, let z'_n be its orbit, and let
  e_n = z'_n − z_n.  Put τ_n = |h| Σ_{k≤n} 1/φ(u_k), so that τ_n ≤ 1/100 for n ≤ n₂ by (3c).

  *Claim:* |e_n| ≤ τ_n φ(u_n) for 1 ≤ n ≤ n₂.

  *Proof, by induction.*  e_1 = h, so the case n = 1 holds.  For the step, e_{n+1} = (2z_n + e_n)e_n + h.
  Since z_n ≥ 0, (3b) gives |2z_n + e_n| ≤ 1 + 2u_n + φ(u_n) = φ(u_{n+1})/φ(u_n).  Hence
  |e_{n+1}| ≤ φ(u_{n+1})(τ_n + |h|/φ(u_{n+1})) = τ_{n+1} φ(u_{n+1}).  Numerically the bound is sharp.

  Consequences:
  - At n₂ we have z_{n₂} ≥ 1.9 and φ(u_{n₂}) ≤ 12, so |z'_{n₂}| ≥ 1.78.  Then
    |z'_{n₂+1}| ≥ 1.78² − |c'| > 2, so c' escapes and c' ∉ M.
  - For n ≤ n₁, |z'_n| ≤ 1 + 1/100·φ(u_n) ≤ 1.1.  Using g_c(z) ≤ 1.2 for |z| ≤ 2 and |c| ≤ 1/2,
    and g_M(c') = 2^{1−n} g_{c'}(z'_n), we get g_M(c') ≤ 2.4 · 2^{−n₁}.
- *(3e) The horn.*  For ε ≤ 1/32, the rectangle H_ε = {1/4 + x + iy : ε ≤ x ≤ 2ε, |y| ≤ ε^{3/2}/1000} applies
  (3d) with base point 1/4 + x, since ε^{3/2} ≤ x^{3/2}.  It lies outside M, satisfies
  g_M ≤ δ(ε) := 2.4 · 2^{2−1/√(2ε)}, and has area ε^{5/2}/500.

**Step 4 (conclusion).**  Since log(1/δ(ε)) ≍ (ln 2)/√(2ε), Step 3 gives A(δ) ≥ c₁ (log 1/δ)^{-5}.
With δ = M^{-2}, Step 2 then gives T(M) ≥ c (log M)^{-5}. ∎

The constants are far from sharp: 10 could be π/2, and 1/1000 could be about 1/(50π).  The exponent 5
comes from the horn's area, ε^{5/2}, where the width ε^{3/2} is exactly the width at which shadowing
still works.  The heuristic in §3 predicts the true rate is 1/log N.

**Formal version (Lean, `ray` branch `cusp-tail`).**  The argument is formalized in `ray`, where
Grönwall's theorem for M already existed as `multibrot_volume_sum`.

- `Cusp.lean` proves the passage lemma, the key sum via the potential H(u) = u/(2s²φ) + arctan(u/s)/(2s³)
  (whose derivative is 1/φ²), shadowing, and escape.
- `CuspArea.lean` shows these parameters lie outside M with potential ≥ exp(−2/2^m).
- `Shell.lean` applies Grönwall at radius R to get the shell area.
- `Tail.lean` assembles the pieces and proves `areaTail_pow_two_ge`:

      ∀ᶠ k in atTop, 1 / (2·10⁸ (2k+1)⁶) ≤ ∑_{n>2^k} n |b_n|²

  It also proves `areaTail_pow_two_gt_rpow`: for every a > 0, eventually the tail is greater than
  (2^k)^-a.

The formal version uses a disk of parameters instead of the horn, which costs one power of log:
the exponent is 6 rather than 5.  It depends only on the axioms `propext`, `Classical.choice` and
`Quot.sound`.

**Stronger versions, not yet proved.**
- At the satellite root −3/4, the exterior gap between the cardioid and the period-2 disk has area
  about s³ within distance s, and escape there takes about π/s iterations (Boll's π; proved on the
  vertical line by Klebanoff, Fractals 9 (2001) — citation not checked).  This would give T(N) ≳ (log N)^-3.
  Its exterior parameters are not real, though, so the shadowing proof above does not apply; it would
  need the parabolic normal form of f² made uniform over the gap.
- A proof of the 1/log N rate would need the neutral-boundary shell of §3 at every component, including
  boundary points with irrational internal angles.  That is well beyond the elementary methods here.

## 5. Digits: fattened sets from escape times, bias from the 1/log law

Grönwall's theorem at radius e^δ gives the area of the fattened set {c : g_M(c) < δ} exactly:

  F(δ) = π e^{2δ} − π Σ n |b_n|² e^{−2nδ}.

The weights e^{−2nδ} decay, so the 2^27 coefficients determine F(2^-k) to about 10^-8 for k ≤ 19
(`fattened`).

The same sets are cheap to measure by Monte Carlo far past what the series can reach.  If the orbit of
c has not escaped after m iterations, then g_M(c) ≤ 2^-m · O(1).  So deciding whether g < 2^-k takes
about k iterations, plus cardioid, disk and cycle tests for the interior of M.  The tool is
`escape_area`: a jittered grid with 8 independent points per cell, which serve as replicas.

**Validation.**  With 3 × 2·10⁹ samples, the Monte Carlo F(2^-k) for k = 8, 12, 16, 18, 19 matches the
exact Grönwall values to within 1σ, which is about 10^-6.  This checks both the classifier and the
coefficients.

**Results.**  The run goes to k = 2^20, which corresponds to 2^(10^6) Böttcher terms.  The effective
1/k coefficient drifts slowly, from 3.8 at k = 2^10 to 2.8 at k = 2^19; at k ≈ 20 the Böttcher fits give
πK/ln 2 ≈ 3.9.

Two tail models fitted on k ≥ 2^14 agree to 2·10^-8: μ + a/k^{1.045}, and μ + a/(k·√log k).  With 9 seeds
of 2·10⁹ samples each:

  **μ = 1.5065931 ± 0.0000006 (stat) ± 0.0000015 (rounding) ± ~0.0000005 (tail model)**

**Checks.**
- *Rounding.*  On 4·10⁸ samples, double and double-double classifications flip in both directions at
  nearly equal rates.  The net bias is ≲ 1.5·10^-6 up to k = 2^17.
- *Agreement.*  The result matches the pixel-counting estimate 1.5065918849 to about 10^-6.
- *Precision gained.*  Böttcher extrapolation alone gave μ ≈ 1.507 ± 0.003, so the escape-time
  measurement adds about three digits.

### 5.1 Toward two more digits (September 2026)

Target: error ~3·10^-11, two digits beyond the published 1.5065918849.

- *Certified adaptive tree* (`escape_tree`).  Cells whose center carries a Koebe distance certificate are
  decided exactly; the rest split, down to leaves sampled with m points.  Leaf samples are jittered in
  2 × 2 strata (1.64× in error² × time over iid points).  Error ∝ cost^-0.69.  Leaves are collected and
  sampled in batches.  The whole pipeline (`tree`, on `engine`) runs on CPU threads or the GPU (`--cuda`).
- *Failed ideas*: control variates from known components, smoothing, roulette, multilevel over thresholds,
  linearization jumps and BLA (orbits leave the linear regime at once), and sampling only the exterior
  shell {2^-K ≤ g < 2^-19} against the exact Grönwall F(2^-19).  The last one fails because at leaf scale
  the level curve g = 2^-19 is as wiggly as ∂M: var(shell) ≈ var(A(19)) + var(A(K)).
- *Float is biased.*  `escape_tree --prec compare` classifies every leaf sample in float and double.  At
  k = 2^20 float overcounts by +1.2·10^-5 ± 8·10^-8, growing with k.  The flips are float orbits wrongly
  declared non-escaping: Brent cycles 65%, Newton certificates 27%, max_iter 9%.  These are finite-state
  and tolerance artifacts, so float is out; H200 double is only 2× slower anyway.
- *Double rounding.*  For slow escapers, double and double-double escape steps differ by about the step
  count itself (sd ~7·10^4 at 2^16 steps).  Rounding re-randomizes the orbit, so it acts like a
  ~10^-16 jitter of c.  Averaging the indicator over a jitter kernel preserves its integral, so the bias is
  second order (c-dependence of the kernel) plus non-shadowing artifacts like the float ones.  This is an
  argument, not a proof.  A paired double vs double-double tree run can check it only to ~10^-9.
- *H200 pipeline (2026-09-28).*  The whole tree runs on the GPU (engine.h: persistent warp-synchronous
  kernels with scrambled claims; Newton certificates deferred and settled in lockstep rounds; 32-bit orbit
  state and register budgets).  Results are bit-identical to CPU runs.  Throughput: bare z² + c loop
  2.4·10^12 it/s, Orbit::run 1.2·10^12, leaf sampling about 6·10^11 overall (rounds at 1.0·10^12), centers
  1.4·10^11.  Error vs H200 time at base 1000, 16 samples/leaf, strata 2:

  | depth | error at k = 2^20 | seconds |
  |---|---|---|
  | 6 | 1.33·10^-7 | 4.1 |
  | 8 | 2.02·10^-8 | 37.7 |
  | 9 | 7.88·10^-9 | 131 |
  | 10 | 3.09·10^-9 | 478 |

  Error ∝ cost^-0.72 between depths 9 and 10 (0.85 at 6→8).  At ~0.7, 3·10^-11 needs roughly 3.5–4
  H200-days.  Depth 10 also gives A(2^18) − A(2^20) = 8.5761·10^-6 ± 1.6·10^-9.
- *Budget (earlier, CPU).*  The tree reaches 3.4·10^-7 in 44 s on an M5 Pro (3.2·10^9 iterations/s).  3·10^-11 then
  needs ~335 laptop-days.  An H200 does ~34 TFLOPS in double, about 3·10^12 iterations/s at peak; at
  30–50% of peak that is 300–500× the laptop, so 0.7–1.1 GPU-days.  Tree traversal (escape_de at cell
  centers, 15% of CPU time) would then dominate and needs to move to the GPU too.  The tail extrapolation
  in k (the tail beyond 2^20 is ~2.7·10^-6, needed to ~10^-5 relative) is the other open risk.

### 5.2 Science risks: the tail in k and rounding (2026-09-28)

*The tail in k.*  A(k) − μ = π Σ n|b_n|² (1 − e^{−2n·2^−k}), a 1/log N tail at N ≈ 2^k.  Measured differences
D = A(k) − A(k') at half-octaves from k = 2^10 to 2^24 (depth 8–9) and octaves to 2^28 (depth 6–7) rule out every
pure power series in 1/k with k ≥ 2^16 (χ²/dof ≥ 26), and over k ≥ 2^16 plausible models moved μ by 1e-7 to
3e-7.  Restricting fits to k ≥ 2^20, with data to 2^26, the well-fitting models agree:

| model (k ≥ 2^20) | χ²/dof | μ from A(2^24) |
|---|---|---|
| a/k + b/k² + c/k³ | 0.64 | 1.5065918672 |
| a/k + b log k/k² + c/k² | 0.75 | 1.5065918684 |
| a/(k log k) + b/k | 1.09 | 1.5065918697 |
| a/k | 15.5 | 1.5065918580 |

(A(2^24) itself has statistical error 2e-8.)  The effective coefficient K(j) = (A − μ) j ln2/π, with the Böttcher
bounds U_j at j ≤ 27 and escape areas at j* = k − 1 − γ/ln2, is flat at 0.865 for j = 10–27 and falls slowly with
ln j for j = 10³–4·10⁶ (0.81 → 0.59 at μ = 1.50659188); the large-j values depend on μ (ΔK = Δμ j ln2/π), and a
single model fit on k ≥ 2^16 does not predict the Böttcher values (off by 1.7–2.5×).  So the asymptotic form is
not known, and the robust route is to measure D directly to k ≈ 2^30, where the remaining tail is ~2.5e-9.

*Rounding.*  Paired runs (the same samples in double and in double rounded to b bits by Veltkamp splitting, same
algorithm and tolerances), depth 8, max_iter 2^20: at k = 2^20 the rounded − double difference is +2.8e-7 ± 4.5e-9
at 30 bits, and −0.9e-9, +1.6e-9, +3.0e-9 (± 4e-9) at 36, 42, 48 bits.  The 30-bit bias grows with k, as expected if
slowly escaping orbits near parabolic points stall once their increments (~1/N² at cusps) drop below the rounding
step; scaling N ~ √(C/ε) from 30 bits suggests double could stall near k ~ 2^27–2^28.  At max_iter 2^26 (depth 6),
42 and 48 bits agree with double to ≤ 1.3σ through k = 2^25 and are 1.7σ low at 2^26.  Float's earlier bias
(+1.2e-5) came from its looser cycle and Newton tolerances, not from rounding as such.  Double-double orbits
(--prec comparedd) give a direct reference at large k.

*Brent false positives.*  Orbits past 2^24 steps (max_iter 2^28) were mostly "certified" by Brent right after
the 2^24 checkpoint, but Newton there finds repelling cycles (|λ| from 1.000007 to 1.4): slow exterior orbits
returning within 1e-13.  Newton now has to confirm Brent's cycle, and an unconfirmed exact (bitwise) return
stops the orbit as capped, since the computed orbit is then periodic and would reach max_iter anyway.  At depth 6
and max_iter 2^28 the estimates are bit-identical to before, so the false positives did not matter there.  In
orbit_census at 2^28 the exact-return rule cuts work 72×.  (A bug in the first version, squaring a −1 default
tolerance to +1, certified exterior orbits and moved A(2^20) at depth 8 by 8e-5; tree_newton_matches_reference now
checks the Newton schedule against escape().)

*Cost of large k.*  Depth 8: max_iter 2^20 → 2^24 costs 37.7 → 84.4 s.  Depth 6: 2^28 costs 19× 2^20: orbits that
neither escape nor certify run to max_iter.

### 5.3 Local approximations of long orbits: measured headroom (2026-09-28)

Could Koenigs linearization near repelling cycles, Fatou coordinates through parabolic gates, renormalization
near baby copies, per-component tails, or two-phase importance sampling cut the work?  `orbit_regimes` samples
leaf-like points (as `orbit_census`), replays orbits longer than `--long` in windows, and at each window start
finds the smallest period q ≤ 1024 returning within 1e-2, Newton-solves for that q-cycle, and classifies the
window.  Runs of consecutive windows near one cycle (identified up to phase) bound what a jump through the
cycle's local dynamics could skip: everything but a run's first and last windows.  On the GPU node's 22 CPUs:

| max_iter | samples | long-orbit share of work | skippable (window 4096 / 512) | speedup bound |
|---|---|---|---|---|
| 2^20 | 160k | 69% | 13.7% / 17.3% | 1.16–1.21× |
| 2^24 | 4M | 51% | 11.3% / 11.3% | 1.13× |
| 2^28 | 1M | 39% | 9.8% / 10.0% | 1.11× |

Long-orbit steps are 71–73% near repelling cycles, but that is mostly the triviality that Julia-set points are
near high-period repelling cycles; they leave them quickly.  Multipliers are spread over | |λ| − 1 | from 1e-7
to 1, with no dominant near-parabolic population (strictly parabolic, |λ| ≈ 1 at rational rotation: ≤ 1%).
Renormalizable windows (approaches to 0 only at multiples of p ≥ 2): 3.7–3.9%.  **Verdict: Koenigs jumps,
Fatou gates, and renormalization together are bounded by ~1.1–1.2×**, before counting their per-jump cost
and certification complexity.  Not worth building.

*Two-phase importance sampling.*  Flag a 16-sample box by 8 of its samples surviving past T1; measure how many
of the other 8 samples' survivors past T2 fall in flagged boxes; Neyman allocation then scales the tail
variance per sample by (√(φψ) + √((1−φ)(1−ψ)))².  Best case T1 = 2^14: 28% of boxes flagged hold 90–94% of the
deep survivors, variance × 0.52–0.59, stable from T2 = 2^16 to 2^22.  Flagging deeper is worse (T1 = 2^16:
× 0.82–0.86), since deep survivors are not predicted by other deep survivors (1.01–1.04 per box, as random).
So ≤ 2× on the tail's variance, and less on the total; the tree's leaves already localize the boundary.

*Two-phase allocation across leaves (`--leaf-stats`).*  Class leaves by their first P samples, measure each
class's variance and cost on held-out samples 8..15, and allocate a second phase by cost-weighted Neyman
(n ≥ 1 per leaf, estimating from the second phase only, so unbiased).  Depth 8, A(2^20): cost at equal
variance 0.99 / 0.97 / 1.12 of uniform for P = 2 / 4 / 8.  97.7% of leaves have all 8 pilot samples outside,
but they still carry 8.5% of the variance (and all-inside leaves 5.4%) and 37% of the cost each, so they cannot
be starved; without the pilot's cost the gain would be 1.6×, which a free predictor (the leaf center's orbit)
might partly recover.

*Stratification is the same as depth.*  Depth 8 at A(2^20) (std err, total time, err²·t):

| run | finest cell | groups | std err | time | err²·t |
|---|---|---|---|---|---|
| d8 m16 s2 | h_9 | 4 | 2.02e-8 | 35.7 s | 1.46e-14 |
| d8 m32 s1 | h_8 | 32 | 1.83e-8 | 61.4 s | 2.06e-14 |
| d8 m32 s2 | h_9 | 8 | 1.43e-8 | 61.4 s | 1.26e-14 |
| d9 m16 s2 | h_10 | 4 | 7.88e-9 | 131 s | 8.1e-15 |
| d8 m32 s4 | h_10 | 2 | 1.11e-8 | 61.4 s | 7.6e-15 |
| d7 m128 s8 | h_10 | 2 | 1.11e-8 | 70.4 s | 8.7e-15 |
| d8 m64 s4 | h_10 | 4 | 7.88e-9 | 118 s | 7.3e-15 |
| d10 m16 s2 | h_11 | 4 | 3.09e-9 | 478 s | 4.6e-15 |
| d9 m32 s4 | h_11 | 2 | 4.37e-9 | 211 s | 4.0e-15 |
| d8 m128 s8 | h_11 | 2 | 4.37e-9 | 232 s | 4.4e-15 |

(h_d is the depth-d cell size.)  The variance is set by the finest stratum size and the samples per stratum,
not by how the cells got there: runs with the same stratum size agree to three digits (d7 m128 s8 = d8 m32 s4,
d8 m128 s8 = d9 m32 s4), and at fixed stratum size variance ∝ 1/groups (d9 m16 s2 = d8 m32 s4 / 2).  So
err²·t improves about 1.9× per halving of the stratum size however it is reached; strata vs depth is a ≤ 12%
cost trade, and one sample per stratum (G = 1, with a collapsed-strata variance estimate) would not change
err²·t.  Projection from d9 m32 s4 at 3.43× cost and 6.45× variance per level: statistical error 3e-11 at
about 5.4 more levels, ~1.5e5 H200-seconds (~43 GPU-hours), before the tail extrapolation.

*Per-component tails.*  The period of the last hugged cycle of escaping survivors is broad (9–256, mode
33–64) and the distribution does not change from length 2^12 to 2^22: no finite set of components carries the
tail, so a per-component analytic tail would need infinitely many components with self-similar weights —
which is the tail-fit problem again, not a shortcut.

*Log-periodic test of the tail.*  If tuned copies scale the tail as T(k) ⊃ |ρ_W|² T(k/p_W) (Green's function
g_M(c) ≈ C g_M(χ_W(c))^p near a period-p copy), the tail's Mellin transform satisfies
T̂(s)(1 − Σ_W |ρ_W|² p_W^s) = T̂_0(s), and complex roots would give log-periodic terms k^−σ cos(ω ln k).  Fitting
the 46 distinct differences of §5.2 with smooth bases plus (c cos ω ln k + d sin ω ln k)/k, ω ∈ [0.6, 12] (at
least one period over ln k ∈ [6.9, 17.3]):
- the best smooth model is a/k + b/(k ln k) + c/k², χ²/dof 1.31 from k ≥ 2^14 (35 dof) and 1.46 from 2^16;
- adding the oscillation gains Δχ² ≤ 4 (k ≥ 2^14) or ≤ 11 (k ≥ 2^16, ω ≈ 2.4–3.6) for 3 extra parameters over
  ~20 effective trial frequencies: not significant;
- 2σ amplitude bounds relative to a/k: < 2% at ω = 0.6 (degenerate with smooth terms), < 0.4% at 1.2, < 0.2% for
  ω ≥ 1.8, < 0.05% for ω ≥ 4.2.  At k = 2^30 (T ≈ 2.3e-9) that is ≤ 5e-12 for ω ≥ 1.8 but up to 5e-11 for the
  slowest periods, which only data at larger k (or the D(s) roots from the component areas) can separate from
  smooth drift.
Across well-fitting models and ranges (k ≥ 2^14 … 2^20), T(2^24) = 1.567–1.608e-7 (the ±2e-9 μ spread) and the
extrapolated T(2^30) = 2.27–2.50e-9: measuring to 2^30 leaves a model spread of ~2e-10, still ~7× the goal.

*The tail as a sum over roots (`rootsum`).*  Prior evidence (§2, `octscan`'s inverse-tuning test) is against
copies reproducing the whole tail at k/p, so instead of a renewal equation we tested an additive model
S_j ≈ Σ_W w_W Ψ_type(j/p_W): each root switches on at j ~ p and then decays with its local law.  For octaves
j = 12–26, `rootsum` sums the energy within ±4 bins of both rays of every root of period p ≤ j − 7 (low periods
first, so each bin counts once; the windows cover 6–7% of the angles):
- roots of period ≤ j − 7 hold ~70% of S_j on 7% of the angles: 35% at satellite roots (stable in j) and
  10% → 38% at primitive roots (growing with j); ~28% is left for higher periods and non-root angles;
- each period's energy E_p(j) is modulated with period p in j (the rays return under p doublings: j → j + p),
  so envelopes are averages over p consecutive octaves;
- satellite envelopes: Ē_p(j) = W_p j^−γ with γ = 3.93 (theory 4), rms 0.08 in log over p = 2–10; the weights
  w_p = W_p p^−γ give p³ w_p ≈ 2.1–2.9 for p = 4–10 (6.6, 3.8 at p = 2, 3; composites a little high), i.e.
  w_p ~ p^−3, and satellites alone would give S_j ~ j^−2, T ~ 1/k: the observed law;
- primitive envelopes: γ ≈ 6.5 (theory 6) but rms 0.3, and extrapolated weights p³ w_p grow (3.7 → 14 over
  p = 3–10).  Their energy sits near the frontier p ≈ j/1.4 … j where Ψ is not measured, and their share grows.
So satellites explain the 1/k law cleanly, but the primitive frontier, which the Böttcher data (j ≤ 26) cannot
follow to the escape-time range (j ~ 10³–10⁶), increasingly carries the tail and plausibly the drift of K_eff
there.  A principled tail form would need the frontier behaviour of primitive roots at large period.

*Tuned copies carry p× the tail per unit area (escape side, `--box`).*  Near a period-p copy,
g_M(c) ≈ C g_M(χ(c))^p, so the copy's fattened sets should be M's at threshold ~k/p, and its tail per unit
area p times M's (a/(k/p) = p·a/k), where plain area scaling predicts 1.  Depth 8, max_iter 2^22, octave
differences D(k) = A(k) − A(2k), ratio R = (D/A)_box / (D/A)_M with A the area at 2^22:

| k | M: k·D/A | period-3 copy, box [−1.7930, −1.7485] × [0, 0.024] | period-4 copy, box [−1.94290, −1.94045] × [0, 0.0013] |
|---|---|---|---|
| 2^8 | 1.315 | 2.878 | 4.129 |
| 2^12 | 1.201 | 2.853 | 3.988 |
| 2^16 | 1.060 | 2.857 | 3.996 ± 0.002 |
| 2^20 | 0.957 | 2.826 ± 0.003 | 3.938 ± 0.009 |
| 2^21 | 0.938 | 2.813 ± 0.005 | 3.920 ± 0.013 |

(Copy areas 4.999e-4 and 1.449e-6, doubled boxes.)  The enhancement is the period, to 5–6% for p = 3 and 1–2%
for p = 4, flat over k = 2^8 … 2^21, so the multiplicative scaling holds on the escape side; the Böttcher-side
rejection (§2) was pre-asymptotic (j/p ≤ 13).  With M's own drift, T(k/p)/T(k) predicts R = p·K(k/p)/K(k) ≈
3.11 and 4.19 at 2^20, so the effective tail weights are 0.91 and 0.94 of the area ratios (the straightening
distorts the boundary region differently from the bulk).  So the renewal equation
T(k) = T_0(k) + Σ_W r_W T((k − c_W)/p_W) over maximal copies is quantitatively right, and D(s) = Σ_W r_W p_W^s
governs the tail's form.

*D(s) from the tuning structure (`tuning`).*  `hyperbolic 16 1024 1 components.txt` dumps every component's
center and area; `tuning` classifies each Lavaurs root as non-renormalizable or tuned (`maximal_tuning` in
angles.h: both angles are concatenations of a smaller component's two words), matches roots to centers by
tracing the lower ray to near the root (a bijection over all 65242 roots, second-nearest center ≥ 2.9× farther;
0.6 s), and sums D(s) = Σ_{non-renormalizable W} (area(W)/area(cardioid)) p_W^s.  Through period 16:

| p | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| p³ a_nr_p | 1.57 | 1.53 | 0.79 | 1.64 | 0.28 | 1.73 | 0.89 | 1.19 | 0.55 | 1.84 | 0.63 | 1.87 | 0.68 | 0.87 | 0.96 |

(primes are all non-renormalizable; composites lose most of their area to copies).  D(1) = 0.676 and the real
root of D(s) = 1 is s* = 1.33 at P' = 16, still falling with P'; extending with p³ a_nr_p = 0.7–1.6 beyond 16
gives D(1) = 0.71–0.76 and s* = 1.18–1.25.  All complex roots have Re s > 2.9.  So D(1) < 1: the 1/k term
comes from the non-copy part T_0, amplified by 1/(1 − D(1)) ≈ 4, and renormalization adds a correction
k^−s* with non-integer s* ≈ 1.2 (complex roots only at Re s > 2.9, so no visible log-periodicity, as §5.3
found).  Independently, fitting the tail differences with a/k + b k^−s + c/k² prefers s = 1.20–1.25 from
k ≥ 2^16 (χ² 40.7 at s = 1.20 vs 45.1 for a/k + b/(k ln k) + c/k², 31 dof) and s ≈ 1.15 from k ≥ 2^14.  The
extrapolated T(2^30) is 2.36–2.40e-9 at s = 1.20–1.25 against 2.27e-9 for the log model.

A cheap shortcut to the tuning classification, closest returns of the center's critical orbit (renormalizable
with period d iff d is a closest return and all later ones are multiples of d), never misses a tuned component
but calls 13866 of 65242 non-renormalizable ones tuned (21%, mostly p ≥ 12): necessary, not sufficient.

*Tail weights, and D in closed form.*  More `--box` copies (depth 8, max_iter 2^22), as tail weight r_W over
area(W)/area(cardioid), using R/(p K(k/p)/K(k)) and the copy's area: primitive copies 0.875 (period 3),
0.928 (period 4), 0.893 and 0.890 (period-5 copies at −1.8608 and −1.6254); the period-2 satellite copy
(box [−1.56, −0.75] × [0, 0.32], R = 1.48 flat in k) 0.627.  (The period-5 copy at −1.9854 is lost: its
area, 5e-9, underflows escape_tree's %.10f output.)  Satellites hold 97–99% of the non-renormalizable area
in every period, so D is essentially the cardioid bulbs', and since a satellite of any other component
lies in that component's copy, the non-renormalizable satellites are exactly the cardioid's bulbs.  With
Σ_r area(r/q bulb) ≈ π(φ(q) − μ(q))/(2q⁴) (§1) and a1 = 3π/8, their Dirichlet series is
  D(s) ≈ (4w/3) (ζ(3 − s) − 1) / ζ(4 − s),
which gives D(0) = 0.2489, D(1) = 0.715, s* = 1.238 at w = 1 (enumeration through 16 plus extrapolation:
0.248, 0.72–0.75, 1.22–1.24).  s* depends on the satellite weight: 1.31 at w = 0.89, 1.50 at w = 0.63.  The
tail fits' s ≈ 1.20–1.25 matches w ≈ 1, not the single satellite calibration, so either other bulbs'
copies carry more tail than the period-2 copy, or the period-2 box mixes in non-scaling tail (the root
horn at −3/4, a tangency rather than the image of M's cusp); the agreement is not yet established.  D has
poles where ζ(4 − s) = 0, i.e. at Re s = 7/2 on the critical line, and at s = 2 from ζ(3 − s).

*Satellite tail weights.*  Box runs on satellite copies (depth 8, max_iter 2^22), as r_W / (area(W)/a1):
1/2 bulb 0.628 (box to the tip −1.5437; a taller box gives the same R = 1.480, and cutting the root region at
−0.80 lowers R to 1.41, so the root horn carries more tail per area, not less); 1/3 bulb 0.662 (box
[−0.235, 0.005] × [0.64952, 0.95629], from the cardioid's top to the limb's principal branch point, which is
the copy's tip); the 1/4 bulb box [0.25, 0.372] × [0.5, 0.598] cuts off part of the bulb (copy area only
0.72 of scaled) and is unreliable.  Both clean satellite copies have copy area 0.871 of the scaled bulb and
tail ratio R 0.72–0.76 of pure scaling, flat from k = 2^12; primitive copies (now five, including the
period-5 copy at −1.9854 at 5.0e-9 area) stay at 0.87–0.93.  With satellite weight 0.65, D predicts
s* ≈ 1.45, and the tail fits cannot tell: with a log term, a/k + b/(k ln k) + d k^−s + c/k² has χ² 40.4–41.1
(k ≥ 2^16) and 45.4–46.0 (k ≥ 2^14) for every s in 1.2–1.6, with T(2^30) = 2.32–2.38e-9 (k ≥ 2^16).  So the
s ≈ 1.2 agreement above was fortuitous; D gives the correction exponent (1.45 ± 0.05 if larger bulbs weigh
like these two) but the tail data do not yet constrain it.

*The non-copy tail T_0.*  Subtracting the copies' share with the renewal equation in octave differences,
D(k) = D_0(k) + Σ_W r_W D(k/p_W) (shifts c_W neglected; D at non-dyadic points by local cubic interpolation
of log(k D) in log k; weights 0.65 satellite, 0.90 primitive, periods to 16 plus a p^−3 tail), using M's
depth-8 octave differences for k = 2^12 … 2^21 and the §5.2 octaves beyond:

| k | 2^12 | 2^14 | 2^16 | 2^18 | 2^20 | 2^22 | 2^23 |
|---|---|---|---|---|---|---|---|
| k·D_M | 1.809 | 1.698 | 1.597 | 1.509 | 1.441 | 1.381 | 1.341 |
| k·D_0 | 0.903 | 0.839 | 0.788 | 0.747 | 0.720 | 0.693 | 0.668 |
| copies' share | 0.501 | 0.506 | 0.506 | 0.505 | 0.501 | 0.498 | 0.502 |

(With all weights 1 the share is 0.77, equally flat.)  The copies' share is constant, so they inherit M's drift
rather than causing it: the slow departure from a/k lives in T_0, the neutral-boundary shells (§3) and
decorations outside all maximal copies.  Simple forms for it do not fit: a/(k (ln k)^ν) + c/k² is best at
ν = 0.4 with χ² 67 (k ≥ 2^16, 31 dof; 400 from 2^14), and adding b/k drives ν to 1, i.e. back to
a/k + b/(k ln k) + c/k² (χ² 45).  So the renormalization structure is quantitatively understood, but the
tail's functional form reduces to T_0's, which is not yet known.

*Systematics checks (depth 8, GPU).*
- **Absolute accuracy against Grönwall.**  `fattened f-k27.npy 8 22` gives the exact F(2^−k) for k = 8…22 from
  the 2^27 coefficients (omitted weight e^{−2^(28−k)} ≤ 1.6e-28; 4 s).  The tree at the same k differs by
  −2.4e-8 … +2.2e-8 with errors 6e-9 … 1.7e-8, all within 2.2σ, χ² 17.1 over 15 (ignoring correlations
  between thresholds): the whole pipeline (tree, certificates, leaf sampling, estimator) is unbiased at the
  1e-8 level, 100× tighter than the §5 validation of the older Monte Carlo.
- **Certification margin.**  Safety 2, 4 and 8 give bit-identical A(2^14) and A(2^20) and identical errors,
  although leaves number 1.08e9, 1.29e9 and 1.56e9: a cell certified at one margin but not another has its
  certified disk covering it, so its 16 samples are unanimous and agree with the certificate, adding the
  same area and zero variance.  The margin never changed the answer, only the cost: safety 2 runs in
  31.2 s against 36.5 s (1.17×).
- **Leaf Newton settings.**  With the same leaves and samples, margin 1e-4 and tolerance 1e-12 (tighter than
  the defaults 1e-6, 1e-10) reproduce A(2^20) bit for bit; margin 1e-9 and tolerance 1e-8 (looser) raise it
  by 2e-10, i.e. a few false interior certificates.  The defaults sit two orders of magnitude inside the
  safe range.

*Careful universality test (`--tiles`).*  Per-tile octave differences from single runs (depth 8,
max_iter 2^22), with each tile's share s_t(k) = D_t(k)/Σ_t D_t(k) fitted as ln s_t = α + β log₂ k:
- Whole M, 16 × 16 tiles (0.16 × 0.075), k = 2^14 … 2^21: 61 tiles, share-weighted rms β = 0.0014 per octave,
  largest tiles within ±2.4σ, χ² 173 over 61 (a mild excess, possibly from errors treated as independent
  across k).  At this scale the shares are constant to ~0.1% per octave.
- Golden-mean cardioid box, 8 × 8 tiles (0.002): shares move strongly, but the drift dies away:
  rms β 0.33 per octave over k = 2^11 … 2^15 and 0.016 over 2^17 … 2^21; the biggest tiles fall from
  +0.03 … +0.57 to |β| ≤ 0.007 (still ~10σ), small tiles to −0.01 … −0.09.
So universality is asymptotic, not exact at finite k: pre-asymptotic corrections are large at small scales
and decay with k (by ~20× over six octaves at scale 0.002), consistent with faster-decaying local pieces
such as bulb-root horns (T ~ k^−3) inside small tiles.  The weak-convergence form (U) survives; an exact
finite-k form F(k) · w(U) does not.

*Deep-run calibration (64-bit orbit counters, `-DMANDELBROT_ORBIT64`).*  One H200:
- The 64-bit build costs 9% at depth 8 (39.1 s vs 35.8 s) and reproduces the 32-bit estimates bit for bit.
- Depth 5, max_iter 2^33, octaves 2^20 … 2^32: 1776 s for 6.6e11 iterations, i.e. 3.7e8 it/s, three orders
  below normal throughput.  Deep runs are latency-bound: each of the 8 batches waits ~220 s for its longest
  orbits, single lanes stepping ~9e6 it/s through 2^31-step orbits.  The octave counts at 2^28 … 2^32
  (11, 3, 0 samples) match the 1/k extrapolation from lower octaves (≈ 7, 3.6, 1.8), so no sign of double
  orbits stalling, but depth 5 is far too coarse (one sample weighs 3.7e-10) to measure the deep tail.
- Double-double vs double (paired, depth 4, max_iter 2^31): consistent at every k through 2^30 (errors
  2.5e-9 … 1.5e-8, 144 flips of 1.9e8 samples), extending the rounding check from 2^28 to 2^30 at this
  precision.
- Projection: a depth-8 run to 2^32 needs ≥ 10 batches (the 2^31-sample limit), each waiting ~220–1000 s for
  its longest orbit, so 1–3 hours of latency regardless of GPU speed, unless batches overlap (feeding the
  next batch's items to idle lanes while long orbits finish).  With the absolute-error budget of the
  previous paragraphs, two more digits by direct measurement stays ~10× over the 1–2 GPU-day budget; one
  more digit (±3e-10, measuring to 2^26–2^28 where double is validated) fits in about a day.

*Deep runs: batch overlap and a CPU tail.*  Two engine changes, results bit-identical throughout:
- `--overlap K` runs K batches at once in host threads with their own non-blocking streams (balanced batches
  after a probe).  Alone it gave 1.09× (depth 8, 2^20) and 1.23× (depth 8, 2^26), but made the deep run 2×
  slower: its time is single long orbits, which get slower still when they share SMs.
- The CPU tail: once at most MANDELBROT_CUDA_CPU_TAIL (default 1024) parked orbits remain, they are finished
  on the host's CPU threads (a core steps a sequential orbit ~40× faster than a lone GPU lane) and written
  back by a finish kernel.

| run | before | CPU tail | CPU tail + overlap 4 |
|---|---|---|---|
| depth 5, max_iter 2^33, octaves to 2^32 | 1776 s | 95.5 s | 57.5 s |
| depth 8, max_iter 2^26 | 203.4 s | 157.2 s | 158.0 s |
| depth 8, max_iter 2^20 | 38.4 s | 38.0 s | |

Deep runs are no longer latency-bound, which changes the two-digit budget.  The cheapest split is
μ = A(2^20) + [A(2^32) − A(2^20)] − T(2^32), with the middle term paired on the same samples (its variance
comes only from samples between the thresholds: σ ≈ 3.5e-9 at depth 8, against 2.0e-8 for A itself), and
T(2^32) ≈ 6e-10 from a model to ~5%.  A(2^20) to 2e-11 needs ~43 H200-hours as before; the paired deep
difference to 1.5e-11 needs ~6 more depth levels over depth 8, ~1400 × a depth-8 deep run (a few minutes)
≈ 100 H200-hours.  So two digits is ~6–7 GPU-days (down from ~12), still above the 1–2 day budget, while
one more digit (±3e-10) is a few GPU-hours.

*GPU engine: utilization (2026-09-29).*  Depth 8 at 2^20 went 24.9 → 14.1 s and depth 9 84.8 → 47.3 s, with
estimates bit-identical throughout (1.5065945640, 1.5065945785).  Sampling now runs at ~1.3e12 iterations/s,
~35% of the H200's FP64 peak (9 flops per step, FMA counted as 2); `Orbit::run` alone reaches 21.3 Tflop/s
(63%), against 27 for a bare FMA loop.  The changes, by effect:
- The atom-domain candidate is recomputed at the first Newton step (from z_0, same arithmetic), not tracked
  by selects every step (27% of a step); `Orbit::run` became a tight loop of 8-step fast blocks with the rare
  checks behind one branch.
- Refills come from a per-warp shared-memory queue of started samples (32 started at once): a lane's start
  had run with the rest of its warp idle.
- Freed device buffers stay in the stream-ordered pool: by default it returned them at each sync, so every
  batch mapped gigabytes afresh (leaf reductions went 771 → 67 ms).
- Parked orbits settle sorted by Newton's period (a counting sort), and leaf Newton gives up once
  |(f^p)'(w)| > √2 (74% of attempts converge to repelling cycles, 85% of Newton iterations).
- interior_distance takes the minimal period from the converged root instead of a Newton solve per divisor
  (in 13 of 400k centers the old loop used a non-minimal period, where the Koebe bound does not hold).
- Overlapped batches (`--overlap 2`, now the default) with adaptive batch sizes.

Tried and dropped: Newton from the best-returning orbit point (no gain on real leaves); shorter bursts; a
refill queue for resumed orbits (hoards work at round ends); cutting resume rounds short when warps go idle
(a round's tail is its longest orbits' serial remainder, so carrying them adds rounds); center hints for an
earlier first Newton (−40% center iterations, no time: center levels are bound by per-run latency, which
overlap hides instead).  The largest experiment: a persistent engine with each lane's orbit a state machine
and settles grouped on device queues (segmented fetch-and-add MPMC queues, stress-tested), removing the host
rounds.  It was correct but slower (depth 8 19.6 s against 15.1 s): each queue transition is a chain of
dependent L2 round trips that stalls its warp, and at the occupancy the registers allow (~6 warps per
scheduler) the stalls are not hidden, where the rounds do the same transitions in bulk.  Lock-free queues
also needed five fixes (CAS contention, hoisted spin loads, wraparound and segment-sizing deadlocks, a stale
nonempty bit) before the stress test passed.  The code is in git (0310d03 to dd1fbc1).

*Russian roulette for deep runs (2026-09-30).*  In the paired deep difference nearly all the cost is long
orbits, while the variance comes from the first octaves past the split.  So at Newton decision steps
d_j = first_newton · 2^j + 8, an unsettled leaf orbit survives with probability 1/2 (a hash of c and n) and
stops otherwise (`--roulette-from`, `--roulette-stride`, `--roulette-log2`).  Only escapes are reweighted: a
sample's value at threshold k past the reference r (the last threshold before any decision) is
b_r − W (b_r − b_k), where W = 2^(decisions survived), or 0 if killed.  Killed samples count as not escaped,
so interior orbits add no variance, and the sums stay integers.  With escapes per octave falling like 2^-j,
keeping 1/2 every octave adds constant variance per octave (as observed); every other octave (stride 2)
makes both the added variance and the cost fall like 2^-j/2.  Base 200, first Newton 512, A(2^14) − A(2^24)
on CPU (base cost 5.13e9 iterations; the columns are the cost past 2^14):

| roulette | iterations past 2^14 | std err | err² · cost vs plain |
|---|---|---|---|
| none | 6.57e9 | 1.30e-6 | 1 |
| from 2^14, stride 1 | 0.51e9 | 4.15e-6 | 1/1.3 |
| from 2^14, stride 2 | 1.06e9 | 2.25e-6 | 1/2.1 |
| from 2^15, stride 2 | 1.54e9 | 1.80e-6 | 1/2.2 |
| from 2^14, stride 4 | 1.94e9 | 1.92e-6 | 1/1.55 |

Over 48 seeds (from 2^15, stride 2) the mean differs from plain by −5.8e-8 ± 1.8e-7, and the seed spread
(1.89e-6) matches the reported error.  The gain grows with the number of octaves (plain cost is linear in
them, roulette's bounded), so over 2^20 … 2^32 it should be larger than 2.2×.  In the shallow part
roulette does not pay: its error comes from the early octaves.  `--prec dd` runs leaf orbits in
double-double alone, so deep runs can combine the two.  On CPU, double-double is 6.3× slower than double
here, 1.45× of it from more iterations: double settles slow interior orbits once they fall into an exact
floating-point cycle, which takes double-double longer to reach.

On one H200 (depth 6, `--prec dd`, octaves 2^20 … 2^32, roulette from 2^21 with stride 2): 855.6 s → 148.7 s
(5.75×), with iterations 3.83e12 → 1.81e12 and throughput 2.3e9 → 7.2e9 it/s, since the long orbits that
left lanes idle are mostly gone.  The paired deep difference A(2^20) − A(2^32) has std err 1.57e-8 → 2.16e-8
(summing octave variances: a sample escapes in one octave only), so variance 1.88× and err² · time 3.1× better.
The two estimates, 2.718e-6 and 2.697e-6, differ by 1.5σ of the added noise.  Past 2^28 single escapes weigh
2^4 … 2^6 samples, so at depth 6 some octaves have 0 or 1 escape and their std errors mean little; the
production run's ~4^8× more samples leave plenty.

*Double vs double-double at depth 8 (`--prec comparedd --flip-stats`, max_iter 2^32, two seeds, 4 shards
each, ~5 H200-hours per seed).*  Of 2.07e10 samples, 2.59e7 (0.13%) flip, and all but ~140 per seed escape in
both precisions at different steps: rounding re-randomizes the escape times of long orbits completely (log2
of the escape-time ratio spreads flat out to ±4 octaves; ~77% of the orbits escaping near 2^20 cross a
different threshold).  The ratio histogram is symmetric: 12,796,019 earlier vs 12,798,619 later (+0.5σ,
χ² 5.8/16 over mirrored quarter-octave bins) for seed 1, and 12,799,700 vs 12,797,156 (−0.5σ, χ² 7.1/16) for
seed 2.  Seed 1's A_dd − A_double was +(2–9)e-9 at 2^20 … 2^26 (2–3σ, correlated across k), but seed 2 gives
−(0.2–5)e-9 there, so it was a fluctuation.  The two-seed means are consistent with zero at every k: largest
+9.7e-10 ± 3.6e-10 at 2^26 (2.7σ, the largest of 19 correlated values), −(1–2)e-10 ± (1–2)e-10 at 2^27 …
2^31.  Stalls cancel too: double escaped where double-double hit max_iter 56 and 74 times, the reverse 62
and 65.  One sample per seed was interior in double and escaped in double-double (5.7e-12 each).  Interior
there includes exact floating-point cycles, which say nothing about c, so these are double orbits stuck in
a cycle rather than false Newton certificates.  So there is no evidence of a precision bias, but the bound is
only ≈ 3e-9 at 2^20 and 3e-10 at 2^30 (2σ), far above 3e-11: two digits rest on the symmetry, and a direct
check at that level is out of reach.

*Deep part: depth scaling with double-double and roulette.*  Single H200s, `--prec dd`, octaves 2^22 …
2^32, roulette from 2^23 with stride 2 (the paired difference A(2^22) − A(2^32), summing octave variances):

| depth | time | iterations | samples | A(2^22) − A(2^32) | err² · time |
|---|---|---|---|---|---|
| 8 | 1235 s | 2.94e13 | 2.07e10 | 6.540e-7 ± 2.68e-9 | 8.9e-15 |
| 9 | 4493 s | 1.14e14 | 6.93e10 | 6.523e-7 ± 1.33e-9 | 7.9e-15 |
| 10 | 15411 s | 4.44e14 | 2.34e11 | 6.519e-7 ± 6.6e-10 | 6.8e-15 |

err² · time improves only 1.12× and 1.17× per level (1.9× for the shallow part): variance falls 4× per level,
but time rises 3.4–3.6×, in proportion to iterations.  So time is no longer latency: it is every leaf sample's
orbit below 2^22 (1400–1900 steps on average, in double-double at ~2.9e10 it/s), although only samples
escaping past 2^22 contribute.  Deep leaf samples are plain Monte Carlo: the tree certifies nothing near
them, so depth buys only more samples.  At depth 10's efficiency, 1.5e-11 would take ~8,000 H200-hours.

*Levers for the deep part (2026-09-30).*
- Leaf-level importance sampling, measured offline: `--dump` of every sample in three base rows at depth 10
  (37.5M leaves, double, to 2^32), samples 0–7 as the pilot and 8–15 as the outcome, cost-weighted Neyman
  allocation across leaf classes, cross-validated on held-out leaves.  Best gain 1.3× (classes by the pilot
  samples alive at 2^14–2^18); finer classes overfit the 5,417 deep escapes and lose.  As before (two-phase
  importance sampling, above), deep escapes are scattered at random among near-boundary leaves.
- Double and double-double carry no shared information on long orbits (their escape times are re-randomized),
  so double cannot serve as a control variate for a double-double deep run.  But if double's bias is small
  enough, the deep run is unnecessary: the shallow run already computes every sample to 2^22 in double, and
  continuing its survivors to 2^32 with roulette adds 0.24× its iterations (0.96× without roulette; from the
  dump).  The deep difference's error falls 4× per level in variance, reaching ~1e-11 at depth 16.
- The fused double-double step (41 FP64 ops against 84; double 6) halves double-double's arithmetic, should
  it be needed.

*Precision ladder (depth 5, base 1000, 5.9e8 samples, octaves 2^14 … 2^26, CPU).*  Double rounded to b bits
against double, by outcome (A_b(k) − A_double(k) contributions):

| bits | escaped → interior | interior → escaped | escaped ↔ max_iter | re-timed escapes at 2^20 |
|---|---|---|---|---|
| 24 (= float, identical) | 34,288 (+1.26e-5) | 4,655 | 120 / 79 | −3.8e-7 |
| 27 | 7,258 (+2.66e-6) | 1,008 | 103 / 81 | −1.8e-7 |
| 30 | 1,558 (+5.7e-7) | 235 | 73 / 73 | −4.3e-8 |
| 36 | 76 (+2.8e-8) | 10 | 68 / 77 | +2.7e-8 (noise) |

Two mechanisms, both vanishing fast with precision:
- Stuck orbits: the low-precision orbit falls into an exact floating-point cycle (Brent), so an exterior point
  counts as interior.  Counts fall 2^-0.73 per bit over 24 … 36 bits; extrapolated to 53 bits, ~0.5 per
  2.07e10 samples, i.e. +3e-12 in area (the double-double runs saw 0 and 1 per seed).
- Re-timed escapes: as a fraction r of D(k), the bias depends only on log2(k) − bits, i.e. on k·ε (curves for
  24, 27, 30 and 36 bits collapse): r = −1.8% at −10, −4…−8% at −9, rising to a −25…−35% plateau by −4, and
  |r| ≲ 0.5% (zero within errors) from −13 to −22.  Double sits at log2(k) − 53 = −31 (2^22) … −21 (2^32), so
  even decaying only 2× per octave past −22, the bias at 2^22 is < 3e-12.
Symmetric outcomes (escaped ↔ max_iter) cancel at every precision.  So double's precision bias is a few
1e-12, inside the 3e-11 target, and two digits need no double-double.

*Two-digit budget for one double run to 2^32 with roulette.*  Depth scaling on CPU, the same four base rows
(shard 7/250 of base 1000), base-only runs to 2^22 against full runs to 2^32 (roulette from 2^23, stride 2):

| depth | A(2^22) err | base time | err² · t | full / base: time, iterations | deep difference err |
|---|---|---|---|---|---|
| 9 | 5.13e-10 | 22.7 s | 5.97e-18 | 1.90, 1.24 | |
| 10 | 2.00e-10 | 84.8 s | 3.39e-18 (1.76×) | 1.37, 1.20 | 4.3e-11 |
| 11 | 7.80e-11 | 322.7 s | 1.96e-18 (1.73×) | 1.37, 1.23 | 2.17e-11 |

So the shallow gain holds at ~1.75× per level (GPU: 2.05×, 1.82×), and continuing to 2^32 costs 1.2–1.4×.
But the deep difference converges more slowly: its variance falls 3.9× per level against 6.55× for A(2^22),
so by depth 15 it carries ~40% of the variance.  From the GPU's depth 10 (A(2^22) ± 3.08e-9 in ~245
GPU-s; deep difference ± 6.8e-10 scaled from the CPU rows), time ×3.7 per level, one run at effective depth
d = depth + log2(base / 1000):

| statistical target | effective depth | H200-hours (overhead 1.24–1.37) | if A's gain tapers to 1.6× |
|---|---|---|---|
| 2.7e-11 | 15.34 (depth 15, base 1266) | 91–101 | 120–132 |
| 3.0e-11 | 15.21 (depth 15, base 1159) | 77–85 | 101–112 |

With the tail model (±6e-12 … 1.2e-11) and precision (a few 1e-12), the total error is ~3e-11.

*Taylor-series ideas (shelved).*
- Skipping iterations by linearizing each leaf's orbit (series approximation): samples of a leaf share at most
  the prefix before the first of them ends, 22–28% of iterations, so ≤ 1.36× even before accuracy limits.
- Control variates from the linearized multiplier map: an interior sample's λ and dλ/dc model its hyperbolic
  component as {c : |λ_i + λ'_i (c − c_i)| < 1}.  For sample j, ĥ_{−j} is the union of the other samples'
  models, and mean_j h_j + β · mean_j (∫ĥ_{−j} − ĥ_{−j}(c_j)) (∫ over j's stratum) is unbiased for any
  fixed β, since c_j is uniform in its stratum and independent of the others.  On 50k mixed leaves of depth 10 (two independent replicas each; these leaves carry
  91% of A(2^22)'s variance): variance ratio 0.564 at β = 1, 0.490 at the optimal β = 0.72 (2.0×).  Only 4.8%
  of samples are mispredicted, but the model itself is noisy (which small components 16 samples hit), and
  trusting models only near their samples is worse (τ = 0.5: 0.597).  The deep difference (~43% of the
  variance at depth 15) is untouched, so the overall gain is ~1.4× at best before the cost of evaluating models
  (~4000 ops per sample) and reworking error bars.  Pooling models across neighbouring leaves might reach
  3–5× on A(2^22), still capped at ~2.3× overall.
- Leaf first Newton step: 8192 → 2048 saves 7.7% of iterations (5% CPU time) at depth 10, estimates
  identical; on the GPU (depth 15 pilot below) it is 23% slower.

*Engine work for the depth-15 run (2026-10-01).*  The CPU-based budget missed GPU latencies; each was found
by a pilot (base rows 0 and 633 of base 1266 at depth 15, the production command, one H200) and fixed with
results bit-identical (unit tests, and identical pilot estimates throughout: 4.2540360431e-03 ± 1.45e-12):
- Deep queue (`--deep-from 2^23+8`): continuing to 2^32 made each batch wait ~11.5 s for its few long orbits
  (depth 8: 10.3× the time for 1.23× the iterations).  Samples unsettled at deep_from are suspended; their
  leaves' words are saved and zeroed in the batch, and deep passes of 2^23 states from many batches finish
  them at full width, then reduce the saved leaves with the same integer sums.
- Memory: tree levels may reach 2^31 cells at ~90 bytes per cell (center states 160 bytes, parked twice), so
  `--max-level-cells 2^27` splits batches, and a level too large even for one base cell (2.6e8 cells on the
  real axis) has its cell list split into independent subtrees continued from that level.  The device pool
  is trimmed and the allocation retried once on out-of-memory.
- Tree overhead: with every batch building all 16 levels, tree time equalled sampling time (1958 vs 2242
  worker-seconds over 1281 batches; the centers' 4e13 iterations are ~40 s of GPU work, so this is latency:
  ~95 ms per level step).  `--split-depth` builds the levels below it once for runs of base cells and batches
  the cells there.
- Storage: the PVC is ReadWriteOnce, so jobs on different nodes wait for each other while holding GPUs;
  production shards use node-local /data (the conda environment builds in ~15 s) and return results in their
  logs, merged locally.

| pilot (2 rows, 3.12e10 leaves, 1.28e15 leaf iterations) | time | tree (worker-s) | peak GPU memory |
|---|---|---|---|
| batch 2^25, overlap 2 | 2336.7 s | 1958 | 48 GB |
| batch 2^26, overlap 3 | 2103.5 s | 3141 | 113 GB |
| split depth 8 | 2027.4 s | 1001 | 41 GB |
| split depth 11 | 1861.5 s | 689 | 38 GB |
| split depth 11, first Newton 2048 | 2296.2 s | 744 | 38 GB |
| split depth 11, batch 2^26 | 2025.1 s | 1367 | 72 GB |
| split depth 11, overlap 3 | **1633.7 s** | 1072 | 57 GB |
| split depth 12, overlap 3 | 1656.7 s | 1168 | 55 GB |
| split depth 13, overlap 3 | 2071.7 s | 2267 | 53 GB |

The pilot rows hold 2.1× the average leaves (production: ~9.4e12), so the whole run is ~1634 / 2.1 × 633 ≈
4.9e5 s ≈ 136 H200-hours, sampling at ~7.8e11 it/s overall (~21% of FP64 peak).  Engine timing shows the main
sampling pass using ~55% of warp lanes, the next place to look.  The two rows alone give A(2^32) ± 1.45e-12,
~2.5e-11 for the whole run.

*Result (2026-10-05).*  The production run `prod-d15` (depth 15, base 1266, seed 3, double, to 2^32 with roulette
and the deep queue, split depth 11, overlap 3; 32 shards of 4.6–5.2 h on single H200s, ~158 GPU-hours;
1.1e13 leaves, 1.77e14 samples, 4.4e17 iterations):

| k | A(k) | std err |
|---|---|---|
| 2^14 | 1.5067932403 | 1.47e-10 |
| 2^20 | 1.5065945847 | 2.60e-11 |
| 2^22 | 1.5065925360 | 2.14e-11 |
| 2^26 | 1.5065919225 | 2.21e-11 |
| 2^32 | 1.5065918842378 | 2.33e-11 |

k·D(k) for k = 2^22 … 2^31: 1.3876 1.3655 1.3463 1.3293 1.3150 1.3024 1.2921 1.2815 1.2736 1.2711 (errors up to
0.0044).  Fitting a + b/j and a quadratic over the last 6, 8 and 10 octaves gives T(2^32) = 5.71e-10 … 5.97e-10,
median 5.843e-10, model spread ±1.3e-11, fit error ±7e-13 (`analysis/final_tail.py`).  So

  μ = 1.506591883654 ± 0.000000000027 (statistics 2.3e-11, tail model 1.3e-11, precision a few 1e-12)

against the published 1.5065918849: 1.25e-9 lower, within its stated uncertainty, with ~100× smaller error.
The dominant systematic is the tail extrapolation, which backtests on these octaves put at 1–2%.

*Tail backtest on the production octaves* (`analysis/tail_backtest.py`): cut at 2^c, fit k·D(k) over the n
octaves before it, and predict A(2^c) − A(2^32), measured to 3–8e-12.  At every cut a + b/j underpredicts and
the quadratic overpredicts, the bracket widening with n; for n = 6, a + b/j is −1.00, −0.99, −0.95, −0.81,
−0.78% and the quadratic +0.50, +0.42, +0.35, +0.34, +0.22% at c = 26 … 30.  Correcting the n = 6 fits at 2^32
by these biases (a + b/j by ~+0.75%, the quadratic by ~−0.2%) gives T(2^32) = 5.85–5.87e-10, so the tail is
known to ~±5e-12 and μ = 1.506591883652(24).

*Direct tail check* (run tail-d13: depth 13, seed 4, max_iter 2^36, 8 shards, ~26 GPU-hours;
`analysis/tail_direct.py`).  Independent samples give A(2^32) = 1.5065918842169 ± 1.4e-10, 2.1e-11 (0.15σ) from production, and k·D(k) agrees with
production on every overlapping octave 2^26 … 2^31.  Past 2^32 it measures k·D(k) = 1.234(24), 1.284(50),
1.173(67), 1.319(142) at 2^32 … 2^35, so the octaves to 2^36 sum to 5.436e-10 ± 9.9e-12; with ~3.6e-11 beyond
2^36 (f̄ = 1.22 … 1.27), T(2^32) = 5.80e-10 ± 1.0e-11 directly, against 5.86e-10 ± 5e-12 extrapolated (0.5σ).
Combined, T(2^32) = 5.847e-10 ± 4.5e-12, and

  **μ = 1.506591883653(24)**  (statistics 2.3e-11, tail 4.5e-12, precision a few 1e-12), i.e.
  μ = 1.506591883653 ± 4.7e-11 at 95% confidence (the convention of Förstemann's and Lo's estimates).

*Comparison with earlier estimates.*  Förstemann (2012, GPU pixel counting, 20 random grids of 2097152², max
iterations ~8.6e9): 1.5065918849 ± 2.8e-9 as a 95% interval, i.e. σ ≈ 1.4e-9; he counted unresolved pixels
(~8e-10 of area, incl. ~3e-10 of exterior points escaping after his limit) as members, so corrected
≈ 1.5065918846, 0.7σ from ours.  Hsing Lo (Oct 2025, github.com/hsingtism/mandelbrot-area, one random point
per grid cell, 30.4e12 points): 1.5065918902 ± 5.4e-9 (95%), +1.9σ from ours when its two groups are combined by
inverse variance (group A 12 grids +0.7σ, group B 99 grids +2.1σ).  Its membership test is loose (consecutive
iterates within 2^-16, though the site says 2^-48; Floyd cycle test within 2^-32; 2^34 iterations).  Measured
against ours (`analysis/lo/`): its false members are ~2e-10 of area (exterior bands along the main cardioid and the
period-2 disk, sampled decade by decade in |λ| − 1; nothing near the cusp or elsewhere above ~1e-11), while its
undecided points (6.5e-10 of samples, 3.7e-9 of area, ≥ 96% interior since exterior survivors past 2^34 are only
~1.5e-10) are dropped from the denominator, biasing it low by ~2.8e-9.  Corrected, Lo ≈ 1.5065918917 ± 2.9e-9:
2.8σ above ours and 2.2σ above corrected Förstemann, so most likely a fluctuation in its group B.

## 6. Reproducing

```
meson compile -C build/release
./build/release/hyperbolic 16 1024                      # areas, ∫|c'|², perimeters; ~2 min on 18 cores
./build/release/octscan f-k27.npy 26                    # local exponents; ~10 s
./build/release/census f-k27.npy 26 20 4                # root census and wake ownership; ~12 s
./build/release/renorm angle_map.npy                    # wake periodicity and pullback tests
./build/release/rootsum f-k27.npy 8 26 20 4 7           # octave energy near roots, by period
./build/release/hyperbolic 16 1024 1 components.txt     # per-component centers and areas
./build/release/tuning components.txt 16                # tuning structure and D(s)
./build/release/fattened f-k27.npy 8 22                 # exact fattened-set areas for k ≤ 22
./build/release/escape_area 16000 1048576 SEED 16384 ... 1048568   # fattened areas; ~22 s per 2e9 samples
./build/release/escape_tree --prec compare 64 1024 16384 262144 1048568   # certified tree, float vs double; ~85 s
```

The area result is reproduced from the committed shard results by `python analysis/tail_direct.py` (see
`analysis/README.md`); the runs themselves are `cluster/idle_shards.py analysis/results/*/run.json` at tag
`area-run-7`.

Unit tests cover everything the analyses use: `meson test -C build/release angles octaves numpy tests`.
They check Lavaurs' algorithm against A000740 and the satellite counts through period 18, the fast
algorithm against a quadratic reference, tuning and untuning round trips, wake ownership against brute
force, octave spectra against a direct DFT, and `.npy` round trips.

## Open directions

1. Formalize §4 in Lean, in `ray`, which already has Grönwall's theorem for M (`multibrot_volume_sum`).
   Then strengthen it to (log N)^-3 via the satellite root −3/4.
2. Compute C, the lingering constant, from the dynamics just outside neutral boundaries.  Then Q predicts
   K with no fitting, which would be a real test of the 1/log N mechanism.
3. Use the 1/log N law with measured K to sharpen μ, e.g. by fitting only a and the parity oscillation
   with μ tied to K.  Check against Siudem & Świątek's independent ~2^25 coefficients.
4. Push the hyperbolic lower bound to periods 17–18 (about 5 and 20 minutes).
