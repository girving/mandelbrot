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
(`scratch/fattened_exact.py`).

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

## 6. Reproducing

```
meson compile -C build/release
./build/release/hyperbolic 16 1024                      # areas, ∫|c'|², perimeters; ~2 min on 18 cores
./build/release/octscan f-k27.npy 26                    # local exponents; ~10 s
./build/release/census f-k27.npy 26 20 4                # root census and wake ownership; ~12 s
./build/release/renorm angle_map.npy                    # wake periodicity and pullback tests
./build/release/escape_area 16000 1048576 SEED 16384 ... 1048568   # fattened areas; ~22 s per 2e9 samples
```

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
