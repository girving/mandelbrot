// The Lavaurs model at the 1/2 root of the cardioid: phase-plane components and their areas
//
// In w = z + 1/2 the map z² - 3/4 is F(w) = -w + w² (fixed point 0 with multiplier -1); nearby parameters are
// F + δ.  Its Fatou coordinate Φ(F(w)) = Φ(w) + 1/2 has the asymptotic series
//   Φ(w) = 1/(4w²) + 1/(4w) + (11/8) L(w) + Σ_{j≥1} a_j w^j,
// L = log on one petal of each pair and log(-·) on its partner (exact rational a_j, computed here).  Φ_a is the
// attracting coordinate (petals around the real axis), Ψ₊ the repelling parametrization of the upper petal
// (Ψ₊(ζ + 1) = F²(Ψ₊(ζ))), Ψ₋(ζ) = F(Ψ₊(ζ - 1/2)).  The Lavaurs map g_σ = Ψ_s(Φ_a(·) + σ) exits through the repelling
// petal s of the same sign as the entering attracting one.
//
// A component with r transits and excursion n is a σ where the return map R_σ = F^n ∘ g_σ^r ∘ F of the critical point
// w = 1/2 has an attracting fixed point; its area comes from tracing the multiplier map's inverse σ(μ), |μ| = 1, and
// area = π Σ k |a_k|² for σ(μ) - σ(0) = Σ a_k μ^k.  The phase relates to the M family phase σ_M (δ = iπ/(2k + σ_M)) by
// σ = -σ_M/2 + 3πi/8 + 1/2 (excursion n = j - 1), so a family's constant C = lim k^4 area_M is (π²/4) area_σ.
#pragma once

#include "complex.h"
#include "expansion.h"
#include <string>
#include <vector>
namespace mandelbrot {

using std::vector;

template<class S> struct LavaursModel {
  typedef Complex<S> C;
  vector<S> a;  // a_1, a_2, ...: the series coefficients beyond the log term
  S c_log;      // 11/8

  explicit LavaursModel(const int terms = 30);

  // Φ, Φ', Φ'' of the series at w (attracting or repelling branch of L)
  void series(const C w, const bool attracting, C& s, C& d, C& dd) const;

  // Attracting coordinate with derivatives, and the sign of the attracting petal the orbit enters (false if it does
  // not reach a petal within max_steps or escapes)
  bool phi_a(C w, C& s, C& d, C& dd, int& petal, const int max_steps = 1 << 20) const;

  // Repelling parametrization Ψ_petal with derivatives (false if Re ζ is too far out)
  bool psi(const C zeta, const int petal, C& w, C& d, C& dd) const;
};

struct LavaursResult {
  bool ok;
  Complex<double> center;  // σ
  Expansion<2> area;       // In σ; the family constant is (π²/4) area
  double conv;             // Relative change of the area from the N/2 subrule
};

// The component with r transits and excursion n whose center Newton finds from guess
LavaursResult lavaurs_area(const int r, const int n, const Complex<double> guess, const int N = 64);

// Combinatorial labels (double precision).  The dynamical plane is cut by the imaginary axis in z = w - 1/2 (each half
// maps bijectively onto C minus (-∞, -3/4]), so a point's address is the sequence of sides (L: Re z < 0, R: Re z > 0)
// of its forward orbit.  A single-transit center σ (excursion n) is labeled by the address of x = Ψ(ζ0 + σ), whose
// n-th image is the critical point: n symbols.
std::string lavaurs_address(const int n, const Complex<double> sigma);

// Two-transit centers as island preimages: Θ(σ) = H(ζ0 + σ) + σ - ζ0, H = Φ_a ∘ Ψ the horn map; at a single-transit
// center σ_u (excursion n_u) Θ(σ_u) = σ_u - (n_u + 1)/2 and Θ'(σ_u) = 1, Θ is 2:1 on the island, and a two-transit
// center is a solution of Θ(σ) = σ_c + j/2 with σ_c a single-transit center (final excursion n_c - j).
// More generally Θ_r(σ) = p_r - ζ0 with p_1 = ζ0 + σ, p_{i+1} = H(p_i) + σ: the r-transit centers are the solutions of
// Θ_r(σ) = σ_c + j/2, and at an (r-1)-transit center (excursion n) Θ_r' = 1 and Θ_r = σ - (n + 1)/2 again, so they
// are tracked from (r-1)-transit centers the same way.
bool lavaurs_theta(const Complex<double> sigma, Complex<double>& theta, Complex<double>& dtheta,
                   Complex<double>& d2theta, const int r = 2);
// Θ(σ) = target by a homotopy in the target from the island center, on branch ±1 (the quadratic's two roots)
bool lavaurs_island(const Complex<double> center, const Complex<double> target, const int branch,
                    Complex<double>& sigma, const int r = 2);
// The island containing σ: Newton on H'(ζ0 + σ) = 0 (each island has exactly one critical point of the horn map, at
// its single-transit center: Φ_a' = 0 at the critical preimage and Ψ' ≠ 0); σ is replaced by that center
bool lavaurs_island_center(Complex<double>& sigma, const int r = 2);
// The branch label of a point σ in the island of a center with excursion n_u: the side of F^{n_u}(Ψ(ζ0 + σ)) (the
// orbit's passage by the critical point, where F is 2:1)
char lavaurs_island_side(const int n_u, const Complex<double> sigma);

}  // namespace mandelbrot
