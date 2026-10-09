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

}  // namespace mandelbrot
