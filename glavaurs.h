// The Lavaurs model at the root of the cardioid's p/q bulb (double precision; lavaurs.h is the q = 2 special case)
//
// At c_0 = λ/2 - λ²/4, λ = e^{2πi p/q}, the fixed point α_0 = λ/2 has multiplier λ; in w = z - α_0 the map is
// f(w) = λ w + w², and f^q(w) = w + A w^{q+1} + … has q attracting and q repelling petals, permuted by f.  The Fatou
// coordinate Φ(f(w)) = Φ(w) + 1/q has the formal series Φ = Σ_{j=-q}^{N} a_j w^j + β L(w), L = log(w^q)/q on a
// per-petal branch (only a_0 is a gauge: a_0 = 0); the a_j are complex, solved in acb.  Φ_a: iterate into an attracting
// petal; Ψ_k: the repelling parametrization of petal k (Ψ_k(ζ + 1) = f^q(Ψ_k(ζ))).  The transit is g_σ(w) =
// Ψ_s(Φ_a(w) + σ), s = exit_petal(entering attracting petal) (a neighbour on the given side).  An r-transit component
// with final excursion n is a σ where R_σ = f^n ∘ g_σ^r ∘ f fixes the critical point w = -λ/2 with an attracting cycle.
//
// Calibration (checked at q = 2, 3 against M): for c in the limbs [CF(p/q), k] with α_c's multiplier λ e^{2πiε},
// σ = (1/ε + q² k)/q² + b, and a family's constant C = lim k⁴ area_M is (4π² sin²(πp/q)/q⁴) area_σ.
#pragma once

#include "complex.h"
#include <vector>
namespace mandelbrot {

struct GeneralLavaurs {
  typedef Complex<double> C;
  int p, q;
  C lam, A, v, crit;          // λ, f^q(w) = w + A w^{q+1} + …, critical value and point
  std::vector<C> a;           // a_j for j = -q .. N (index j + q); a_0 = 0
  C beta;
  int N;
  int side;                   // the gate: an orbit entering attracting petal k exits through the repelling petal whose
                              // axis is at k's axis + side π/q (the two sides of the root; q = 2 limbs k/(2k+1): +1)

  GeneralLavaurs(const int p, const int q, const int side = 1, const int N = 16);
  // The repelling petal paired with attracting petal k
  int exit_petal(const int k) const;

  // Φ, Φ', Φ'' by the series with L(w) = (log(w^q) + 2πi branch)/q
  void series(const C w, const int branch, C& s, C& d, C& dd) const;
  // The petal (kind = -1 attracting, +1 repelling) whose sector contains w, or -1; and the branch of L there
  int petal(const C w, const int kind) const;
  int branch(const C w, const int kind, const int k) const;
  // Attracting coordinate with derivatives and the entering petal (false if the orbit does not enter)
  bool phi_a(C w, C& s, C& d, C& dd, int& petal, const int max_steps = 1 << 20) const;
  // Repelling parametrization of petal k with derivatives
  bool psi(const C zeta, const int k, C& w, C& d, C& dd) const;

  // Center of the r-transit component with excursion n, from a guess; area_σ by multiplier boundary tracing, conv,
  // cusp = |σ'(1)|/|a_1| (≈ 0 primitive, O(1) satellite)
  bool center(const int r, const int n, C& sigma) const;
  bool area(const int r, const int n, const C guess, C& center, double& area, double& conv, double& cusp,
            C& a1, const int Nb = 64) const;
  // Multiplier-map point σ(μ) of a component
  bool multiplier_point(const int r, const int n, const C center, const C mu, C& sigma) const;
  // Θ_r(σ) = p_r - ζ0 (p_1 = ζ0 + σ, p_{i+1} = H(p_i) + σ), its σ-derivative and Π_{i<r} H'(p_i)
  bool theta(const C sigma, const int r, C& th, C& dth, C& d2th, C& hprod) const;
  C zeta0() const;
};

}  // namespace mandelbrot
