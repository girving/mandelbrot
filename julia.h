// Areas of filled Julia sets K(c) from the area transfer operator
//
// For f(z) = z^2 + c, the operator
//   L g(z) = Σ_{f(w) = z} g(w) / |f'(w)|^2 = (g(√(z - c)) + g(-√(z - c))) / (4 |z - c|)
// pulls area back: the area of f^-n(E) is ∫_E L^n 1.  On an annulus A = {r1 < |z| < r2} with
//   r1 - |c| > r1^2    (so f maps V = {|z| ≤ r1} into itself, inside K, and preimages of A avoid V)
//   r2^2 - r2 > |c|    (so points with |z| ≥ r2 escape, and preimages of A lie inside |z| < r2)
// every point of A either escapes, through X = {z ∈ A : |f(z)| ≥ r2}, or stays in K.  So
//   area K(c) = π r2^2 - ∫_X h,   h = (1 - L)^-1 1,
// and since f^-1(A) ⊂⊂ A, L is compact on functions analytic near A and h is analytic: the fractal Julia set
// enters only through L's spectrum, so discretizations converge exponentially.  Both conditions need |c| < 1/4.
//
// Discretization: collocation at Chebyshev points (first kind) in s = log r times equispaced angles (odd count),
// evaluating h at preimages by barycentric interpolation.  Both preimages of z have radius √|z - c|, so they share
// one Chebyshev row and their angular rows add: L is stored as N × nr and N × nt factors.  (1 - L) h = 1 is solved
// by iterative refinement: GMRES in double for corrections, residuals in the working precision S.  Points,
// weights and quadrature nodes come from arb, rounded to S.  (A two-grid Atkinson-Brakhage preconditioner cut
// GMRES iterations ~3x near c = 1/4 but cost more than it saved: L only halves frequencies per application.)
#pragma once

#include <cstdint>
#include <vector>
namespace mandelbrot {

using std::vector;

struct JuliaParams {
  double cx = 0, cy = 0;  // c, exactly these doubles
  double r1 = 0, r2 = 1.5;  // Annulus; r1 = 0 means max(2|c|, 1/4)
  int nr = 32, nt = 101;  // Chebyshev points in log r, angles (odd)
  double grade = 0;       // Angular grading a ∈ [0, 1): points pack near θ = 0 by (1 + a)/(1 - a)
  int nq = 0, ng = 0;     // Area quadrature: angles (0: 2 nt + 1), Gauss-Legendre points in r (0: nr + 16)
  int pre = -1;            // Parabolic preconditioner (1 - L+)^-1: 1 on, 0 off, -1 when |f'(q)| is near 1
  int stencil = 6;         // Its local interpolation width, on a grid oversampled by
  int oversample = 2;
  bool fast = false;       // Corrections with a fast approximate L (NUFFT style), O(N (nr + nt)) rather than O(N^2)
  int fast_width = 16;     // Its kernel width at 2x oversampling
  int max_refine = 20;
  bool eig = false;  // Also estimate L's leading eigenvalue (power iteration in double)
  bool verbose = false;
};

template<class S> struct JuliaResult {
  S area;
  double residual = 0;  // Final max |1 - (1 - L) h| at the collocation points
  int refinements = 0, gmres_iters = 0;
  double rho = 0;  // L's leading eigenvalue e^{P(2)}, if eig
  double secs = 0;
};

template<class S> JuliaResult<S> julia_area(const JuliaParams& p);

// Barycentric interpolation rows, exposed for testing.  cheb_row: at s in [s0, s1] with nr first kind points.
// trig_row: at angle θ given as η = e^{iθ/2} (either sign), on θ_j = 2πj/nt for odd nt.
template<class S> void cheb_row(const vector<S>& nodes, const vector<S>& weights, S s, S* row);
template<class S> void trig_row(const vector<S>& rho_x, const vector<S>& rho_y, S ex, S ey, S* row);

}  // namespace mandelbrot
