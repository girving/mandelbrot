// Areas of satellite bulbs, batched on CPU threads or the GPU
//
// A job is the p/q satellite of a parent hyperbolic component W of period P, given by W's center.  W's multiplier
// map c_W(λ) is continued from the center to the root λ0 = e^{2πip/q} (explicit for the cardioid); the child's center
// solves f_c^{qP}(crit) = crit by Newton from c_r + λ0 c_W'(λ0)/q^2; its boundary c(μ), |μ| = 1, is traced by
// continuation in double and polished pointwise in Expansion<2> (simplified Newton: residual in E, Jacobian in
// double), and the area is π Σ k |a_k|² for the Taylor coefficients of c(μ) - c0 = Σ a_k μ^k (c0 the center, against
// cancellation), the a_k by DFT of the boundary values, which is spectral since c is analytic.  Cardioid parents use cusp coordinates ζ = z - 1/2, δ = c - 1/4 (ζ ↦ ζ² + ζ + δ), so bulbs
// near the cusp keep full relative precision.  Continuation is predictor-corrector (Euler along (z', c'), then
// Newton) with adaptive step splitting.  F = area q^4 / (π |c_W'(λ0)|^2) normalizes by the parent.
//
// Each job is one thread; the arithmetic is identical on CPU and GPU (no transcendentals on the device: twiddles
// come from the host, and -ffp-contract=off), so their results agree bit for bit.
#pragma once

#include "complex.h"
#include "expansion.h"
#include <cstdint>
#include <vector>
namespace mandelbrot {

using std::vector;
typedef Expansion<2> E2;

struct BulbJob {
  int P;                   // Parent period
  int p, q;                // Child rotation number p/q
  int shift;               // 1 for a cardioid parent (P = 1): cusp coordinates
  Complex<double> center;  // Parent center, in shifted coordinates if shift
  Complex<E2> lam0;        // e^{2πi p/q}
};

enum BulbStatus { bulb_ok = 0, bulb_parent = 1, bulb_center = 2, bulb_period = 3, bulb_area = 4 };

struct BulbResult {
  int status;
  Complex<double> center;  // Child center (unshifted)
  E2 area, F;
  double w;                // |c_W'(λ0)|^2
  double conv;             // Relative change from the N/2 subrule (every other point)
};

struct BulbParams {
  int N = 64;              // Boundary points
  int radial_steps = 4;    // Initial pieces from the center to the first boundary point (split adaptively)
  int substeps = 1;        // Initial pieces between boundary points (split adaptively)
  int polish = 3;          // Maximum Expansion<2> Newton steps per boundary point (adaptive)
  double accept = 1e-9;    // Final double Newton step accepted when roundoff stops it short of 1e-15
  bool cuda = false;
};

// Areas for all jobs (in job order)
vector<BulbResult> bulb_areas(const vector<BulbJob>& jobs, const BulbParams& params);

// A job for the p/q bulb of the parent with period P and (unshifted) center c
BulbJob bulb_job(const int P, const Complex<double> center, const int p, const int q);

}  // namespace mandelbrot
