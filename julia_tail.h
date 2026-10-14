// Escape-time tails of filled Julia sets K(c) at and near parabolic parameters, on CPU threads or the GPU
//
// Monte Carlo over starting points z: iterate z → z^2 + c until |z| > 2 (escaped; histogram the octave of the
// escape step) or z is certified in K, or max_iter.  Certification:
//   at a parabolic c: c has a cycle z_0, ..., z_{P-1} whose multiplier is a primitive q-th root of unity, and
//     near each point f^{Pq}(z_k + w) = z_k + w + a_k w^{q+1} + ..., so the attracting petals Re(-1/(q a_k w^q)) > R
//     (where f^{Pq} acts on the Fatou coordinate as u ↦ u + 1 + ...) lie in K;
//   near: c inside the main cardioid, |z - α| < (1 - |λ|)/2 maps into itself (α the attracting fixed point,
//     λ = 2α), so lies in K.
// Sampling:
//   global: uniform in the square [-2, 2]^2 (points outside |z| ≤ 2 escape at once);
//   local: square shells max(|x|, |y|) ∈ (s_{b+1}, s_b], s_b = r0 2^-b, around each cycle point, down to the scale
//     whose escape time reaches max_iter / 4, so every scale gets equal effort; shells weight their samples by area.
// Sample i is a pure function of (seed, i) and the arithmetic is the same on CPU and GPU (no transcendentals, no
// contraction), so their histograms agree exactly.
#pragma once

#include <complex>
#include <cstdint>
#include <vector>
namespace mandelbrot {

using std::vector;

enum class TailMode { global, local, near };

struct TailParams {
  TailMode mode = TailMode::global;
  std::complex<double> c;
  std::complex<double> z_guess;  // A point near the parabolic cycle (unused by near)
  int P = 1, q = 1;              // Cycle period and petals (unused by near)
  int64_t samples = 1 << 20;
  int64_t max_iter = 1 << 20;
  double r0 = 0;                 // Local: outermost shell half-side (0: 0.3 times the germ scale |a|^{-1/q})
  double R = 50;                 // Petal certification threshold
  uint64_t seed = 1;
  bool cuda = false;
};

struct TailResult {
  vector<std::complex<double>> cycle, a;  // The cycle and its germ coefficients
  int bins = 1;                           // Local shells (1 otherwise)
  vector<double> shell;                   // Local: half-sides s_0 > s_1 > ... > s_bins
  vector<int64_t> counts;                 // bins × 66: [0, 64) escape octave, 64 certified in K, 65 undecided
  vector<double> weight;                  // Area per sample, by bin
  double secs = 0;
  int64_t iters = 0;

  // Area escaping in octave j (j = 64: certified, 65: undecided), summed over bins
  double area(const int j) const;
  // Its standard error
  double error(const int j) const;
};

TailResult julia_tail(const TailParams& p);

}  // namespace mandelbrot
