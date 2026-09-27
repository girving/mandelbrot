// Escape-time classification of Mandelbrot parameters, with the Green's function of M
#pragma once

#include <cmath>
#include <cstdint>
namespace mandelbrot {

// Whether c is in the main cardioid or the period 2 disk (both inside M)
bool in_cardioid_or_disk(const double x, const double y);

// Result of iterating z → z^2 + c from z = c for at most max_iter steps
struct Escape {
  int64_t steps;  // Steps until |z| > 2^32, or -1 if it never escaped (or an attracting cycle was found)
  double log2g;   // log2 of the Green's function g_M(c) = lim 2^-n log|z_n| if escaped, else -inf
  int period = 0; // Minimal period of the attracting cycle if one was found with period ≤ 32, else 0
  int64_t iters = 0;  // Iterations performed
  double g() const { return std::exp2(log2g); }
};

// Iterate with an attracting-cycle check.  If the orbit escapes after m steps with z_m = f^(m-1)(c),
// g_M(c) = log|z_m| / 2^(m-1) up to a relative error of about |c| / |z_m|^2.  We return log2 of that,
// since 2^(m-1) overflows for long orbits.
Escape escape(const double x, const double y, const int64_t max_iter);

// Classification with distance estimates, for certifying whole cells.  Exterior: log2 of the Green's function
// and a lower bound on dist(c, M) from Koebe, |z_n| log|z_n| / (2 |dz_n/dc|).  Interior (attracting cycle
// found by Newton): the Koebe lower bound (1 - |A|^2) / (4 |D + C B / (1 - A)|) on the distance to the
// component's boundary, where A, B, C, D are derivatives of f^p at the periodic point.  Zero if unknown.
struct EscapeDE {
  Escape e;
  double dist = 0;  // Exterior: lower bound on dist(c, M).  Interior: lower bound on distance to ∂(component).
};
EscapeDE escape_de(const double x, const double y, const int64_t max_iter);

// Whether g_M(c) < 2^-k, decided from an Escape computed with max_iter ≥ k + 8
static inline bool below(const Escape& e, const int k) {
  return e.steps < 0 || e.log2g < -k;
}

}  // namespace mandelbrot
