// Escape-time classification of Mandelbrot parameters, with the Green's function of M
#pragma once

#include "orbit.h"
#include <cmath>
#include <cstdint>
namespace mandelbrot {

// Iterate with an attracting-cycle check.  If the orbit escapes after m steps with z_m = f^(m-1)(c),
// g_M(c) = log|z_m| / 2^(m-1) up to a relative error of about |c| / |z_m|^2.  We return log2 of that,
// since 2^(m-1) overflows for long orbits.
Escape escape(const double x, const double y, const int64_t max_iter);

// Classification with distance estimates, for certifying whole cells.  Exterior: log2 of the Green's function
// and the Koebe lower bound on dist(c, M), (1 - e^-g) / (4 |∇g|) ≈ |z_n| log|z_n| / (4 |dz_n/dc|).  Interior (attracting cycle
// found by Newton): the Koebe lower bound (1 - |A|^2) / (4 |D + C B / (1 - A)|) on the distance to the
// component's boundary, where A, B, C, D are derivatives of f^p at the periodic point.  Zero if unknown.
struct EscapeDE {
  Escape e;
  double dist = 0;  // Exterior: lower bound on dist(c, M).  Interior: lower bound on distance to ∂(component).
};
EscapeDE escape_de(const double x, const double y, const int64_t max_iter);

// Whether g_M(c) < 2^-k, decided from an Escape computed with max_iter ≥ k + 8
__host__ __device__ static inline bool below(const Escape& e, const int k) {
  return e.steps < 0 || e.log2g < -k;
}

}  // namespace mandelbrot
