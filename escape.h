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

// Classification with distance estimates (see EscapeDE in orbit.h)
EscapeDE escape_de(const double x, const double y, const int64_t max_iter);

// Whether g_M(c) < 2^-k, decided from an Escape computed with max_iter ≥ k + 8
__host__ __device__ static inline bool below(const Escape& e, const int k) {
  return e.steps < 0 || e.log2g < -k;
}

}  // namespace mandelbrot
