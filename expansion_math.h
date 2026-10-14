// Elementary functions for expansion arithmetic (sqrt, exp, log, sin/cos, atan2), host and device
//
// Each starts from the double result and refines: sqrt and log by Newton iteration (precision doubles per step), exp
// and sin/cos by range reduction (multiples of ln 2 or π/2, then halving s times) and a Taylor series, then
// squaring / double angles back; atan2 by rotating the point by the double angle and taking the (tiny) remaining
// angle's series.  Accurate to a few ulps of the expansion (tested against arb in expansion_math_test.cc).  The
// double overloads forward to libm, so code templated on the scalar runs unchanged in double.
#pragma once

#include "expansion_arith.h"
#include <cmath>
namespace mandelbrot {

namespace expansion_math_detail {
// π and ln 2 as nonoverlapping expansions (exact splits of the true values)
template<int n> __host__ __device__ static inline Expansion<n> pi_const() {
  constexpr double p[4] = {3.141592653589793, 1.2246467991473532e-16, -2.9947698097183397e-33, 1.1124542208633653e-49};
  Expansion<n> y;
  for (int i = 0; i < n; i++) y.x[i] = p[i];
  return y;
}
template<int n> __host__ __device__ static inline Expansion<n> ln2_const() {
  constexpr double p[4] = {0.6931471805599453, 2.3190468138462996e-17, 5.707708438416212e-34, -3.5824322106018114e-50};
  Expansion<n> y;
  for (int i = 0; i < n; i++) y.x[i] = p[i];
  return y;
}
// Taylor terms needed for |x| ≤ 2^-e at n doubles of precision
__host__ __device__ static constexpr int taylor_terms(const int n, const int e) { return (53 * n + 8) / e + 2; }
}  // namespace expansion_math_detail

template<class S> __host__ __device__ static inline S gl_pi();
template<> __host__ __device__ inline double gl_pi<double>() { return M_PI; }
template<> __host__ __device__ inline Expansion<2> gl_pi<Expansion<2>>() { return expansion_math_detail::pi_const<2>(); }
template<> __host__ __device__ inline Expansion<3> gl_pi<Expansion<3>>() { return expansion_math_detail::pi_const<3>(); }
template<> __host__ __device__ inline Expansion<4> gl_pi<Expansion<4>>() { return expansion_math_detail::pi_const<4>(); }

// double versions
__host__ __device__ static inline double gl_sqrt(const double x) { return ::sqrt(x); }
__host__ __device__ static inline double gl_exp(const double x) { return ::exp(x); }
__host__ __device__ static inline double gl_log(const double x) { return ::log(x); }
__host__ __device__ static inline void gl_sincos(const double t, double& s, double& c) { s = ::sin(t); c = ::cos(t); }
__host__ __device__ static inline double gl_atan2(const double y, const double x) { return ::atan2(y, x); }

template<int n> __host__ __device__ static inline Expansion<n> gl_sqrt(const Expansion<n> x) {
  const double d = double(x);
  if (!(d > 0)) return Expansion<n>(d > 0 ? d : 0.0);
  Expansion<n> y(::sqrt(d));
  // Newton y ← y + (x - y²)/(2y), precision doubling: 53 → 106 → 212
  for (int i = 1; i < 2 * n; i *= 2) y = y + (x - y * y) / twice(y);
  return y;
}

template<int n> __host__ __device__ static inline Expansion<n> gl_exp(const Expansion<n> x) {
  using namespace expansion_math_detail;
  const double d = double(x);
  if (!(d > -700 && d < 700)) return Expansion<n>(::exp(d));
  // x = k ln2 + r, |r| ≤ ln2/2; then r/2^s, Taylor, square s times
  const double k = ::nearbyint(d / 0.6931471805599453);
  const Expansion<n> r = x - ln2_const<n>() * k;
  const int s = 10;
  const Expansion<n> u = ldexp(r, -s);
  // carry e = exp(u) - 1, squaring as (1 + e)² - 1 = 2e + e², which keeps the relative precision (squaring the
  // value itself would lose s bits)
  const int terms = taylor_terms(n, s + 1);
  Expansion<n> e = u, t = u;
  for (int j = 2; j <= terms; j++) { t = t * u / Expansion<n>(int32_t(j)); e = e + t; }
  for (int j = 0; j < s; j++) e = twice(e) + e * e;
  return ldexp(e + Expansion<n>(1.0), int(k));
}

template<int n> __host__ __device__ static inline Expansion<n> gl_log(const Expansion<n> x) {
  const double d = double(x);
  if (!(d > 0)) return Expansion<n>(::log(d));
  // Newton on exp(y) = x: y ← y + x e^{-y} - 1 (precision doubles per step)
  Expansion<n> y(::log(d));
  for (int i = 1; i < 2 * n; i *= 2) y = y + x * gl_exp(-y) - Expansion<n>(1.0);
  return y;
}

template<int n> __host__ __device__ static inline void gl_sincos(const Expansion<n> t, Expansion<n>& sn, Expansion<n>& cs) {
  using namespace expansion_math_detail;
  // t = k π/2 + r, |r| ≤ π/4; sin and cos of r/2^s by Taylor, then double angles
  // reduce with one more part of π than the working precision, so results near zeros keep their relative accuracy
  constexpr double hp[4] = {1.5707963267948966, 6.123233995736766e-17, -1.4973849048591698e-33, 5.562271104316826e-50};
  const double k = ::nearbyint(double(t) / hp[0]);
  Expansion<n> r = t;
  for (int i = 0; i < (n < 4 ? n + 1 : 4); i++) r = r - Expansion<n>(hp[i]) * k;
  const int s = 8;
  const Expansion<n> u = ldexp(r, -s), u2 = u * u;
  const int terms = taylor_terms(n, s);
  // carry s = sin a and v = 1 - cos a (no cancellation): sin 2a = 2 s (1 - v), 1 - cos 2a = 2 s²
  Expansion<n> sv = u, vv = half(u2), ts = u, tc = half(u2);
  for (int j = 1; j <= terms; j++) {
    ts = -(ts * u2) / Expansion<n>(int32_t((2 * j) * (2 * j + 1)));
    tc = -(tc * u2) / Expansion<n>(int32_t((2 * j + 1) * (2 * j + 2)));
    sv = sv + ts;
    vv = vv + tc;
  }
  for (int j = 0; j < s; j++) {
    const Expansion<n> s2 = twice(sv - sv * vv), v2 = twice(sv * sv);
    sv = s2; vv = v2;
  }
  const Expansion<n> cv = Expansion<n>(1.0) - vv;
  const int q = ((int(k) % 4) + 4) % 4;
  if (q == 0) { sn = sv; cs = cv; }
  else if (q == 1) { sn = cv; cs = -sv; }
  else if (q == 2) { sn = -sv; cs = -cv; }
  else { sn = -cv; cs = sv; }
}

template<int n> __host__ __device__ static inline Expansion<n> gl_atan2(const Expansion<n> y, const Expansion<n> x) {
  // θ0 = the double angle; rotate (x, y) by -θ0 to (u, v) with v/u tiny, θ = θ0 + atan(v/u)
  const double t0 = ::atan2(double(y), double(x));
  Expansion<n> s0, c0;
  gl_sincos(Expansion<n>(t0), s0, c0);
  const Expansion<n> u = x * c0 + y * s0, v = y * c0 - x * s0;
  const Expansion<n> t = v / u, t2 = t * t;
  // atan t = t - t³/3 + t⁵/5 - ..., |t| ≲ 2^-52: n terms suffice
  Expansion<n> sum = t, p = t;
  for (int j = 1; j < n; j++) { p = -(p * t2); sum = sum + p / Expansion<n>(int32_t(2 * j + 1)); }
  return Expansion<n>(t0) + sum;
}

}  // namespace mandelbrot
