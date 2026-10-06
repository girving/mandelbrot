// Orbits over Expansion<2> (double-double, about 106 bits), as a precision reference for double, and over
// Expansion<3> (triple-double, about 159 bits) for high precision dynamics
//
// Orbit<Expansion<2>> runs double's algorithm with double's tolerances, so paired differences against double
// (escape_tree --prec comparedd) measure double's rounding alone.
#pragma once

#include "expansion_arith.h"
#include "orbit.h"
namespace mandelbrot {

// Same algorithm and tolerances as double; the parameter c and constants are exact doubles, which mixed
// expansion-double arithmetic handles more cheaply than promoted expansions
template<> struct OrbitTol<Expansion<2>> : public OrbitTol<double> {};
template<> struct OrbitParam<Expansion<2>> { typedef double type; };

// Exact doubling, componentwise
__host__ __device__ inline Expansion<2> orbit_twice(const Expansion<2> a) { return twice(a); }

// The step's fused multiply-adds: exact products and sums to about 2^-104
__host__ __device__ inline Expansion<2> orbit_fma(const Expansion<2> a, const Expansion<2> b, const Expansion<2> c) {
  return a * b + c;
}
__host__ __device__ inline Expansion<2> orbit_fma(const Expansion<2> a, const Expansion<2> b, const double c) {
  return a * b + c;
}

// The whole step z → z^2 + c, fused: 41 FP64 operations against 84 for the generic expansion arithmetic above
// (double's step is 6).  With z = (a + al) + i (b + bl) and c = x + i y exact doubles,
//   Re z' = a^2 - b^2 + x + 2 (a al - b bl),   Im z' = 2ab + y + 2 (a bl + al b),
// where the high parts' squares and products are exact (two_prod by fma), their sums with each other and c exact
// (two_sum), and the cross terms and errors summed in double; the dropped al^2, bl^2, al bl terms are ~2^-106
// |z|^2.  Per-step error ≤ 2 · 2^-104 (|z|^2 + |c|) (tested).  zy2 is not needed (carried for the generic
// signature), and r2 = |z'|^2 is computed from the high parts only, since it serves the escape test and the
// escape radius.
template<> __host__ __device__ inline void orbit_step<Expansion<2>>(Expansion<2>& zx, Expansion<2>& zy,
                                                                    Expansion<2>& zy2, Expansion<2>& r2,
                                                                    const double x, const double y) {
  const double a = zx.x[0], al = zx.x[1], b = zy.x[0], bl = zy.x[1];
  // Re: a^2 - b^2 + x exactly as s2 + (f + g + e1 - e2), plus the cross terms
  const double p1 = a * a, e1 = fma(a, a, -p1), p2 = b * b, e2 = fma(b, b, -p2);
  const double s = p1 - p2, sb = s - p1, f = (p1 - (s - sb)) + (-p2 - sb);
  const double s2 = s + x, xb = s2 - s, g = (s - (s2 - xb)) + (x - xb);
  const double lx = (f + g) + (e1 - e2) + 2 * fma(a, al, -(b * bl));
  // Im: 2ab + y exactly as u + (eq + h), plus the cross terms
  const double ta = a + a, q = ta * b, eq = fma(ta, b, -q);
  const double u = q + y, yb = u - q, h = (q - (u - yb)) + (y - yb);
  const double ly = (h + eq) + 2 * fma(a, bl, al * b);
  // Renormalize (fast two_sum: the high parts dominate unless they cancel, where the error is still ~2^-106 |z|^2)
  const double xh = s2 + lx, yh = u + ly;
  zx = Expansion<2>(xh, lx - (xh - s2), nonoverlap);
  zy = Expansion<2>(yh, ly - (yh - u), nonoverlap);
  zy2 = zy;
  r2 = Expansion<2>(fma(xh, xh, yh * yh));
}

// Error-free transforms: s + e = a + b exactly (orbit_two_sum), and also when |a| ≥ |b| (orbit_fast_two_sum)
__host__ __device__ inline void orbit_two_sum(const double a, const double b, double& s, double& e) {
  s = a + b;
  const double t = s - a;
  e = (a - (s - t)) + (b - t);
}
__host__ __device__ inline void orbit_fast_two_sum(const double a, const double b, double& s, double& e) {
  s = a + b;
  e = b - (s - a);
}

// Expansion<3> mirrors Expansion<2>: double's tolerances, exact double parameters and constants
template<> struct OrbitTol<Expansion<3>> : public OrbitTol<double> {};
template<> struct OrbitParam<Expansion<3>> { typedef double type; };
__host__ __device__ inline Expansion<3> orbit_twice(const Expansion<3> a) { return twice(a); }
__host__ __device__ inline Expansion<3> orbit_fma(const Expansion<3> a, const Expansion<3> b, const Expansion<3> c) {
  return a * b + c;
}
__host__ __device__ inline Expansion<3> orbit_fma(const Expansion<3> a, const Expansion<3> b, const double c) {
  return a * b + c;
}

// The whole step z → z^2 + c over Expansion<3>, fused.  With z = (a + a1 + a2) + i (b + b1 + b2) and c = x + i y
// exact doubles, sort the terms by size relative to Z = |z|^2 + |c|:
//   order 0:  a^2 - b^2 + x,  2ab + y                              exact (two_prod by fma, two_sum)
//   order 1:  2 (a a1 - b b1),  2 (a b1 + a1 b), and the errors      exact products, summed with two_sum
//   order 2:  a1^2 - b1^2 + 2 (a a2 - b b2),  2 (a b2 + a2 b + a1 b1), and the errors   in double
// The dropped terms (a1 a2, a2^2, ...) are ≤ 2^-157 Z, and the order 2 roundings are ~2^-159 Z each.  Per-step
// error ≤ 2^-156 Z (tested: worst 0.6 · 2^-156).  Renormalization is a two_sum of the order 0 and 1 parts
// (which can cancel), then a fast_two_sum with order 2, whose failures cost only ~2^-159 Z.  The result is
// nonoverlapping unless a component cancels below ~2^-50 Z, where it can overlap but stays accurate.  As for
// Expansion<2>, zy2 is unused and r2 = |z'|^2 comes from the high parts.
template<> __host__ __device__ inline void orbit_step<Expansion<3>>(Expansion<3>& zx, Expansion<3>& zy,
                                                                    Expansion<3>& zy2, Expansion<3>& r2,
                                                                    const double x, const double y) {
  const double a = zx.x[0], a1 = zx.x[1], a2 = zx.x[2], b = zy.x[0], b1 = zy.x[1], b2 = zy.x[2];
  const double ta = a + a, tb = b + b;
  double s, f, s0, g, u1, v1, u2, v2, u3, v3, w, v4, t1, v5, xh, r, xm, xl;
  // Re, order 0: a^2 - b^2 + x = s0 + (f + g + e1 - e2) exactly
  const double p1 = a * a, e1 = fma(a, a, -p1), p2 = b * b, e2 = fma(b, b, -p2);
  orbit_two_sum(p1, -p2, s, f);
  orbit_two_sum(s, x, s0, g);
  // Order 1: 2 (a a1 - b b1) = (m - n) + (me - ne) exactly; with the errors above, six order 1 terms summed
  // exactly to t1 + (v1 + ... + v5)
  const double m = ta * a1, me = fma(ta, a1, -m), n = tb * b1, ne = fma(tb, b1, -n);
  orbit_two_sum(e1, -e2, u1, v1);
  orbit_two_sum(m, -n, u2, v2);
  orbit_two_sum(f, g, u3, v3);
  orbit_two_sum(u1, u2, w, v4);
  orbit_two_sum(w, u3, t1, v5);
  // Order 2
  const double l2 = ((v1 + v3) + (v4 + v5)) + ((v2 + (me - ne)) + fma(ta, a2, -(tb * b2))) + fma(a1, a1, -(b1 * b1));
  orbit_two_sum(s0, t1, xh, r);
  orbit_fast_two_sum(r, l2, xm, xl);

  double u0, h, c1, d1, c2, d2, ty1, d3, yh, ry, ym, yl;
  // Im, order 0: 2ab + y = u0 + (h + eq) exactly
  const double q = ta * b, eq = fma(ta, b, -q);
  orbit_two_sum(q, y, u0, h);
  // Order 1: 2 (a b1 + a1 b) = (k1 + k2) + (k1e + k2e) exactly; four order 1 terms summed exactly
  const double k1 = ta * b1, k1e = fma(ta, b1, -k1), k2 = a1 * tb, k2e = fma(a1, tb, -k2);
  orbit_two_sum(h, eq, c1, d1);
  orbit_two_sum(k1, k2, c2, d2);
  orbit_two_sum(c1, c2, ty1, d3);
  // Order 2
  const double l2y = ((d1 + d2) + (d3 + (k1e + k2e))) + fma(ta, b2, fma(tb, a2, a1 * (b1 + b1)));
  orbit_two_sum(u0, ty1, yh, ry);
  orbit_fast_two_sum(ry, l2y, ym, yl);

  zx = Expansion<3>(xh, xm, xl, nonoverlap);
  zy = Expansion<3>(yh, ym, yl, nonoverlap);
  zy2 = zy;
  r2 = Expansion<3>(fma(xh, xh, yh * yh));
}

}  // namespace mandelbrot
