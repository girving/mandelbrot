// Orbits over Expansion<2> (double-double, about 106 bits), as a precision reference for double
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

}  // namespace mandelbrot
