// Escape-time classification of Mandelbrot parameters

#include "escape.h"
namespace mandelbrot {

bool in_cardioid_or_disk(const double x, const double y) {
  // Period 2 disk: |c + 1| ≤ 1/4
  const double y2 = y * y;
  if ((x + 1) * (x + 1) + y2 <= 1.0 / 16) return true;
  // Main cardioid: q (q + x - 1/4) ≤ y^2 / 4 with q = (x - 1/4)^2 + y^2
  const double q = (x - 0.25) * (x - 0.25) + y2;
  return q * (q + (x - 0.25)) <= 0.25 * y2;
}

namespace {

// Newton's method for an attracting p-cycle of z → z^2 + c near w.  Returns true if Newton converges to a
// periodic point whose multiplier |(f^p)'(w)| < 1, which certifies that c is in a hyperbolic component.
bool attracting_cycle(const double x, const double y, double wx, double wy, const int p) {
  for (int it = 0; it < 30; it++) {
    // F(w) = f^p(w) - w, F'(w) = (f^p)'(w) - 1
    double zx = wx, zy = wy, dx = 1, dy = 0;
    for (int k = 0; k < p; k++) {
      const double ndx = 2 * (zx * dx - zy * dy), ndy = 2 * (zx * dy + zy * dx);
      dx = ndx; dy = ndy;
      const double t = zx * zx - zy * zy + x;
      zy = 2 * zx * zy + y;
      zx = t;
      if (zx * zx + zy * zy > 16) return false;
    }
    const double fx = zx - wx, fy = zy - wy, gx = dx - 1, gy = dy;
    const double den = gx * gx + gy * gy;
    if (!(den > 0)) return false;
    const double sx = (fx * gx + fy * gy) / den, sy = (fy * gx - fx * gy) / den;
    wx -= sx; wy -= sy;
    if (sx * sx + sy * sy < 1e-28 * (1 + wx * wx + wy * wy)) {
      // Converged: the multiplier at the periodic point decides
      double mx = 1, my = 0, zx2 = wx, zy2 = wy;
      for (int k = 0; k < p; k++) {
        const double nmx = 2 * (zx2 * mx - zy2 * my), nmy = 2 * (zx2 * my + zy2 * mx);
        mx = nmx; my = nmy;
        const double t = zx2 * zx2 - zy2 * zy2 + x;
        zy2 = 2 * zx2 * zy2 + y;
        zx2 = t;
      }
      return mx * mx + my * my < 1 - 1e-9;
    }
  }
  return false;
}

}  // namespace

Escape escape(const double x, const double y, const int64_t max_iter) {
  const double inside = -INFINITY;
  if (in_cardioid_or_disk(x, y)) {
    const bool disk = (x + 1) * (x + 1) + y * y <= 1.0 / 16;
    return {-1, inside, disk ? 2 : 1, 0};
  }
  // Iterate z_1 = c, z_{n+1} = z_n^2 + c, escaping at |z| > 2^32 so that log|z| is accurate
  const double R2 = std::ldexp(1.0, 64);
  double zx = x, zy = y;
  // Attracting-cycle detection (Brent): compare against a checkpoint refreshed at powers of two
  double cx = zx, cy = zy;
  int64_t next_check = 16;
  // Atom domains: the step where |z_n| reaches a new minimum is a candidate period
  double min_r2 = zx * zx + zy * zy;
  int64_t candidate = 1, next_newton = 64;
  for (int64_t n = 1; n <= max_iter; n++) {
    const double r2 = zx * zx + zy * zy;
    if (r2 > R2) return {n, std::log2(0.5 * std::log(r2)) - double(n - 1), 0, n};
    const double t = zx * zx - zy * zy + x;
    zy = 2 * zx * zy + y;
    zx = t;
    const double r2n = zx * zx + zy * zy;
    if (r2n < min_r2) { min_r2 = r2n; candidate = n + 1; }  // (zx, zy) is now z_{n+1}
    if (n == next_newton) {
      next_newton *= 2;
      if (candidate <= 4096 && attracting_cycle(x, y, zx, zy, int(candidate))) {
        // Report the minimal period if it is small
        int period = 0;
        if (candidate <= 32)
          for (int q = 1; q <= int(candidate); q++)
            if (int(candidate) % q == 0 && attracting_cycle(x, y, zx, zy, q)) { period = q; break; }
        return {-1, inside, period, n};
      }
    }
    const double dx = zx - cx, dy = zy - cy;
    if (dx * dx + dy * dy < 1e-26) {
      // Converged to an attracting cycle: find its minimal period, if small
      double wx = zx, wy = zy;
      for (int p = 1; p <= 32; p++) {
        const double t2 = wx * wx - wy * wy + x;
        wy = 2 * wx * wy + y;
        wx = t2;
        const double ex = wx - zx, ey = wy - zy;
        if (ex * ex + ey * ey < 1e-20) return {-1, inside, p, n};
      }
      return {-1, inside, 0, n};
    }
    if (n == next_check) { cx = zx; cy = zy; next_check *= 2; }
  }
  return {-1, inside, 0, max_iter};
}

}  // namespace mandelbrot
