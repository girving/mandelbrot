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

Escape escape(const double x, const double y, const int64_t max_iter) {
  const double inside = -INFINITY;
  if (in_cardioid_or_disk(x, y)) {
    const bool disk = (x + 1) * (x + 1) + y * y <= 1.0 / 16;
    return {-1, inside, disk ? 2 : 1};
  }
  // Iterate z_1 = c, z_{n+1} = z_n^2 + c, escaping at |z| > 2^32 so that log|z| is accurate
  const double R2 = std::ldexp(1.0, 64);
  double zx = x, zy = y;
  // Attracting-cycle detection (Brent): compare against a checkpoint refreshed at powers of two
  double cx = zx, cy = zy;
  int64_t next_check = 16;
  for (int64_t n = 1; n <= max_iter; n++) {
    const double r2 = zx * zx + zy * zy;
    if (r2 > R2) return {n, std::log2(0.5 * std::log(r2)) - double(n - 1)};
    const double t = zx * zx - zy * zy + x;
    zy = 2 * zx * zy + y;
    zx = t;
    const double dx = zx - cx, dy = zy - cy;
    if (dx * dx + dy * dy < 1e-26) {
      // Converged to an attracting cycle: find its minimal period, if small
      double wx = zx, wy = zy;
      for (int p = 1; p <= 32; p++) {
        const double t2 = wx * wx - wy * wy + x;
        wy = 2 * wx * wy + y;
        wx = t2;
        const double ex = wx - zx, ey = wy - zy;
        if (ex * ex + ey * ey < 1e-20) return {-1, inside, p};
      }
      return {-1, inside, 0};
    }
    if (n == next_check) { cx = zx; cy = zy; next_check *= 2; }
  }
  return {-1, inside};
}

}  // namespace mandelbrot
