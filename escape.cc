// Escape-time classification of Mandelbrot parameters

#include "escape.h"
#include <complex>
namespace mandelbrot {

namespace {

// Interior distance lower bound at an attracting p-cycle near w (p must be the minimal period, since the
// Koebe bound needs the multiplier map to be univalent), or 0 if Newton does not certify one
double interior_distance_exact(const double x, const double y, double wx, double wy, const int p) {
  typedef std::complex<double> C;
  const C c(x, y);
  C w(wx, wy);
  for (int it = 0; it < 30; it++) {
    C z = w, dz = 1;
    for (int k = 0; k < p; k++) { dz = 2.0 * z * dz; z = z * z + c; }
    const C step = (z - w) / (dz - 1.0);
    w -= step;
    if (std::norm(step) < 1e-28 * (1 + std::norm(w))) {
      // Derivatives of F = f^p at the periodic point w: A = F_z, B = F_c, Cz = F_zz, D = F_zc
      C z2 = w, A = 1, B = 0, Cz = 0, D = 0;
      for (int k = 0; k < p; k++) {
        const C nA = 2.0 * z2 * A, nB = 2.0 * z2 * B + 1.0;
        const C nC = 2.0 * (A * A + z2 * Cz), nD = 2.0 * (A * B + z2 * D);
        A = nA; B = nB; Cz = nC; D = nD;
        z2 = z2 * z2 + c;
      }
      const double a2 = std::norm(A);
      if (!(a2 < 1 - 1e-9)) return 0;
      return (1 - a2) / (4 * std::abs(D + Cz * B / (1.0 - A)));
    }
  }
  return 0;
}

// Interior distance at the minimal period dividing p for which Newton finds an attracting cycle
double interior_distance(const double x, const double y, double wx, double wy, const int p) {
  if (!attracting_cycle(x, y, wx, wy, p)) return 0;
  for (int q = 1; q <= p; q++)
    if (p % q == 0 && attracting_cycle(x, y, wx, wy, q)) return interior_distance_exact(x, y, wx, wy, q);
  return 0;
}

}  // namespace

EscapeDE escape_de(const double x, const double y, const int64_t max_iter) {
  EscapeDE r;
  if (in_cardioid_or_disk(x, y)) {
    // Distance from the cardioid/disk boundary is at least the Newton interior bound; recompute it
    const bool disk = (x + 1) * (x + 1) + y * y <= 1.0 / 16;
    r.e = {-1, -INFINITY, disk ? 2 : 1, 0};
    // Start Newton from an orbit point near the attracting cycle (not a fixed guess, which can converge to
    // a repelling cycle instead)
    double zx = x, zy = y;
    for (int n = 0; n < 256; n++) {
      const double t = zx * zx - zy * zy + x;
      zy = 2 * zx * zy + y;
      zx = t;
    }
    r.dist = interior_distance(x, y, zx, zy, disk ? 2 : 1);
    return r;
  }
  // Iterate z and dz/dc together.  dz/dc grows like 2^n on escaping orbits, so rescale it and carry a
  // binary exponent to avoid overflow.
  double zx = x, zy = y, dx = 1, dy = 0;
  int64_t dexp = 0;
  const double R2 = std::ldexp(1.0, 64);
  double min_r2 = zx * zx + zy * zy;
  int64_t candidate = 1, next_newton = 64;
  double cx = zx, cy = zy;  // Brent checkpoint
  int64_t check_n = 1, next_check = 16;
  for (int64_t n = 1; n <= max_iter; n++) {
    const double r2 = zx * zx + zy * zy;
    if (r2 > R2) {
      const double lz = 0.5 * std::log(r2);
      r.e = {n, std::log2(lz) - double(n - 1), 0, n};
      // Koebe: dist(c, M) ≥ (1 - e^-g) / (4 |∇g|), with g = log|z_n| / 2^(n-1) and
      // |∇g| = |dz_n/dc| / (|z_n| 2^(n-1)).  So dist ≥ (1 - e^-g) 2^(n-1) |z_n| / (4 |dz_n/dc|), where
      // (1 - e^-g) 2^(n-1) = lz for small g.
      const double g = std::exp2(std::log2(lz) - double(n - 1));
      const double scale = g > 1e-8 ? -std::expm1(-g) / g : 1 - g / 2;  // (1 - e^-g) / g
      r.dist = std::ldexp(scale * std::sqrt(r2) * lz / (4 * std::hypot(dx, dy)), int(-std::min<int64_t>(dexp, 100000)));
      return r;
    }
    const double ndx = 2 * (zx * dx - zy * dy) + std::ldexp(1.0, int(-std::min<int64_t>(dexp, 2000))),
                 ndy = 2 * (zx * dy + zy * dx);
    dx = ndx; dy = ndy;
    if (dx * dx + dy * dy > 1e200) { dx = std::ldexp(dx, -256); dy = std::ldexp(dy, -256); dexp += 256; }
    const double t = zx * zx - zy * zy + x;
    zy = 2 * zx * zy + y;
    zx = t;
    const double r2n = zx * zx + zy * zy;
    if (r2n < min_r2) { min_r2 = r2n; candidate = n + 1; }
    if (n == next_newton) {
      next_newton *= 2;
      if (candidate <= 4096) {
        const double b = interior_distance(x, y, zx, zy, int(candidate));
        if (b > 0) { r.e = {-1, -INFINITY, 0, n}; r.dist = b; return r; }
      }
    }
    // Brent fallback: converged to a cycle; recover its period by iterating until the orbit returns
    const double ex = zx - cx, ey = zy - cy;
    if (ex * ex + ey * ey < 1e-26) {
      const int64_t lag = n - check_n;
      double wx = zx, wy = zy;
      for (int64_t q = 1; q <= std::min<int64_t>(lag, 1 << 16); q++) {
        const double t2 = wx * wx - wy * wy + x;
        wy = 2 * wx * wy + y;
        wx = t2;
        const double fx = wx - zx, fy = wy - zy;
        if (fx * fx + fy * fy < 1e-20) {
          r.e = {-1, -INFINITY, 0, n};
          r.dist = interior_distance(x, y, zx, zy, int(q));
          return r;
        }
      }
      r.e = {-1, -INFINITY, 0, n};
      return r;
    }
    if (n == next_check) { cx = zx; cy = zy; check_n = n; next_check *= 2; }
  }
  r.e = {-1, -INFINITY, 0, max_iter};
  return r;
}

Escape escape(const double x, const double y, const int64_t max_iter) {
  Orbit<double> o;
  if (!o.start(x, y))
    o.run(max_iter, max_iter);
  return o.e;
}

}  // namespace mandelbrot
