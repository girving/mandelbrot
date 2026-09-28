// Resumable escape-time orbits, for batched iteration on CPU lanes and GPU threads
//
// Orbit<T> is escape() as a state machine: start() at a parameter, then call step() until it reports the
// result.  One step is one iteration of z → z^2 + c plus the (rare) attracting-cycle checks, so CPU code can
// interleave many independent orbits for instruction-level parallelism, and GPU threads can refill finished
// orbits with new samples without waiting for the rest of their warp.  T is double or float; tolerances
// scale with the precision.
#pragma once

#include "complex.h"
#include "cutil.h"
#include <cmath>
#include <cstdint>
namespace mandelbrot {

// Whether c is in the main cardioid or the period 2 disk (both inside M)
__host__ __device__ static inline bool in_cardioid_or_disk(const double x, const double y) {
  // Period 2 disk: |c + 1| ≤ 1/4
  const double y2 = y * y;
  if ((x + 1) * (x + 1) + y2 <= 1.0 / 16) return true;
  // Main cardioid: q (q + x - 1/4) ≤ y^2 / 4 with q = (x - 1/4)^2 + y^2
  const double q = (x - 0.25) * (x - 0.25) + y2;
  return q * (q + (x - 0.25)) <= 0.25 * y2;
}

// Result of iterating z → z^2 + c from z = c for at most max_iter steps
struct Escape {
  int64_t steps;  // Steps until |z| > 2^32, or -1 if it never escaped (or an attracting cycle was found)
  double log2g;   // log2 of the Green's function g_M(c) = lim 2^-n log|z_n| if escaped, else -inf
  int period = 0; // Minimal period of the attracting cycle if one was found with period ≤ 32, else 0
  int64_t iters = 0;  // Iterations performed
  double r2 = 0;      // |z_steps|^2 if escaped
  double g() const { return std::exp2(log2g); }
};

// Whether an escaped orbit has g_M(c) < 2^-k, from its escape step and |z|^2 alone, without logarithms.
// log2 g = log2(log(r2) / 2) - (steps - 1), and r2 ∈ (2^64, 2^128] at escape, so log(r2) / 2 ∈ (22.1, 44.4]:
// with d = steps - 1 - k, g < 2^-k iff log(r2) / 2 < 2^d, which holds for d ≥ 6, fails for d ≤ 4, and for
// d = 5 means r2 < e^64.
__host__ __device__ static inline bool escaped_below(const int64_t steps, const double r2, const int k) {
  const int64_t d = steps - 1 - k;
  return d >= 6 || (d == 5 && r2 < 6.235149080811617e27);
}

// Rare, heavy paths (Newton certificates, distance bounds) stay out of line, so that the hot iteration loops
// that call them keep few registers on GPUs
#define ORBIT_COLD __attribute__((noinline))

// Tolerances by precision (squared distances, relative where noted)
template<class T> struct OrbitTol;
template<> struct OrbitTol<double> {
  static constexpr double cycle = 1e-26;      // Orbit returned to the Brent checkpoint: |Δz|^2
  static constexpr double period = 1e-20;     // Period recovery: |f^q(z) - z|^2
  static constexpr double newton = 1e-28;     // Newton convergence: |step|^2 relative to 1 + |w|^2
  static constexpr double multiplier = 1e-9;  // Attracting if |λ|^2 < 1 - multiplier
};
template<> struct OrbitTol<float> {
  static constexpr double cycle = 1e-12;
  static constexpr double period = 1e-10;
  static constexpr double newton = 1e-12;
  static constexpr double multiplier = 1e-5;
};

// Newton's method for an attracting p-cycle of z → z^2 + c near w.  Returns true if Newton converges to a
// periodic point whose multiplier |(f^p)'(w)| < 1, which certifies that c is in a hyperbolic component.
template<class T> ORBIT_COLD __host__ __device__ bool attracting_cycle(const T x, const T y, T wx, T wy, const int p,
                                                         const int iters = 30) {
  typedef OrbitTol<T> Tol;
  for (int it = 0; it < iters; it++) {
    // F(w) = f^p(w) - w, F'(w) = (f^p)'(w) - 1
    T zx = wx, zy = wy, dx = 1, dy = 0;
    for (int k = 0; k < p; k++) {
      const T ndx = 2 * (zx * dx - zy * dy), ndy = 2 * (zx * dy + zy * dx);
      dx = ndx; dy = ndy;
      const T t = zx * zx - zy * zy + x;
      zy = 2 * zx * zy + y;
      zx = t;
      if (zx * zx + zy * zy > 16) return false;
    }
    const T fx = zx - wx, fy = zy - wy, gx = dx - 1, gy = dy;
    const T den = gx * gx + gy * gy;
    if (!(den > 0)) return false;
    const T sx = (fx * gx + fy * gy) / den, sy = (fy * gx - fx * gy) / den;
    wx -= sx; wy -= sy;
    if (sx * sx + sy * sy < T(Tol::newton) * (1 + wx * wx + wy * wy)) {
      // Converged: the multiplier at the periodic point decides
      T mx = 1, my = 0, zx2 = wx, zy2 = wy;
      for (int k = 0; k < p; k++) {
        const T nmx = 2 * (zx2 * mx - zy2 * my), nmy = 2 * (zx2 * my + zy2 * mx);
        mx = nmx; my = nmy;
        const T t = zx2 * zx2 - zy2 * zy2 + x;
        zy2 = 2 * zx2 * zy2 + y;
        zx2 = t;
      }
      return mx * mx + my * my < T(1 - Tol::multiplier);
    }
  }
  return false;
}

// Orbit and OrbitDE keep step counts in 32 bits to save GPU registers, so max_iter must be below 2^30
template<class T> struct Orbit {
  T x, y, zx, zy;
  T cx, cy;        // Brent checkpoint, refreshed at powers of two
  T min_r2;        // Atom domains: the step where |z_n| reaches a new minimum is a candidate period
  int32_t n, next_check, candidate, next_newton;  // 32 bits to save registers: max_iter < 2^31
  int max_period;  // Largest atom-domain period candidate that Newton tries
  int newton_iters;  // Newton iterations per certificate attempt
  bool logs;       // Compute e.log2g on escape (otherwise only e.steps and e.r2, for escaped_below)
  Escape e;        // Result, once done

  // Start at c = x + iy, with Newton certificate attempts at step first_newton and then at each doubling, for
  // period candidates up to max_period.  Returns true if already decided (the cardioid or period 2 disk).
  __host__ __device__ bool start(const double x_, const double y_, const int64_t first_newton = 64,
                                 const int max_period = 4096, const bool logs = true, const int newton_iters = 30) {
    this->max_period = max_period;
    this->newton_iters = newton_iters;
    this->logs = logs;
    x = T(x_); y = T(y_);
    zx = x; zy = y; cx = x; cy = y;
    min_r2 = zx * zx + zy * zy;
    n = 1; next_check = 16; candidate = 1;
    next_newton = int32_t(first_newton < (int64_t(1) << 30) ? first_newton : int64_t(1) << 30);  // ≥ 2^30: never
    if (in_cardioid_or_disk(x_, y_)) {
      const bool disk = (x_ + 1) * (x_ + 1) + y_ * y_ <= 1.0 / 16;
      e = {-1, -INFINITY, disk ? 2 : 1, 0};
      return true;
    }
    return false;
  }

  // Iterate at most `budget` steps.  Returns true when done, with the result in e.  State lives in locals
  // during the loop so that it stays in registers.  Each step does the escape test and z → z^2 + c, carrying
  // the squares, and tracks the atom-domain minimum while n < max_period (later candidates are never used).
  // The rarer checks (Newton, Brent, checkpoint) run at steps that are multiples of 8: they cost as much as
  // the iteration itself on GPUs, and tying them to absolute step numbers keeps results independent of how
  // the orbit is split into bursts.
  __host__ __device__ bool run(const int64_t max_iter, const int64_t budget) {
    typedef OrbitTol<T> Tol;
    T zx = this->zx, zy = this->zy, cx = this->cx, cy = this->cy, min_r2 = this->min_r2;
    T zx2 = zx * zx, zy2 = zy * zy, r2 = zx2 + zy2;
    int32_t n = this->n, next_check = this->next_check, candidate = this->candidate;
    const int32_t end = int32_t(n + budget < max_iter + 1 ? n + budget : max_iter + 1);
    bool done = true;
    while (n < end) {
      const int32_t block = (n | 7) + 1 < end ? (n | 7) + 1 : end;  // Up to the next multiple of 8
      for (; n < block; n++) {
        // Escape at |z| > 2^32 so that log|z| is accurate
        if (r2 > T(18446744073709551616.0)) {
          e = {n, logs ? std::log2(0.5 * std::log(double(r2))) - double(n - 1) : NAN, 0, n, double(r2)};
          goto finish;
        }
        const T xy = zx * zy;
        zx = zx2 - zy2 + x;
        zy = xy + xy + y;
        zx2 = zx * zx; zy2 = zy * zy; r2 = zx2 + zy2;
        if (n < max_period && r2 < min_r2) { min_r2 = r2; candidate = n + 1; }  // (zx, zy) is now z_{n+1}
      }
      if (n & 7) continue;  // Partial block at the end of the burst
      // (zx, zy) is z_n, n a multiple of 8
      if (n > next_newton) [[unlikely]] {
        while (next_newton < n) next_newton *= 2;
        if (candidate <= max_period && attracting_cycle(x, y, zx, zy, int(candidate), newton_iters)) {
          // Report the minimal period if it is small
          int period = 0;
          if (candidate <= 32)
            for (int q = 1; q <= int(candidate); q++)
              if (int(candidate) % q == 0 && attracting_cycle(x, y, zx, zy, q)) { period = q; break; }
          e = {-1, -INFINITY, period, n};
          goto finish;
        }
      }
      {
        const T dx = zx - cx, dy = zy - cy;
        if (dx * dx + dy * dy < T(Tol::cycle)) [[unlikely]] {
          // Converged to an attracting cycle: find its minimal period, if small
          T wx = zx, wy = zy;
          int period = 0;
          for (int p = 1; p <= 32; p++) {
            const T t2 = wx * wx - wy * wy + x;
            wy = 2 * wx * wy + y;
            wx = t2;
            const T ex = wx - zx, ey = wy - zy;
            if (ex * ex + ey * ey < T(Tol::period)) { period = p; break; }
          }
          e = {-1, -INFINITY, period, n};
          goto finish;
        }
      }
      if (n > next_check) { cx = zx; cy = zy; next_check *= 2; }
    }
    if (n > max_iter) e = {-1, -INFINITY, 0, max_iter};
    else done = false;
  finish:
    this->zx = zx; this->zy = zy; this->cx = cx; this->cy = cy; this->min_r2 = min_r2;
    this->n = n; this->next_check = next_check; this->candidate = candidate;
    return done;
  }
};

// Distance-estimating orbits (escape_de): z and dz/dc, with Koebe distance bounds on exit

// Interior distance lower bound at an attracting p-cycle near w (p must be the minimal period, since the
// Koebe bound needs the multiplier map to be univalent), or 0 if Newton does not certify one
ORBIT_COLD __host__ __device__ static double interior_distance_exact(const double x, const double y, const double wx,
                                                                 const double wy, const int p) {
  typedef Complex<double> C;
  const C c(x, y), one(1);
  C w(wx, wy);
  for (int it = 0; it < 30; it++) {
    C z = w, dz = one;
    for (int k = 0; k < p; k++) { dz = 2.0 * (z * dz); z = z * z + c; }
    const C step = (z - w) * inv(dz - one);
    w -= step;
    if (sqr_abs(step) < 1e-28 * (1 + sqr_abs(w))) {
      // Derivatives of F = f^p at the periodic point w: A = F_z, B = F_c, Cz = F_zz, D = F_zc
      C z2 = w, A = one, B, Cz, D;
      for (int k = 0; k < p; k++) {
        const C nA = 2.0 * (z2 * A), nB = 2.0 * (z2 * B) + one;
        const C nC = 2.0 * (A * A + z2 * Cz), nD = 2.0 * (A * B + z2 * D);
        A = nA; B = nB; Cz = nC; D = nD;
        z2 = z2 * z2 + c;
      }
      const double a2 = sqr_abs(A);
      if (!(a2 < 1 - 1e-9)) return 0;
      return (1 - a2) / (4 * abs(D + Cz * B * inv(one - A)));
    }
  }
  return 0;
}

// Interior distance at the minimal period dividing p for which Newton finds an attracting cycle
ORBIT_COLD __host__ __device__ static double interior_distance(const double x, const double y, const double wx,
                                                           const double wy, const int p) {
  if (!attracting_cycle(x, y, wx, wy, p)) return 0;
  for (int q = 1; q <= p; q++)
    if (p % q == 0 && attracting_cycle(x, y, wx, wy, q)) return interior_distance_exact(x, y, wx, wy, q);
  return 0;
}

// Classification with distance estimates, for certifying whole cells.  Exterior: log2 of the Green's function
// and the Koebe lower bound on dist(c, M), (1 - e^-g) / (4 |∇g|) ≈ |z_n| log|z_n| / (4 |dz_n/dc|).  Interior
// (attracting cycle found by Newton): the Koebe lower bound (1 - |A|^2) / (4 |D + C B / (1 - A)|) on the
// distance to the component's boundary, where A, B, C, D are derivatives of f^p at the periodic point.
// Zero if unknown.
struct EscapeDE {
  Escape e;
  double dist = 0;  // Exterior: lower bound on dist(c, M).  Interior: lower bound on distance to ∂(component).
};

// Result for an orbit escaping at step n with |z|^2 = r2 and dz/dc = 2^dexp (dx + i dy)
ORBIT_COLD __host__ __device__ static EscapeDE escaped(const int32_t n, const double r2, const double dx,
                                                       const double dy, const int32_t dexp) {
  EscapeDE r;
  const double lz = 0.5 * std::log(r2);
  r.e = {n, std::log2(lz) - double(n - 1), 0, n, r2};
  // Koebe: dist(c, M) ≥ (1 - e^-g) / (4 |∇g|), with g = log|z_n| / 2^(n-1) and
  // |∇g| = |dz_n/dc| / (|z_n| 2^(n-1)).  So dist ≥ (1 - e^-g) 2^(n-1) |z_n| / (4 |dz_n/dc|), where
  // (1 - e^-g) 2^(n-1) = lz for small g.
  const double g = std::exp2(std::log2(lz) - double(n - 1));
  const double scale = g > 1e-8 ? -std::expm1(-g) / g : 1 - g / 2;  // (1 - e^-g) / g
  r.dist = std::ldexp(scale * std::sqrt(r2) * lz / (4 * std::hypot(dx, dy)), int(-(dexp < 100000 ? dexp : 100000)));
  return r;
}

// escape_de as a resumable state machine, like Orbit
struct OrbitDE {
  double x, y, zx, zy, dx, dy, min_r2, cx, cy;
  int32_t n, dexp, candidate, next_newton, check_n, next_check;  // 32 bits to save registers: max_iter < 2^31
  EscapeDE r;  // Result, once done

  // Start at c = x + iy, with Newton interior certificates attempted at step first_newton and each doubling.
  // Returns true if already decided (the cardioid or period 2 disk, whose interior distance is computed here).
  __host__ __device__ bool start(const double x_, const double y_, const int64_t first_newton = 64) {
    x = x_; y = y_;
    r = EscapeDE();
    if (in_cardioid_or_disk(x, y)) {
      const bool disk = (x + 1) * (x + 1) + y * y <= 1.0 / 16;
      r.e = {-1, -INFINITY, disk ? 2 : 1, 0};
      // Start Newton from an orbit point near the attracting cycle (not a fixed guess, which can converge to
      // a repelling cycle instead)
      double wx = x, wy = y;
      for (int k = 0; k < 256; k++) {
        const double t = wx * wx - wy * wy + x;
        wy = 2 * wx * wy + y;
        wx = t;
      }
      r.dist = interior_distance(x, y, wx, wy, disk ? 2 : 1);
      return true;
    }
    // Iterate z and dz/dc together.  dz/dc grows like 2^n on escaping orbits, so rescale it and carry a
    // binary exponent to avoid overflow.
    zx = x; zy = y; dx = 1; dy = 0; dexp = 0;
    min_r2 = zx * zx + zy * zy;
    candidate = 1;
    next_newton = int32_t(first_newton < (int64_t(1) << 30) ? first_newton : int64_t(1) << 30);  // ≥ 2^30: never
    cx = zx; cy = zy; check_n = 1; next_check = 16;  // Brent checkpoint
    n = 1;
    return false;
  }

  // Iterate at most `budget` steps.  Returns true when done, with the result in r.  As in Orbit::run, squares
  // are carried between steps and the rare checks (Newton, Brent, checkpoint) run at steps divisible by 8.
  __host__ __device__ bool run(const int64_t max_iter, const int64_t budget) {
    double zx = this->zx, zy = this->zy, dx = this->dx, dy = this->dy, min_r2 = this->min_r2;
    double zx2 = zx * zx, zy2 = zy * zy, r2 = zx2 + zy2;
    int32_t n = this->n, dexp = this->dexp, candidate = this->candidate;
    double unit = std::ldexp(1.0, int(-(dexp < 2000 ? dexp : 2000)));  // The +1 in dz/dc, rescaled
    const int32_t end = int32_t(n + budget < max_iter + 1 ? n + budget : max_iter + 1);
    bool done = true;
    while (n < end) {
      const int32_t block = (n | 7) + 1 < end ? (n | 7) + 1 : end;  // Up to the next multiple of 8
      for (; n < block; n++) {
        if (r2 > 18446744073709551616.0) {
          r = escaped(n, r2, dx, dy, dexp);
          goto finish;
        }
        {
          // dz/dc ← 2 z dz/dc + 1, rescaled by 2^-dexp to avoid overflow
          const double ndx = 2 * (zx * dx - zy * dy) + unit, ndy = 2 * (zx * dy + zy * dx);
          dx = ndx; dy = ndy;
          if (dx * dx + dy * dy > 1e200) [[unlikely]] {
            dx = std::ldexp(dx, -256); dy = std::ldexp(dy, -256); dexp += 256;
            unit = std::ldexp(1.0, int(-(dexp < 2000 ? dexp : 2000)));
          }
        }
        {
          const double xy = zx * zy;
          zx = zx2 - zy2 + x;
          zy = xy + xy + y;
          zx2 = zx * zx; zy2 = zy * zy; r2 = zx2 + zy2;
          if (r2 < min_r2) { min_r2 = r2; candidate = n + 1; }  // (zx, zy) is now z_{n+1}
        }
      }
      if (n & 7) continue;  // Partial block at the end of the burst
      // (zx, zy) is z_n, n a multiple of 8
      if (n > next_newton) [[unlikely]] {
        while (next_newton < n) next_newton *= 2;
        if (candidate <= 4096) {
          const double b = interior_distance(x, y, zx, zy, int(candidate));
          if (b > 0) { r.e = {-1, -INFINITY, 0, n}; r.dist = b; goto finish; }
        }
      }
      {
        // Brent fallback: converged to a cycle; recover its period by iterating until the orbit returns
        const double ex = zx - cx, ey = zy - cy;
        if (ex * ex + ey * ey < 1e-26) [[unlikely]] {
          const int32_t lag = n - check_n;
          double wx = zx, wy = zy;
          r.e = {-1, -INFINITY, 0, n};
          for (int32_t q = 1; q <= (lag < 65536 ? lag : 65536); q++) {
            const double t2 = wx * wx - wy * wy + x;
            wy = 2 * wx * wy + y;
            wx = t2;
            const double fx = wx - zx, fy = wy - zy;
            if (fx * fx + fy * fy < 1e-20) { r.dist = interior_distance(x, y, zx, zy, int(q)); break; }
          }
          goto finish;
        }
      }
      if (n > next_check) { cx = zx; cy = zy; check_n = n; next_check *= 2; }
    }
    if (n > max_iter) r.e = {-1, -INFINITY, 0, max_iter};
    else done = false;
  finish:
    this->zx = zx; this->zy = zy; this->dx = dx; this->dy = dy; this->min_r2 = min_r2;
    this->n = n; this->dexp = dexp; this->candidate = candidate;
    return done;
  }
};

}  // namespace mandelbrot
