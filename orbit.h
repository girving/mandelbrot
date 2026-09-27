// Resumable escape-time orbits, for batched iteration on CPU lanes and GPU threads
//
// Orbit<T> is escape() as a state machine: start() at a parameter, then call step() until it reports the
// result.  One step is one iteration of z → z^2 + c plus the (rare) attracting-cycle checks, so CPU code can
// interleave many independent orbits for instruction-level parallelism, and GPU threads can refill finished
// orbits with new samples without waiting for the rest of their warp.  T is double or float; tolerances
// scale with the precision.
#pragma once

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
  double g() const { return std::exp2(log2g); }
};

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
template<class T> __host__ __device__ bool attracting_cycle(const T x, const T y, T wx, T wy, const int p) {
  typedef OrbitTol<T> Tol;
  for (int it = 0; it < 30; it++) {
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

template<class T> struct Orbit {
  T x, y, zx, zy;
  T cx, cy;        // Brent checkpoint, refreshed at powers of two
  T min_r2;        // Atom domains: the step where |z_n| reaches a new minimum is a candidate period
  int64_t n, next_check, candidate, next_newton;
  Escape e;        // Result, once done

  // Start at c = x + iy.  Returns true if already decided (the cardioid or period 2 disk).
  __host__ __device__ bool start(const double x_, const double y_) {
    x = T(x_); y = T(y_);
    zx = x; zy = y; cx = x; cy = y;
    min_r2 = zx * zx + zy * zy;
    n = 1; next_check = 16; candidate = 1; next_newton = 64;
    if (in_cardioid_or_disk(x_, y_)) {
      const bool disk = (x_ + 1) * (x_ + 1) + y_ * y_ <= 1.0 / 16;
      e = {-1, -INFINITY, disk ? 2 : 1, 0};
      return true;
    }
    return false;
  }

  // Iterate at most `budget` steps.  Returns true when done, with the result in e.  State lives in locals
  // during the loop so that it stays in registers.
  __host__ __device__ bool run(const int64_t max_iter, const int64_t budget) {
    typedef OrbitTol<T> Tol;
    T zx = this->zx, zy = this->zy, cx = this->cx, cy = this->cy, min_r2 = this->min_r2;
    int64_t n = this->n, next_check = this->next_check, candidate = this->candidate;
    const int64_t end = n + budget < max_iter + 1 ? n + budget : max_iter + 1;
    bool done = true;
    for (; n < end; n++) {
      // Escape at |z| > 2^32 so that log|z| is accurate
      const T r2 = zx * zx + zy * zy;
      if (r2 > T(18446744073709551616.0)) {
        e = {n, std::log2(0.5 * std::log(double(r2))) - double(n - 1), 0, n};
        goto finish;
      }
      const T t = zx * zx - zy * zy + x;
      zy = 2 * zx * zy + y;
      zx = t;
      const T r2n = zx * zx + zy * zy;
      if (r2n < min_r2) { min_r2 = r2n; candidate = n + 1; }  // (zx, zy) is now z_{n+1}
      if (n == next_newton) [[unlikely]] {
        next_newton *= 2;
        if (candidate <= 4096 && attracting_cycle(x, y, zx, zy, int(candidate))) {
          // Report the minimal period if it is small
          int period = 0;
          if (candidate <= 32)
            for (int q = 1; q <= int(candidate); q++)
              if (int(candidate) % q == 0 && attracting_cycle(x, y, zx, zy, q)) { period = q; break; }
          e = {-1, -INFINITY, period, n};
          goto finish;
        }
      }
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
      if (n == next_check) { cx = zx; cy = zy; next_check *= 2; }
    }
    if (n > max_iter) e = {-1, -INFINITY, 0, max_iter};
    else done = false;
  finish:
    this->zx = zx; this->zy = zy; this->cx = cx; this->cy = cy; this->min_r2 = min_r2;
    this->n = n; this->next_check = next_check; this->candidate = candidate;
    return done;
  }
};

}  // namespace mandelbrot
