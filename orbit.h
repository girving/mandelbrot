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
#include <cstring>
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
__host__ __device__ static inline bool escaped_below(const int64_t steps, const double r2, const int64_t k) {
  const int64_t d = steps - 1 - k;
  return d >= 6 || (d == 5 && r2 < 6.235149080811617e27);
}

// Rare, heavy paths (Newton certificates, distance bounds) stay out of line, so that the hot iteration loops
// that call them keep few registers on GPUs
#define ORBIT_COLD __attribute__((noinline))

// The type of an orbit's parameter c and of the constants it compares or combines with T: T itself, except for
// expansions (orbit_expansion.h), whose constants are exact doubles, so that cheaper mixed expansion-double
// arithmetic applies.  (Writing constants as double would instead promote float arithmetic to double.)
template<class T> struct OrbitParam { typedef T type; };

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

// Newton certificate settings.  Negative tol2 or margin mean the precision's defaults (OrbitTol).
struct NewtonOptions {
  int iters = 30;             // Newton iterations per attempt
  double close2 = INFINITY;   // Give up after one step unless |f^p(w) - w|^2 < close2
  double tol2 = -1;           // Converged when |step|^2 < tol2 (1 + |w|^2)
  double margin = -1;         // Attracting when |λ|^2 < 1 - margin
  bool best_return = false;   // If the atom-domain candidate fails, also try the best return period (costly on
                              // GPUs, where most Newton attempts are on exterior orbits and fail)
};

// Newton's method for an attracting p-cycle of z → z^2 + c near w.  Returns true if Newton converges to a
// periodic point whose multiplier |(f^p)'(w)| < 1, which certifies that c is in a hyperbolic component.
//
// If close2 is finite, give up after the first iteration unless |f^p(w) - w|^2 < close2: orbit points that
// have not nearly closed up rarely converge, and failures otherwise cost the full iteration count.
template<class T> ORBIT_COLD __host__ __device__ bool attracting_cycle(const typename OrbitParam<T>::type x,
                                                                   const typename OrbitParam<T>::type y, T wx, T wy,
                                                                   const int p,
                                                         const NewtonOptions nw = NewtonOptions()) {
  typedef OrbitTol<T> Tol;
  typedef typename OrbitParam<T>::type P;
  const double tol2 = nw.tol2 < 0 ? Tol::newton : nw.tol2, margin = nw.margin < 0 ? Tol::multiplier : nw.margin,
               close2 = nw.close2;
  for (int it = 0; it < nw.iters; it++) {
    // F(w) = f^p(w) - w, F'(w) = (f^p)'(w) - 1
    T zx = wx, zy = wy, dx = T(1), dy = T(0);
    for (int k = 0; k < p; k++) {
      const T ndx = 2 * (zx * dx - zy * dy), ndy = 2 * (zx * dy + zy * dx);
      dx = ndx; dy = ndy;
      const T t = zx * zx - zy * zy + x;
      zy = 2 * zx * zy + y;
      zx = t;
      if (zx * zx + zy * zy > P(16)) return false;
    }
    const T fx = zx - wx, fy = zy - wy, gx = dx - P(1), gy = dy;
    if (it == 0 && !(double(fx * fx + fy * fy) < close2)) return false;
    const T den = gx * gx + gy * gy;
    if (!(den > P(0))) return false;
    const T sx = (fx * gx + fy * gy) / den, sy = (fy * gx - fx * gy) / den;
    wx -= sx; wy -= sy;
    if (sx * sx + sy * sy < P(tol2) * (P(1) + wx * wx + wy * wy)) {
      // Converged: the multiplier at the periodic point decides
      T mx = T(1), my = T(0), zx2 = wx, zy2 = wy;
      for (int k = 0; k < p; k++) {
        const T nmx = 2 * (zx2 * mx - zy2 * my), nmy = 2 * (zx2 * my + zy2 * mx);
        mx = nmx; my = nmy;
        const T t = zx2 * zx2 - zy2 * zy2 + x;
        zy2 = 2 * zx2 * zy2 + y;
        zx2 = t;
      }
      return mx * mx + my * my < P(1 - margin);
    }
  }
  return false;
}

// Orbit and OrbitDE keep step counts in 32 bits to save GPU registers, so max_iter must be below 2^30.  Building
// with -DMANDELBROT_ORBIT64 makes them 64 bits, for deep runs (max_iter below 2^62), at some cost in registers.
#ifdef MANDELBROT_ORBIT64
typedef int64_t orbit_int;
constexpr int64_t kOrbitNever = int64_t(1) << 62;
#else
typedef int32_t orbit_int;
constexpr int64_t kOrbitNever = int64_t(1) << 30;
#endif

// Fused multiply-add for orbit arithmetic, overloaded per number type (Rounded and Expansion<2> define their own).
// Explicit, so that CPU and GPU agree bit for bit (contraction is off in all builds) while the step uses FMAs.
__host__ __device__ inline double orbit_fma(const double a, const double b, const double c) { return fma(a, b, c); }
__host__ __device__ inline float orbit_fma(const float a, const float b, const float c) { return fmaf(a, b, c); }

// Order key for |z|^2 in the atom-domain minimum: the high 32 bits of the double, which order like the value
// for nonnegative doubles (nan and inf sort above everything).  An integer compare, which GPUs issue alongside
// FP64 work, instead of an FP64 one; ties within 2^-20 relative keep the earlier step.  Only the Newton period
// candidate depends on it, and Newton certifies only true attracting cycles, so classifications cannot.
__host__ __device__ inline int32_t orbit_key(const double r2) {
#ifdef __CUDA_ARCH__
  return __double2hiint(r2);
#else
  uint64_t b;
  memcpy(&b, &r2, sizeof(b));
  return int32_t(b >> 32);
#endif
}
template<class T> __host__ __device__ inline int32_t orbit_key(const T r2) { return orbit_key(double(r2)); }

// Exact doubling (Expansion<2> overloads it with its componentwise twice)
template<class T> __host__ __device__ inline T orbit_twice(const T a) { return a + a; }

// One step z → z^2 + c carrying zy^2 and r2 = |z|^2: 3 FMAs, 2 adds and 1 multiply
template<class T> __host__ __device__ inline void orbit_step(T& zx, T& zy, T& zy2, T& r2,
                                                             const typename OrbitParam<T>::type x,
                                                             const typename OrbitParam<T>::type y) {
  const T t = x - zy2;
  zy = orbit_fma(orbit_twice(zx), zy, y);
  zx = orbit_fma(zx, zx, t);
  zy2 = zy * zy;
  r2 = orbit_fma(zx, zx, zy2);
}

template<class T> struct Orbit {
  typedef typename OrbitParam<T>::type P;
  // State is kept small, since it lives in GPU registers across the persistent loop: the result is encoded in
  // existing fields (see status), and task-wide settings are arguments to run.
  P x, y;         // The parameter c
  T zx, zy;
  T cx, cy;        // Brent checkpoint, refreshed at powers of two.  After escape, cx = |z_n|^2.
  int32_t min_key; // Atom domains: the step where |z_n| reaches a new minimum (orbit_key) is a candidate period
  orbit_int n, next_check, candidate, next_newton;  // max_iter < kOrbitNever
  int32_t status;  // 0 running, 1 escaped at step n, 2 attracting cycle (period in candidate), 3 hit max_iter,
                   // 4 stopped at a Newton step, 5 stopped at an exact return (settle does the rest of the step)

  // Start at c = x + iy, with Newton certificate attempts at step first_newton and then at each doubling.
  // Returns true if already decided (the cardioid or period 2 disk).
  __host__ __device__ bool start(const double x_, const double y_, const int64_t first_newton = 64) {
    // Start at z_0 = 0 (the first step gives z_1 = c exactly), so that 8-step blocks align from the start
    x = P(x_); y = P(y_);
    zx = T(0); zy = T(0); cx = T(x); cy = T(y);
    min_key = INT32_MAX;
    n = 0; next_check = 16; candidate = 1; status = 0;
    next_newton = orbit_int(first_newton < kOrbitNever ? first_newton : kOrbitNever);  // ≥ kOrbitNever: never
    if (in_cardioid_or_disk(x_, y_)) {
      const bool disk = (x_ + 1) * (x_ + 1) + y_ * y_ <= 1.0 / 16;
      status = 2; candidate = disk ? 2 : 1; n = 0;
      return true;
    }
    return false;
  }

  // Stopped for settle
  __host__ __device__ bool pending() const { return status == 4 || status == 5 || status == 8; }
  // Stopped for a settle that should run at once (status 8: a block overflowed), all lanes of a warp together
  __host__ __device__ bool immediate() const { return status == 8; }

  // Run to completion, settling as we go
  __host__ __device__ void finish(const int64_t max_iter, const int max_period = 4096,
                                  const NewtonOptions& nw = NewtonOptions()) {
    for (;;) {
      if (!run(max_iter, max_iter, max_period)) continue;
      if (!pending() || settle(max_iter, max_period, nw)) return;
    }
  }

  // Iterations performed, once done
  __host__ __device__ int64_t iters() const { return status == 3 && n > 1 ? n - 1 : n; }

  // The result as an Escape, once done.  period is the minimal period if found and at most 32, else 0.
  __host__ __device__ Escape result() const {
    switch (status) {
      case 1: return {n, std::log2(0.5 * std::log(double(cx))) - double(n - 1), 0, n, double(cx)};
      case 2: return {-1, -INFINITY, candidate <= 32 ? candidate : 0, n};
      default: return {-1, -INFINITY, 0, n - 1};
    }
  }

  // Brent: the computed orbit has returned exactly (bit for bit) to its checkpoint, so it is periodic and will
  // never escape: below every threshold, exactly as iterating to max_iter would find.  (Near returns prove
  // nothing: slow exterior orbits near nearly neutral repelling cycles make them too.  Converging interior
  // orbits reach an exact floating-point cycle soon after.)  Sets the minimal period if at most 32, else 33.
  __host__ __device__ void exact_cycle() {
    typedef OrbitTol<T> Tol;
    T wx = zx, wy = zy;
    candidate = 33;
    for (int p = 1; p <= 32; p++) {
      const T t2 = wx * wx - wy * wy + x;
      wy = 2 * wx * wy + y;
      wx = t2;
      const T ex = wx - zx, ey = wy - zy;
      if (ex * ex + ey * ey < P(Tol::period)) { candidate = p; break; }
    }
  }

  // The q ≤ max_period minimizing |f^q(z) - z|: the period (or a multiple) once the orbit is near a cycle
  __host__ __device__ int best_return(const int max_period) const { return best_return(zx, zy, max_period); }
  __host__ __device__ int best_return(const T zx, const T zy, const int max_period) const {
    T wx = zx, wy = zy, best = T(INFINITY);
    int q = 0;
    for (int k = 1; k <= max_period; k++) {
      const T t = wx * wx - wy * wy + x;
      wy = 2 * wx * wy + y;
      wx = t;
      const T dx = wx - zx, dy = wy - zy, d = dx * dx + dy * dy;
      if (d < best) { best = d; q = k; }
    }
    return q;
  }

  // Newton at the current point, on the atom-domain candidate and, with nw.best_return, if that fails, on the
  // best return (atom domains need not match components near their boundaries).  If it certifies an attracting cycle, set the
  // minimal period (or 33).
  __host__ __device__ bool newton(const int max_period, const NewtonOptions& nw) {
    if (!(candidate <= max_period && attracting_cycle(x, y, zx, zy, int(candidate), nw))) {
      if (!nw.best_return) return false;
      const int q = best_return(max_period);
      if (!(q && q != candidate && attracting_cycle(x, y, zx, zy, q, nw))) return false;
      candidate = q;
    }
    int period = 0;
    if (candidate <= 32)
      for (int q = 1; q <= int(candidate); q++)
        if (int(candidate) % q == 0 && attracting_cycle(x, y, zx, zy, q)) { period = q; break; }
    candidate = period ? period : 33;  // 33: a period above 32, reported as 0
    return true;
  }

  // Finish a step at which run stopped (status 4 or 5) as a single step would: Newton (if due), then Brent and the
  // checkpoint.  Returns true if the orbit is done; otherwise it can run again.
  __host__ __device__ bool settle(const int64_t max_iter, const int max_period, const NewtonOptions& nw) {
    if (status == 8) {
      // A fast block overflowed: redo it step by step to find the escape step
      status = 0;
      T zy2 = zy * zy, r2 = orbit_fma(zx, zx, zy2);
      const P big = P(18446744073709551616.0);
      const orbit_int block = (n | 7) + 1;
      for (;; n++) {
        if (r2 > big) { cx = r2; status = 1; return true; }
        if (n == block) return false;  // |z|^2 landed exactly on the threshold: keep going (skipping one block's checks)
        orbit_step(zx, zy, zy2, r2, x, y);
        if (n < max_period) {
          const int32_t key = orbit_key(r2);
          if (key < min_key) { min_key = key; candidate = n + 1; }
        }
      }
    }
    const bool due = status == 4;
    status = 0;
    if (due && newton(max_period, nw)) { status = 2; return true; }
    if (zx == cx && zy == cy) { exact_cycle(); status = 2; return true; }
    if (n > next_check) { cx = zx; cy = zy; next_check *= 2; }
    if (n > max_iter) status = 3;
    return status != 0;
  }

  // Iterate at most `budget` steps, trying Newton (options nw) on atom-domain candidates up to max_period.  Returns true when done (see status).  State lives in locals during the loop so that it stays
  // in registers.  Each step does the escape test and z → z^2 + c, carrying the squares, and tracks the
  // atom-domain minimum while n < max_period (later candidates are never used).  The rarer checks (Newton,
  // Brent, checkpoint) run at steps that are multiples of 8: they cost as much as the iteration itself on GPUs,
  // and tying them to absolute step numbers keeps results independent of how the orbit is split into bursts.
  //
  // run stops at Newton steps (status 4) and exact returns (status 5) and leaves them to settle, which callers
  // run next: keeping Newton out of this loop keeps its registers down on GPUs, and lets GPU lanes settle together.
  __host__ __device__ bool run(const int64_t max_iter, const int64_t budget, const int max_period = 4096) {
    T zx = this->zx, zy = this->zy, cx = this->cx, cy = this->cy;
    int32_t min_key = this->min_key;
    T zy2 = zy * zy, r2 = orbit_fma(zx, zx, zy2);
    orbit_int n = this->n, next_check = this->next_check, candidate = this->candidate;
    // End bursts at multiples of 8 (the budget permitting), so that later bursts are whole aligned blocks
    orbit_int end = orbit_int(n + budget < max_iter + 1 ? n + budget : max_iter + 1);
    if (end <= max_iter && (end & ~orbit_int(7)) > n) end &= ~orbit_int(7);
    const P big = P(18446744073709551616.0);  // Escape at |z| > 2^32 so that log|z| is accurate
    while (n < end) {
      const orbit_int block = (n | 7) + 1 < end ? (n | 7) + 1 : end;  // Up to the next multiple of 8
      // Fast path: a whole aligned block with no per-step tests.  Once |z| > 2^32 it only grows (to inf or nan
      // within the block), so a single test at the end finds escapes, and the block is redone step by step to
      // find the exact escape step.  Inside the atom-domain range the block tracks the minimum with selects
      // instead of branches; a block straddling max_period takes the careful path.  Both paths do the same
      // arithmetic, so results do not depend on which ran.
      if (block == n + 8 && !(r2 > big) && (n >= max_period || n + 8 <= max_period)) {
        const T zx0 = zx, zy0 = zy;
        const int32_t min0 = min_key;
        const orbit_int cand0 = candidate;
        const bool track = n < max_period;  // The whole block, since it does not straddle max_period
        // One loop body for every lane: the minimum is tracked with integer selects, so lanes inside and past the
        // atom-domain range do not diverge
#ifdef __CUDA_ARCH__
#pragma unroll
#endif
        for (int s = 0; s < 8; s++) {
          orbit_step(zx, zy, zy2, r2, x, y);
          const int32_t key = orbit_key(r2);
          const bool lower = track && key < min_key;  // (zx, zy) is now z_{n+s+1}
          min_key = lower ? key : min_key;
          candidate = lower ? n + s + 1 : candidate;
        }
        if (r2 < big) {  // False for nan and inf
          n = block;
        } else {
          // It escaped in this block: restore its start and leave the step-by-step search to settle, which the
          // engine runs for all such lanes of a warp together instead of each lane stalling the rest
          zx = zx0; zy = zy0; min_key = min0; candidate = cand0;
          status = 8;
          goto finish;
        }
      }
      for (; n < block; n++) {
        if (r2 > big) {
          cx = r2;
          status = 1;
          goto finish;
        }
        orbit_step(zx, zy, zy2, r2, x, y);
        if (n < max_period) {  // (zx, zy) is now z_{n+1}
          const int32_t key = orbit_key(r2);
          if (key < min_key) { min_key = key; candidate = n + 1; }
        }
      }
      if (n & 7) continue;  // Partial block at the end of the burst
      // (zx, zy) is z_n, n a multiple of 8
      if (n > next_newton) [[unlikely]] {
        while (next_newton < n) next_newton *= 2;
        status = 4;
        goto finish;
      }
      if (zx == cx && zy == cy) [[unlikely]] { status = 5; goto finish; }  // Exact return: settle finds the period
      if (n > next_check) { cx = zx; cy = zy; next_check *= 2; }
    }
    if (n > max_iter) status = 3;
  finish:
    this->zx = zx; this->zy = zy; this->cx = cx; this->cy = cy; this->min_key = min_key;
    this->n = n; this->next_check = next_check; this->candidate = candidate;
    return status != 0;
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
ORBIT_COLD __host__ __device__ static EscapeDE escaped(const int64_t n, const double r2, const double dx,
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
  double x, y, zx, zy, dx, dy, cx, cy;
  orbit_int n, candidate, next_newton, check_n, next_check;  // max_iter < kOrbitNever
  int32_t dexp, min_key;  // min_key: atom-domain minimum of |z|^2 (orbit_key)
  // 0 running, 1 done (result in r), 7 escaped (at step n with |z|^2 = cx; result() computes the distance), and
  // with defer, stopped for settle: 4 Newton step due, 5 Brent fired, 6 cardioid or period 2 disk (whose
  // interior distance settle computes)
  int32_t status;
  EscapeDE r;  // Result, once done

  // Start at c = x + iy, with Newton interior certificates attempted at step first_newton and each doubling.
  // Returns true if already decided: the cardioid or period 2 disk, unless defer (then status 6, for settle).
  __host__ __device__ bool start(const double x_, const double y_, const int64_t first_newton = 64,
                                 const bool defer = false) {
    x = x_; y = y_;
    r = EscapeDE();
    status = 0;
    if (in_cardioid_or_disk(x, y)) {
      status = 6;
      if (defer) return false;
      settle(0);
      return true;
    }
    // Iterate z and dz/dc together.  dz/dc grows like 2^n on escaping orbits, so rescale it and carry a
    // binary exponent to avoid overflow.
    zx = x; zy = y; dx = 1; dy = 0; dexp = 0;
    min_key = orbit_key(zx * zx + zy * zy);
    candidate = 1;
    next_newton = orbit_int(first_newton < kOrbitNever ? first_newton : kOrbitNever);  // ≥ kOrbitNever: never
    cx = zx; cy = zy; check_n = 1; next_check = 16;  // Brent checkpoint
    n = 1;
    return false;
  }

  // The result, once done (status 1 or 7)
  __host__ __device__ EscapeDE result() const { return status == 7 ? escaped(n, cx, dx, dy, dexp) : r; }
  __host__ __device__ int64_t iters() const { return status == 7 ? n : r.e.iters; }

  // The rare checks at a step n divisible by 8, after the Newton step (if due) and Brent (if it fired) were
  // detected by run: Newton, then Brent's period recovery, then the checkpoint.  Also the cardioid/disk start.
  // Returns true if done.
  __host__ __device__ bool settle(const int64_t max_iter, const int max_period = 4096) {
    if (status == 6) {
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
      status = 1;
      return true;
    }
    if (status == 4 && candidate <= max_period) {
      const double b = interior_distance(x, y, zx, zy, int(candidate));
      if (b > 0) { r.e = {-1, -INFINITY, 0, n}; r.dist = b; status = 1; return true; }
    }
    {
      // Brent fallback: converged to a cycle; recover its period by iterating until the orbit returns
      const double ex = zx - cx, ey = zy - cy;
      if (ex * ex + ey * ey < 1e-26) {
        const orbit_int lag = n - check_n;
        double wx = zx, wy = zy;
        r.e = {-1, -INFINITY, 0, n};
        for (int32_t q = 1; q <= int32_t(lag < 65536 ? lag : 65536); q++) {
          const double t2 = wx * wx - wy * wy + x;
          wy = 2 * wx * wy + y;
          wx = t2;
          const double fx = wx - zx, fy = wy - zy;
          if (fx * fx + fy * fy < 1e-20) { r.dist = interior_distance(x, y, zx, zy, int(q)); break; }
        }
        status = 1;
        return true;
      }
    }
    if (n > next_check) { cx = zx; cy = zy; check_n = n; next_check *= 2; }
    if (n > max_iter) { r.e = {-1, -INFINITY, 0, int64_t(max_iter)}; status = 1; return true; }
    status = 0;
    return false;
  }

  // One step of z and dz/dc (rescaled by 2^-dexp, so the +1 is unit): 6 FMAs, 3 adds and 2 multiplies, sharing
  // 2z between z^2 + c and 2 z dz + 1
  __host__ __device__ inline void step(double& zx, double& zy, double& zy2, double& r2, double& dx, double& dy,
                                       const double unit) const {
    const double tzx = zx + zx, tzy = zy + zy;
    const double ndx = fma(tzx, dx, fma(-tzy, dy, unit)), ndy = fma(tzx, dy, tzy * dx);
    dx = ndx; dy = ndy;
    const double t = x - zy2;
    zy = fma(tzx, zy, y);
    zx = fma(zx, zx, t);
    zy2 = zy * zy;
    r2 = fma(zx, zx, zy2);
  }

  // Iterate at most `budget` steps.  Returns true when done, with the result in r, or, with defer, when stopped
  // for settle (status 4 or 5).  As in Orbit::run: whole aligned blocks of 8 steps run without per-step tests
  // (escapes found at the block end and the block redone step by step), the atom-domain minimum is an integer
  // key, and the rare checks (Newton, Brent, checkpoint) run at steps divisible by 8.  dz/dc is rescaled at block
  // ends: from |dz| ≤ 1e100, 8 steps before escape grow it by at most (2^33)^8 < 1e80, and rescaling by exact
  // powers of 2 leaves results independent of when it happens.
  __host__ __device__ bool run(const int64_t max_iter, const int64_t budget, const bool defer = false,
                               const int max_period = 4096) {
    if (status == 6) return true;  // Deferred cardioid/disk start
    double zx = this->zx, zy = this->zy, dx = this->dx, dy = this->dy;
    double zy2 = zy * zy, r2 = fma(zx, zx, zy2);
    orbit_int n = this->n, candidate = this->candidate;
    int32_t dexp = this->dexp, min_key = this->min_key;
    double unit = std::ldexp(1.0, int(-(dexp < 2000 ? dexp : 2000)));  // The +1 in dz/dc, rescaled
    orbit_int end = orbit_int(n + budget < max_iter + 1 ? n + budget : max_iter + 1);
    if (end <= max_iter && (end & ~orbit_int(7)) > n) end &= ~orbit_int(7);  // Bursts end on block boundaries
    const double big = 18446744073709551616.0;
    while (n < end) {
      const orbit_int block = (n | 7) + 1 < end ? (n | 7) + 1 : end;  // Up to the next multiple of 8
      if (block == n + 8 && !(r2 > big)) {
        const double zx0 = zx, zy0 = zy, zy20 = zy2, r20 = r2, dx0 = dx, dy0 = dy;
        const int32_t min0 = min_key;
        const orbit_int cand0 = candidate;
#ifdef __CUDA_ARCH__
#pragma unroll
#endif
        for (int s = 0; s < 8; s++) {
          step(zx, zy, zy2, r2, dx, dy, unit);
          const int32_t key = orbit_key(r2);
          const bool lower = key < min_key;  // (zx, zy) is now z_{n+s+1}
          min_key = lower ? key : min_key;
          candidate = lower ? n + s + 1 : candidate;
        }
        if (r2 < big) {  // False for nan and inf
          n = block;
        } else {
          zx = zx0; zy = zy0; zy2 = zy20; r2 = r20; dx = dx0; dy = dy0; min_key = min0; candidate = cand0;
        }
      }
      for (; n < block; n++) {
        if (r2 > big) {
          // The distance bound is computed by result(), off the hot loop: GPU lanes that finish in the same burst
          // compute it together, instead of each stalling its warp
          cx = r2;
          status = 7;
          goto finish;
        }
        step(zx, zy, zy2, r2, dx, dy, unit);
        const int32_t key = orbit_key(r2);
        if (key < min_key) { min_key = key; candidate = n + 1; }  // (zx, zy) is now z_{n+1}
      }
      if (dx * dx + dy * dy > 1e200) [[unlikely]] {
        dx = std::ldexp(dx, -256); dy = std::ldexp(dy, -256); dexp += 256;
        unit = std::ldexp(1.0, int(-(dexp < 2000 ? dexp : 2000)));
      }
      if (n & 7) continue;  // Partial block at the end of the burst
      // (zx, zy) is z_n, n a multiple of 8
      {
        const bool newton = n > next_newton;
        if (newton) while (next_newton < n) next_newton *= 2;
        const double ex = zx - cx, ey = zy - cy;
        const bool brent = ex * ex + ey * ey < 1e-26;
        if (newton || brent) [[unlikely]] {
          status = newton ? 4 : 5;
          if (defer) goto finish;
          this->zx = zx; this->zy = zy; this->n = n; this->candidate = candidate;
          if (settle(max_iter, max_period)) goto finish;
          continue;  // settle did the checkpoint
        }
      }
      if (n > next_check) { cx = zx; cy = zy; check_n = n; next_check *= 2; }
    }
    if (n > max_iter) { r.e = {-1, -INFINITY, 0, int64_t(max_iter)}; status = 1; }
  finish:
    this->zx = zx; this->zy = zy; this->dx = dx; this->dy = dy; this->min_key = min_key;
    this->n = n; this->dexp = dexp; this->candidate = candidate;
    return status != 0;
  }
};

}  // namespace mandelbrot
