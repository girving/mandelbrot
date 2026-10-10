// The general-root Lavaurs model's evaluation core (glavaurs.h), callable on host and device
//
// A trivially copyable snapshot of the model (λ, A, the Fatou coefficients a_j, β, the gate) with the maps built from
// it: the Fatou series, petals, Φ_a, Ψ_k, the return map R_σ = f^n g_σ^r f with its jet, Newton for centers and
// multiplier points, Θ_r, and the horn map.  GeneralLavaurs (glavaurs.cc) computes the coefficients in acb and
// delegates its evaluation here, so CPU and GPU runs execute the same code.
#pragma once

#include "complex.h"
#include "cutil.h"
#include <cmath>
#include <math.h>
#include <vector>
namespace mandelbrot {

namespace glcore {
typedef Complex<double> Cd;
__host__ __device__ static inline double cabs(const Cd z) { return ::hypot(z.r, z.i); }
__host__ __device__ static inline double carg(const Cd z) { return ::atan2(z.i, z.r); }
__host__ __device__ static inline Cd cdiv(const Cd a, const Cd b) { return a * inv(b); }
__host__ __device__ static inline Cd cpolar(const double r, const double t) { return Cd(r * ::cos(t), r * ::sin(t)); }
__host__ __device__ static inline double wrap(const double t) { return ::remainder(t, 2 * M_PI); }   // to [-π, π]
__host__ __device__ static inline Cd cpow(Cd z, int n) {   // integer power by squaring
  const bool neg = n < 0;
  if (neg) n = -n;
  Cd r(1);
  while (n) { if (n & 1) r = r * z; z = sqr(z); n >>= 1; }
  return neg ? inv(r) : r;
}
}  // namespace glcore

struct GLJet {
  Complex<double> v, w, s, ww, ws;  // value, ∂w, ∂σ, ∂ww, ∂wσ
  __host__ __device__ GLJet apply(const Complex<double> f0, const Complex<double> f1, const Complex<double> f2) const {
    return GLJet{f0, f1 * w, f1 * s, f2 * sqr(w) + f1 * ww, f2 * w * s + f1 * ws};
  }
};

struct GLCore {
  typedef Complex<double> Cd;
  static constexpr int kMaxA = 64;
  static constexpr double kRepel = 400;   // Re ζ' ≤ -max(kRepel, 2|Im ζ'|) for the repelling local inverse
  int p, q, side, kv, N;
  double r0;                              // |w| below which the attracting series is used
  double argA, argmA;                     // arg A, arg(-A)
  Cd lam, A, v, crit, beta, tau;          // τ = 2πiβ/q
  Cd a[kMaxA];                            // a_j for j = -q..N at index j + q

  __host__ __device__ int exit_petal(const int k) const {
    const double att = -argmA / q + 2 * M_PI * k / q + side * M_PI / q;
    const double base = -argA / q;
    const int j = int(::lround((att - base) / (2 * M_PI / q)));
    return ((j % q) + q) % q;
  }

  __host__ __device__ Cd transit_shift(const int k) const {
    const int j = (exit_petal(k) - k) - (exit_petal(kv) - kv);
    return Cd(double(j)) * tau;
  }

  // Φ, Φ', Φ'' by the series with L(w) = log|w| + i (θ + wrap(arg w - θ)), continuous within π of the axis θ (a fixed
  // integer branch of the principal log(w^q) would jump where arg(w^q) = π, which for q = 2 is the repelling axis)
  __host__ __device__ void series(const Cd w, const double ax, Cd& s, Cd& d, Cd& dd) const {
    using namespace glcore;
    const Cd L(::log(cabs(w)), ax + wrap(carg(w) - ax));
    const Cd iw = inv(w), iw2 = sqr(iw);
    Cd S = beta * L, D = beta * iw, DD = -(beta * iw2);
    Cd pw = cpow(w, -q);   // w^j for j = -q
    for (int j = -q; j <= N; j++, pw = pw * w) {
      if (!j) continue;
      const Cd c = a[j + q];
      S += c * pw;
      D += Cd(double(j)) * (c * pw * iw);
      DD += Cd(double(j) * double(j - 1)) * (c * pw * iw2);
    }
    s = S; d = D; dd = DD;
  }

  // The petal (kind = -1 attracting, +1 repelling) whose sector contains w, or -1; and its axis angle
  __host__ __device__ int petal(const Cd w, const int kind) const {
    using namespace glcore;
    const double t = carg(A * glcore::cpow(w, q) * Cd(double(kind)));
    if (::fabs(t) > 0.6 * M_PI) return -1;
    const double base = -(kind > 0 ? argA : argmA) / q;
    const int k = int(::lround((carg(w) - base) / (2 * M_PI / q)));
    return ((k % q) + q) % q;
  }
  __host__ __device__ double axis(const int kind, const int k) const {
    return -(kind > 0 ? argA : argmA) / q + 2 * M_PI * k / q;
  }

  // Attracting coordinate with derivatives and the entering petal
  __host__ __device__ bool phi_a(Cd w, Cd& s, Cd& d, Cd& dd, int& pet, const int max_steps = 1 << 20) const {
    using namespace glcore;
    Cd d1(1), d2(0);
    for (int n = 0; n < max_steps; n++) {
      if (cabs(w) < r0) {
        const int k = petal(w, -1);
        if (k >= 0) {
          Cd ss, sd, sdd;
          series(w, axis(-1, k), ss, sd, sdd);
          s = ss - Cd(double(n) / q);
          d = sd * d1;
          dd = sdd * sqr(d1) + sd * d2;
          pet = k;
          return true;
        }
      }
      const Cd df = lam + twice(w);
      d2 = Cd(2) * sqr(d1) + df * d2;
      d1 = df * d1;
      w = lam * w + sqr(w);
      if (cabs(w) > 10) return false;
    }
    return false;
  }

  // Repelling parametrization of petal k with derivatives
  __host__ __device__ bool psi(const Cd zeta, const int k, Cd& w, Cd& d, Cd& dd) const {
    using namespace glcore;
    if (!(cabs(zeta) < 1e6)) return false;
    const double mm = ::ceil(zeta.r + ::fmax(kRepel, 2 * ::fabs(zeta.i)));
    const long long m = mm > 0 ? (long long)mm : 0;
    if (m > 200000) return false;
    const Cd zl = zeta - Cd(double(m));
    // leading term a_{-q} w^{-q} ≈ ζ': the q-th root in repelling petal k
    const double base = -argA / q + 2 * M_PI * k / q;
    const Cd ratio = cdiv(a[0], zl);
    const Cd r = cpolar(::pow(cabs(ratio), 1.0 / q), carg(ratio) / q);
    Cd best = r;
    double bd = 10;
    for (int t = 0; t < q; t++) {
      const Cd u = r * cpolar(1.0, 2 * M_PI * t / q);
      const double dlt = ::fabs(wrap(carg(u) - base));
      if (dlt < bd) { bd = dlt; best = u; }
    }
    Cd u = best;
    const double br = axis(1, k);
    Cd s, sd, sdd;
    for (int it = 0; it < 60; it++) {
      series(u, br, s, sd, sdd);
      const Cd step = cdiv(s - zl, sd);
      u = u - step;
      if (cabs(step) < 1e-16 * cabs(u)) break;
    }
    series(u, br, s, sd, sdd);
    Cd d1 = inv(sd);
    Cd d2 = -(sdd * d1 * sqr(d1));
    for (long long i = 0; i < q * m; i++) {
      const Cd df = lam + twice(u);
      d2 = Cd(2) * sqr(d1) + df * d2;
      d1 = df * d1;
      u = lam * u + sqr(u);
      if (cabs(u) > 10) return false;
    }
    w = u; d = d1; dd = d2;
    return true;
  }

  // R_σ = f^n g_σ^r f at w with its jet in (w, σ)
  __host__ __device__ bool return_map(const int r, const int n, const Cd w, const Cd sigma, GLJet& x) const {
    x = GLJet{w, Cd(1), Cd(0), Cd(0), Cd(0)};
    x = x.apply(lam * x.v + sqr(x.v), lam + twice(x.v), Cd(2));
    for (int t = 0; t < r; t++) {
      Cd p0, p1, p2, q0, q1, q2;
      int pet;
      if (!phi_a(x.v, p0, p1, p2, pet)) return false;
      x = x.apply(p0, p1, p2);
      x.v = x.v + sigma;
      x.s = x.s + Cd(1);
      if (!psi(x.v + transit_shift(pet), exit_petal(pet), q0, q1, q2)) return false;
      x = x.apply(q0, q1, q2);
    }
    for (int i = 0; i < n; i++) x = x.apply(lam * x.v + sqr(x.v), lam + twice(x.v), Cd(2));
    return true;
  }

  // Newton for (w, σ) with R_σ(w) = w and R_σ'(w) = μ; returns the last step size (∞ on failure)
  __host__ __device__ double newton(const int r, const int n, Cd& w, Cd& sigma, const Cd mu, const double tol,
                                    const int iters = 40) const {
    using namespace glcore;
    double last = INFINITY;
    for (int it = 0; it < iters; it++) {
      GLJet x;
      if (!return_map(r, n, w, sigma, x)) return INFINITY;
      const Cd F1 = x.v - w, F2 = x.w - mu, a_ = x.w - Cd(1), b = x.s, d = x.ww, e = x.ws, det = a_ * e - b * d;
      const Cd dw = cdiv(F1 * e - b * F2, det), ds = cdiv(a_ * F2 - d * F1, det);
      w = w - dw;
      sigma = sigma - ds;
      last = cabs(dw) + cabs(ds);
      if (!(last < 1)) return INFINITY;
      if (last < tol) break;
    }
    return last;
  }

  // The horn map in the consistent coordinate at transit state (p, pet), with H', H'' and the petal entered
  __host__ __device__ bool horn(const Cd p_, const int pet, Cd& h, Cd& dh, Cd& ddh, int& pet2) const {
    Cd q0, q1, q2, a0, a1, a2;
    if (!psi(p_, exit_petal(pet), q0, q1, q2)) return false;
    if (!phi_a(q0, a0, a1, a2, pet2)) return false;
    h = a0 + transit_shift(pet2);
    dh = a1 * q1;
    ddh = a2 * sqr(q1) + a1 * q2;
    return true;
  }

  // Θ_r(σ) = p_r - ζ0 (p_1 = ζ0 + σ, p_{i+1} = H(p_i) + σ), Θ_r', Θ_r'', Π_{i<r} H'(p_i), and the last petal state
  __host__ __device__ bool theta(const Cd sigma, const int r, Cd& th, Cd& dth, Cd& d2th, Cd& hprod,
                                 int* pet_out = nullptr) const {
    Cd s0, ds, dds;
    int pet;
    if (!phi_a(v, s0, ds, dds, pet)) return false;
    Cd pp = s0 + sigma, dp(1), ddp(0), hp(1);
    for (int i = 1; i < r; i++) {
      Cd q0, q1, q2, a0, a1, a2;
      int pet2;
      if (!psi(pp, exit_petal(pet), q0, q1, q2)) return false;   // pp includes transit_shift(pet)
      if (!phi_a(q0, a0, a1, a2, pet2)) return false;
      const Cd h1 = a1 * q1, h2 = a2 * sqr(q1) + a1 * q2;
      ddp = h2 * sqr(dp) + h1 * ddp;
      dp = h1 * dp + Cd(1);
      hp = hp * h1;
      pp = a0 + sigma + transit_shift(pet2);
      pet = pet2;
    }
    th = pp - s0;
    dth = dp;
    d2th = ddp;
    hprod = hp;
    if (pet_out) *pet_out = pet;
    return true;
  }
};

// One children job (glavaurs children): the r-transit children of the source u over the target c, shifts jlo..jhi
struct GLChildJob { int r, nc, jlo, jhi; Complex<double> u, c; double rad; };
// One verified child: job index, shift, start index, center, w = |Θ_r' Π H'|^-2, sat (within 4 radii of the source)
struct GLChild { int job, j, start; Complex<double> s; double w; int sat; };

// The children of all jobs on the GPU (same algorithm and output as the CPU mode); false if CUDA is unavailable
bool glavaurs_gpu_children(const GLCore& core, const std::vector<GLChildJob>& jobs, int starts, std::vector<GLChild>& out);

}  // namespace mandelbrot
