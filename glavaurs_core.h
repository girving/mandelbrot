// The general-root Lavaurs model's evaluation core (glavaurs.h), callable on host and device, at any precision
//
// A trivially copyable snapshot of the model (λ, A, the Fatou coefficients a_j, β, the gate) with the maps built from
// it: the Fatou series, petals, Φ_a, Ψ_k, the return map R_σ = f^n g_σ^r f with its jet, Newton for centers and
// multiplier points, Θ_r, and the horn map.  Templated on the real scalar S: double (GLCore, shared by CPU and GPU
// runs) or Expansion<n> (expansion_math.h), built by gl_make_core<S> from the Fatou coefficients solved in acb.
// Control flow (petal indices, step counts, branch choices) runs in double; values carry S.
#pragma once

#include "complex.h"
#include "cutil.h"
#include "expansion_math.h"
#include <cmath>
#include <math.h>
#include <vector>
namespace mandelbrot {

namespace glcore {
template<class S> __host__ __device__ static inline double dbl(const S x) { return double(x); }
// |z|, arg z, 1/z, polar, wrap: libm in double (as the validated double runs), expansion_math otherwise
__host__ __device__ static inline double cabs(const Complex<double> z) { return ::hypot(z.r, z.i); }
template<int n> __host__ __device__ static inline Expansion<n> cabs(const Complex<Expansion<n>> z) {
  return gl_sqrt(z.r * z.r + z.i * z.i);
}
template<class S> __host__ __device__ static inline S carg(const Complex<S> z) { return gl_atan2(z.i, z.r); }
__host__ __device__ static inline Complex<double> cinv(const Complex<double> z) { return inv(z); }
template<int n> __host__ __device__ static inline Complex<Expansion<n>> cinv(const Complex<Expansion<n>> z) {
  const Expansion<n> m = inv(z.r * z.r + z.i * z.i);
  return Complex<Expansion<n>>(z.r * m, -(z.i * m));
}
template<class S> __host__ __device__ static inline Complex<S> cdiv(const Complex<S> a, const Complex<S> b) {
  return a * cinv(b);
}
template<class S> __host__ __device__ static inline Complex<S> cpolar(const S r, const S t) {
  S s, c;
  gl_sincos(t, s, c);
  return Complex<S>(r * c, r * s);
}
__host__ __device__ static inline double wrap(const double t) { return ::remainder(t, 2 * M_PI); }   // to [-π, π]
template<int n> __host__ __device__ static inline Expansion<n> wrap(const Expansion<n> t) {
  const Expansion<n> tp = twice(gl_pi<Expansion<n>>());
  return t - tp * ::nearbyint(double(t) / (2 * M_PI));
}
template<class S> __host__ __device__ static inline Complex<S> cpow(Complex<S> z, int n) {   // integer power
  const bool neg = n < 0;
  if (neg) n = -n;
  Complex<S> r(1);
  while (n) { if (n & 1) r = r * z; z = sqr(z); n >>= 1; }
  return neg ? cinv(r) : r;
}
template<class S> __host__ __device__ static inline S rpow(const S x, const double e) {   // x^e for x > 0
  return gl_exp(gl_log(x) * S(e));
}
__host__ __device__ static inline double rpow(const double x, const double e) { return ::pow(x, e); }
// The relative step size at which Newton iterations stop
template<class S> __host__ __device__ static inline double eps_of();
template<> __host__ __device__ inline double eps_of<double>() { return 1e-16; }
template<> __host__ __device__ inline double eps_of<Expansion<2>>() { return 1e-31; }
template<> __host__ __device__ inline double eps_of<Expansion<3>>() { return 1e-46; }
template<> __host__ __device__ inline double eps_of<Expansion<4>>() { return 1e-61; }
}  // namespace glcore

template<class S> struct GLJetT {
  typedef Complex<S> C;
  C v, w, s, ww, ws;  // value, ∂w, ∂σ, ∂ww, ∂wσ
  __host__ __device__ GLJetT apply(const C f0, const C f1, const C f2) const {
    return GLJetT{f0, f1 * w, f1 * s, f2 * sqr(w) + f1 * ww, f2 * w * s + f1 * ws};
  }
};
typedef GLJetT<double> GLJet;

template<class S> struct GLCoreT {
  typedef Complex<S> Cd;
  static constexpr int kMaxA = 64;
  static constexpr double kRepel = 400;   // Re ζ' ≤ -max(kRepel, 2|Im ζ'|) for the repelling local inverse
  int p, q, side, kv, N;
  double r0;                              // |w| below which the attracting series is used
  S argA, argmA;                          // arg A, arg(-A)
  Cd lam, A, v, crit, beta, tau;          // τ = 2πiβ/q
  Cd a[kMaxA];                            // a_j for j = -q..N at index j + q

  __host__ __device__ static S pi() { return gl_pi<S>(); }
  __host__ __device__ static S sval(const double x) { return S(x); }

  __host__ __device__ int exit_petal(const int k) const {
    using glcore::dbl;
    const double att = -dbl(argmA) / q + 2 * M_PI * k / q + side * M_PI / q;
    const double base = -dbl(argA) / q;
    const int j = int(::lround((att - base) / (2 * M_PI / q)));
    return ((j % q) + q) % q;
  }

  __host__ __device__ Cd transit_shift(const int k) const {
    const int j = (exit_petal(k) - k) - (exit_petal(kv) - kv);
    return Cd(S(double(j))) * tau;
  }

  // Φ, Φ', Φ'' by the series with L(w) = log|w| + i (θ + wrap(arg w - θ)), continuous within π of the axis θ (a fixed
  // integer branch of the principal log(w^q) would jump where arg(w^q) = π, which for q = 2 is the repelling axis)
  __host__ __device__ void series(const Cd w, const S ax, Cd& s, Cd& d, Cd& dd) const {
    using namespace glcore;
    const Cd L(gl_log(cabs(w)), ax + wrap(carg(w) - ax));
    const Cd iw = cinv(w), iw2 = sqr(iw);
    Cd S_ = beta * L, D = beta * iw, DD = -(beta * iw2);
    Cd pw = cpow(w, -q);   // w^j for j = -q
    for (int j = -q; j <= N; j++, pw = pw * w) {
      if (!j) continue;
      const Cd c = a[j + q];
      S_ += c * pw;
      D += Cd(S(double(j))) * (c * pw * iw);
      DD += Cd(S(double(j) * double(j - 1))) * (c * pw * iw2);
    }
    s = S_; d = D; dd = DD;
  }

  // The petal (kind = -1 attracting, +1 repelling) whose sector contains w, or -1 (decided in double); its axis angle
  __host__ __device__ int petal(const Cd w, const int kind) const {
    using glcore::dbl;
    const Complex<double> wd(dbl(w.r), dbl(w.i)), Ad(dbl(A.r), dbl(A.i));
    const Complex<double> z = Ad * glcore::cpow(wd, q) * Complex<double>(double(kind));
    const double t = ::atan2(z.i, z.r);
    if (::fabs(t) > 0.6 * M_PI) return -1;
    const double base = -(kind > 0 ? dbl(argA) : dbl(argmA)) / q;
    const int k = int(::lround((::atan2(wd.i, wd.r) - base) / (2 * M_PI / q)));
    return ((k % q) + q) % q;
  }
  __host__ __device__ S axis(const int kind, const int k) const {
    return -(kind > 0 ? argA : argmA) / S(double(q)) + twice(pi()) * S(double(k)) / S(double(q));
  }

  // Attracting coordinate with derivatives and the entering petal
  __host__ __device__ bool phi_a(Cd w, Cd& s, Cd& d, Cd& dd, int& pet, const int max_steps = 1 << 20) const {
    using namespace glcore;
    Cd d1(1), d2(0);
    for (int n = 0; n < max_steps; n++) {
      if (::hypot(dbl(w.r), dbl(w.i)) < r0) {
        const int k = petal(w, -1);
        if (k >= 0) {
          Cd ss, sd, sdd;
          series(w, axis(-1, k), ss, sd, sdd);
          s = ss - Cd(S(double(n)) / S(double(q)));
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
      if (::hypot(dbl(w.r), dbl(w.i)) > 10) return false;
    }
    return false;
  }

  // Repelling parametrization of petal k with derivatives
  __host__ __device__ bool psi(const Cd zeta, const int k, Cd& w, Cd& d, Cd& dd) const {
    using namespace glcore;
    const double zr = dbl(zeta.r), zi = dbl(zeta.i);
    if (!(::hypot(zr, zi) < 1e6)) return false;
    const double mm = ::ceil(zr + ::fmax(kRepel, 2 * ::fabs(zi)));
    const long long m = mm > 0 ? (long long)mm : 0;
    if (m > 200000) return false;
    const Cd zl = zeta - Cd(S(double(m)));
    // leading term a_{-q} w^{-q} ≈ ζ': the q-th root in repelling petal k (a double guess; Newton refines)
    const double base = -dbl(argA) / q + 2 * M_PI * k / q;
    const Complex<double> ratio = glcore::cdiv(Complex<double>(dbl(a[0].r), dbl(a[0].i)), Complex<double>(dbl(zl.r), dbl(zl.i)));
    const Complex<double> r = glcore::cpolar(::pow(glcore::cabs(ratio), 1.0 / q), ::atan2(ratio.i, ratio.r) / q);
    Complex<double> best = r;
    double bd = 10;
    for (int t = 0; t < q; t++) {
      const Complex<double> u = r * glcore::cpolar(1.0, 2 * M_PI * t / q);
      const double dlt = ::fabs(glcore::wrap(::atan2(u.i, u.r) - base));
      if (dlt < bd) { bd = dlt; best = u; }
    }
    Cd u(S(best.r), S(best.i));
    const S br = axis(1, k);
    Cd s, sd, sdd;
    const double tol = eps_of<S>();
    for (int it = 0; it < 80; it++) {
      series(u, br, s, sd, sdd);
      const Cd step = cdiv(s - zl, sd);
      u = u - step;
      if (::hypot(dbl(step.r), dbl(step.i)) < tol * ::hypot(dbl(u.r), dbl(u.i))) break;
    }
    series(u, br, s, sd, sdd);
    Cd d1 = cinv(sd);
    Cd d2 = -(sdd * d1 * sqr(d1));
    for (long long i = 0; i < q * m; i++) {
      const Cd df = lam + twice(u);
      d2 = Cd(2) * sqr(d1) + df * d2;
      d1 = df * d1;
      u = lam * u + sqr(u);
      if (::hypot(dbl(u.r), dbl(u.i)) > 10) return false;
    }
    w = u; d = d1; dd = d2;
    return true;
  }

  // R_σ = f^n g_σ^r f at w with its jet in (w, σ)
  __host__ __device__ bool return_map(const int r, const int n, const Cd w, const Cd sigma, GLJetT<S>& x) const {
    x = GLJetT<S>{w, Cd(1), Cd(0), Cd(0), Cd(0)};
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
      GLJetT<S> x;
      if (!return_map(r, n, w, sigma, x)) return INFINITY;
      const Cd F1 = x.v - w, F2 = x.w - mu, a_ = x.w - Cd(1), b = x.s, d = x.ww, e = x.ws, det = a_ * e - b * d;
      const Cd dw = cdiv(F1 * e - b * F2, det), ds = cdiv(a_ * F2 - d * F1, det);
      w = w - dw;
      sigma = sigma - ds;
      last = ::hypot(dbl(dw.r), dbl(dw.i)) + ::hypot(dbl(ds.r), dbl(ds.i));
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

  // Θ_r with σ frozen at sadd in the transits after the first; the slice crosses the fiber at the given speed,
  // p_1 = ζ0 + sadd + speed (σ - sadd).  dfz = dΘ/dσ of the frozen composition, dpar the parameter derivative
  // recursion along the frozen orbit (dp_{i+1} = H' dp_i + 1), and the last state (p_r, petal)
  __host__ __device__ bool theta_frozen(const Cd sigma, const Cd sadd, const int r, Cd& th, Cd& dfz, Cd& dpar, Cd& hprod,
                                        Cd& pr, int& pet_out, const Cd speed = Cd(1)) const {
    Cd s0, ds, dds;
    int pet;
    if (!phi_a(v, s0, ds, dds, pet)) return false;
    Cd pp = s0 + sadd + speed * (sigma - sadd), dz = speed, dq = speed, hp(1);
    for (int i = 1; i < r; i++) {
      Cd q0, q1, q2, a0, a1, a2;
      int pet2;
      if (!psi(pp, exit_petal(pet), q0, q1, q2)) return false;
      if (!phi_a(q0, a0, a1, a2, pet2)) return false;
      const Cd h1 = a1 * q1;
      dz = h1 * dz;
      dq = h1 * dq + Cd(1);
      hp = hp * h1;
      pp = a0 + sadd + transit_shift(pet2);
      pet = pet2;
    }
    th = pp - s0; dfz = dz; dpar = dq; hprod = hp; pr = pp; pet_out = pet;
    return true;
  }
};
typedef GLCoreT<double> GLCore;

// The frozen fiber chain from x: p_1 = x, p_{i+1} = H(p_i) + sadd (i < r), with H' and H'' at each p_i (i < r);
// false if a step fails or r > kMaxChain
static constexpr int kMaxChain = 32;
template<class S> __host__ __device__ static inline bool gl_fiber_chain(const GLCoreT<S>& L, const Complex<S> x,
    const Complex<S> sadd, const int r, Complex<S>* h1, Complex<S>* h2, Complex<S>& pr) {
  typedef Complex<S> Cd;
  if (r > kMaxChain) return false;
  Cd s0, ds, dds;
  int pet;
  if (!L.phi_a(L.v, s0, ds, dds, pet)) return false;
  Cd pp = x;
  for (int i = 1; i < r; i++) {
    Cd h, dh, ddh;
    int pet2;
    if (!L.horn(pp, pet, h, dh, ddh, pet2)) return false;
    h1[i] = dh; h2[i] = ddh;
    pp = h + sadd;
    pet = pet2;
  }
  pr = pp;
  return true;
}

// One children job (glavaurs children): the r-transit children of the source u over the target c, shifts jlo..jhi
struct GLChildJob { int r, nc, jlo, jhi; Complex<double> u, c; double rad; };
// One verified child: job index, shift, start index, center, w = |Θ_r' Π H'|^-2, sat (within 4 radii of the source)
struct GLChild { int job, j, start; Complex<double> s; double w; int sat; };

// The children of all jobs on the GPU (same algorithm and output as the CPU mode); false if CUDA is unavailable
bool glavaurs_gpu_children(const GLCore& core, const std::vector<GLChildJob>& jobs, int starts, std::vector<GLChild>& out);

// The model at precision S (double, Expansion<2>, Expansion<3>): the Fatou coefficients a_j (j = -q..N) and β solved
// in acb at 512 bits and rounded to S, λ, A (with its canonical arg), the gate, r0, and kv (glavaurs.cc)
template<class S> GLCoreT<S> gl_make_core(int p, int q, int side, int N);

// A component's area at precision S by boundary tracing with Nb Fourier points (glavaurs.cc; host only): Newton
// tolerance tol, boundary points accepted below acc; conv = relative change from the half grid, cusp = |σ'(1)|/|a_1|
template<class S> bool gl_area(const GLCoreT<S>& L, int r, int n, Complex<S> guess, Complex<S>& cen, S& area,
                               double& conv, double& cusp, int Nb, double tol, double acc);

}  // namespace mandelbrot
