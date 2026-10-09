// The Lavaurs model at the 1/2 root of the cardioid

#include "lavaurs.h"
#include "acb_cc.h"
#include "arb_cc.h"
#include "debug.h"
#include "expansion_arith.h"
#include "nearest.h"
#include <cmath>
#include <complex>
#include <flint/acb.h>
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpz.h>
namespace mandelbrot {

namespace {

typedef Expansion<2> E2;
typedef Complex<double> Cd;
typedef Complex<E2> Ce;

template<class S> inline Complex<S> cdiv(const Complex<S> a, const Complex<S> b) {
  const S d = sqr(b.r) + sqr(b.i);
  const Complex<S> n = a * conj(b);
  return Complex<S>(n.r / d, n.i / d);
}
inline Cd to_cd(const Cd z) { return z; }
inline Cd to_cd(const Ce z) { return Cd(double(z.r), double(z.i)); }
inline double cabs(const Cd z) { return std::sqrt(z.r * z.r + z.i * z.i); }
template<class S> inline Complex<S> from_cd(const Cd z) { return Complex<S>(S(z.r), S(z.i)); }

// Complex log: libm for double, arb (exact input, correctly rounded output) for Expansion<2>
inline Cd clog_s(const Cd z) {
  const auto l = std::log(std::complex<double>(z.r, z.i));
  return Cd(l.real(), l.imag());
}
inline Ce clog_s(const Ce z) {
  Acb x, y;
  const Arb re = exact_arb(z.r), im = exact_arb(z.i);
  acb_set_arb_arb(x, re.x, im.x);
  acb_log(y, x, 256);
  const auto r = round_nearest<E2>(y, 256);
  slow_assert(r, "clog_s: rounding failed");
  return *r;
}

// Exact series coefficients [C, a_1, ..., a_{M-1}] from orders 1..M of Φ(F(w)) - Φ(w) - 1/2 = 0, with A = B = 1/4:
// order n: (n + 2)/4 + C log(1-w)_n + Σ_j a_j [w^n](F^j - w^j) = 0, [w^n] F^j = (-1)^n binom(j, n - j)
void fatou_coefficients(const int M, fmpq_mat_t X) {
  fmpq_mat_t A, B;
  fmpq_mat_init(A, M, M);
  fmpq_mat_init(B, M, 1);
  fmpz_t bin;
  fmpz_init(bin);
  for (int n = 1; n <= M; n++) {
    fmpq_set_si(fmpq_mat_entry(A, n - 1, 0), -1, n);
    for (int j = 1; j < M; j++) {
      fmpq_zero(fmpq_mat_entry(A, n - 1, j));
      if (n >= j && n - j <= j) {
        fmpz_bin_uiui(bin, j, n - j);
        if (n % 2) fmpz_neg(bin, bin);
        fmpq_set_fmpz(fmpq_mat_entry(A, n - 1, j), bin);
      }
      if (n == j) fmpq_sub_si(fmpq_mat_entry(A, n - 1, j), fmpq_mat_entry(A, n - 1, j), 1);
    }
    fmpq_set_si(fmpq_mat_entry(B, n - 1, 0), -(n + 2), 4);
  }
  slow_assert(fmpq_mat_solve(X, A, B), "fatou_coefficients: singular system");
  fmpz_clear(bin);
  fmpq_mat_clear(A);
  fmpq_mat_clear(B);
}

template<class S> S round_fmpq(const fmpq_t q);
template<> double round_fmpq<double>(const fmpq_t q) { return fmpq_get_d(q); }
template<> E2 round_fmpq<E2>(const fmpq_t q) {
  Arb x;
  arb_set_fmpq(x, q, 256);
  const auto r = round_nearest<E2>(x, 256);
  slow_assert(r, "round_fmpq failed");
  return *r;
}

// Petal tests and radii (in double)
constexpr double kAttract = 0.04;   // |w| below which the attracting series is used
constexpr double kRepel = 160;      // Re ζ' ≤ -max(kRepel, 2|Im ζ'|) for the repelling local inverse: |w| ≲ 0.04
inline bool in_attracting(const Cd w) { return cabs(w) < kAttract && std::fabs(w.i) < 0.7 * std::fabs(w.r); }

}  // namespace

template<class S> LavaursModel<S>::LavaursModel(const int terms) {
  const int M = terms + 6;  // The last unknowns of a truncated system are inaccurate: solve longer, keep `terms`
  fmpq_mat_t X;
  fmpq_mat_init(X, M, 1);
  fatou_coefficients(M, X);
  c_log = round_fmpq<S>(fmpq_mat_entry(X, 0, 0));
  for (int j = 1; j <= terms; j++) a.push_back(round_fmpq<S>(fmpq_mat_entry(X, j, 0)));
  fmpq_mat_clear(X);
}

template<class S> void LavaursModel<S>::series(const C w, const bool attracting, C& s, C& d, C& dd) const {
  const Cd wd = to_cd(w);
  const bool pos = attracting ? wd.r > 0 : wd.i > 0;
  const C L = clog_s(pos ? w : -w);
  const C i1 = cdiv(C(1), w), i2 = sqr(i1), i3 = i2 * i1, i4 = sqr(i2);
  const S q(0.25), h(0.5);
  s = q * i2 + q * i1 + c_log * L;
  d = -(h * i3) - q * i2 + c_log * i1;
  dd = S(1.5) * i4 + h * i3 - c_log * i2;
  C p(1), pm(0);
  for (int j = 1; j <= int(a.size()); j++) {
    if (j >= 2) dd = dd + S(double(j * (j - 1))) * a[j - 1] * pm;
    d = d + S(double(j)) * a[j - 1] * p;
    pm = p;
    p = p * w;
    s = s + a[j - 1] * p;
  }
}

template<class S> bool LavaursModel<S>::phi_a(C w, C& s, C& d, C& dd, int& petal, const int max_steps) const {
  C d1(1), d2(0);
  for (int n = 0; n < max_steps; n++) {
    const Cd wd = to_cd(w);
    if (in_attracting(wd)) {
      C ss, sd, sdd;
      series(w, true, ss, sd, sdd);
      s = ss - C(S(0.5 * n));
      d = sd * d1;
      dd = sdd * sqr(d1) + sd * d2;
      petal = (wd.r > 0 ? 1 : -1) * (n % 2 ? -1 : 1);
      return true;
    }
    const C df = twice(w) - C(1);
    d2 = twice(sqr(d1)) + df * d2;
    d1 = df * d1;
    w = sqr(w) - w;
    if (cabs(wd) > 10) return false;
  }
  return false;
}

template<class S> bool LavaursModel<S>::psi(const C zeta, const int petal, C& w, C& d, C& dd) const {
  if (petal < 0) {
    C v, v1, v2;
    if (!psi(zeta - C(S(0.5)), 1, v, v1, v2)) return false;
    const C df = twice(v) - C(1);
    w = sqr(v) - v;
    d = df * v1;
    dd = twice(sqr(v1)) + df * v2;
    return true;
  }
  const Cd zd = to_cd(zeta);
  const double m_real = std::ceil(zd.r + std::max(kRepel, 2 * std::fabs(zd.i)));
  const int64_t m = std::max(int64_t(0), int64_t(m_real));
  if (m > 100000) return false;
  const C zl = zeta - C(S(double(m)));
  // Local inverse in the upper repelling petal: 1/(4w²) ≈ ζ', Newton on the series
  const auto g = std::sqrt(std::complex<double>(-to_cd(zl).r, -to_cd(zl).i));
  C v = from_cd<S>(to_cd(cdiv(Cd(0, 1), Cd(2 * g.real(), 2 * g.imag()))));
  C s, sd, sdd;
  for (int it = 0; it < 60; it++) {
    series(v, false, s, sd, sdd);
    const C step = cdiv(s - zl, sd);
    v = v - step;
    if (cabs(to_cd(step)) < 1e-33 * cabs(to_cd(v))) break;
  }
  series(v, false, s, sd, sdd);
  C d1 = cdiv(C(1), sd);
  C d2 = -(sdd * d1 * sqr(d1));
  for (int64_t i = 0; i < 2 * m; i++) {
    const C df = twice(v) - C(1);
    d2 = twice(sqr(d1)) + df * d2;
    d1 = df * d1;
    v = sqr(v) - v;
  }
  w = v; d = d1; dd = d2;
  return true;
}

namespace {

// (value, ∂w, ∂σ, ∂ww, ∂wσ) of the return map
template<class S> struct Jet {
  typedef Complex<S> C;
  C v, w, s, ww, ws;
  Jet apply(const C f0, const C f1, const C f2) const {
    return Jet{f0, f1 * w, f1 * s, f2 * sqr(w) + f1 * ww, f2 * w * s + f1 * ws};
  }
};

template<class S> bool return_map(const LavaursModel<S>& L, const int r, const int n, const Complex<S> w,
                                  const Complex<S> sigma, Jet<S>& x) {
  typedef Complex<S> C;
  x = Jet<S>{w, C(1), C(0), C(0), C(0)};
  const auto F = [&]() { x = x.apply(sqr(x.v) - x.v, twice(x.v) - C(1), C(2)); };
  F();
  for (int t = 0; t < r; t++) {
    C p0, p1, p2, q0, q1, q2;
    int petal;
    if (!L.phi_a(x.v, p0, p1, p2, petal)) return false;
    x = x.apply(p0, p1, p2);
    x.v = x.v + sigma;
    x.s = x.s + C(1);
    if (!L.psi(x.v, petal, q0, q1, q2)) return false;
    x = x.apply(q0, q1, q2);
  }
  for (int i = 0; i < n; i++) F();
  return true;
}

// Newton on (R_σ(w) - w, R_σ'(w) - μ); returns the last step size, or infinity on failure
template<class S> double newton(const LavaursModel<S>& L, const int r, const int n, Complex<S>& w, Complex<S>& sigma,
                                const Complex<S> mu, const double tol, const int iters = 40) {
  typedef Complex<S> C;
  double last = INFINITY;
  for (int it = 0; it < iters; it++) {
    Jet<S> x;
    if (!return_map(L, r, n, w, sigma, x)) return INFINITY;
    const C F1 = x.v - w, F2 = x.w - mu, a = x.w - C(1), b = x.s, d = x.ww, e = x.ws, det = a * e - b * d;
    const C dw = cdiv(F1 * e - b * F2, det), ds = cdiv(a * F2 - d * F1, det);
    w = w - dw;
    sigma = sigma - ds;
    last = cabs(to_cd(dw)) + cabs(to_cd(ds));
    if (!(last < 1)) return INFINITY;
    if (last < tol) break;
  }
  return last;
}

template<class S> LavaursResult area_t(const int r, const int n, const Cd guess, const int N) {
  typedef Complex<S> C;
  static const LavaursModel<S> L;
  const double tol = sizeof(S) == sizeof(double) ? 1e-15 : 1e-30, accept = sizeof(S) == sizeof(double) ? 1e-11 : 1e-26;
  LavaursResult res{false, Cd(0), E2(0.0), 0, 0};
  C w = from_cd<S>(Cd(0.5, 0)), s = from_cd<S>(guess);
  if (!(newton(L, r, n, w, s, C(0), tol) < accept)) return res;
  const C center = s;
  // Out to the first boundary point radially, then around the circle in small steps
  const auto tw = [N](const int m) { return nearest_twiddle<S>(m, 2 * N); };  // e^{iπm/N}
  const int radial = 8, sub = 4;
  for (int i = 1; i <= radial; i++)
    if (!(newton(L, r, n, w, s, S(double(i) / radial) * tw(1), tol) < accept)) return res;
  vector<C> pts(N);
  for (int j = 0; j < N; j++) {
    if (j) {
      const Cd a = to_cd(tw(2 * j - 1)), b = to_cd(tw(2 * j + 1));
      for (int t = 1; t < sub; t++) {
        Cd m = a + Cd(double(t) / sub * (b.r - a.r), double(t) / sub * (b.i - a.i));
        m = Cd(m.r / cabs(m), m.i / cabs(m));
        if (!(newton(L, r, n, w, s, from_cd<S>(m), 1e-9) < 1e-6)) return res;
      }
    }
    if (!(newton(L, r, n, w, s, tw(2 * j + 1), tol) < accept)) return res;
    pts[j] = s - center;
  }
  // a_k = (1/N) Σ_j σ_j μ_j^-k, μ_j^-k = conj(tw(k(2j+1) mod 2N)); the N/2 subrule uses even j only
  S sum(0.0), half(0.0);
  Cd deriv(0), a1(0);  // σ'(1) = Σ k a_k (the 1/N normalization cancels in the ratio to a_1)
  for (int k = 1; k < N; k++) {
    C ak(0), hk(0);
    for (int j = 0; j < N; j++) {
      const C t = pts[j] * conj(tw(int((int64_t(k) * (2 * j + 1)) % (2 * N))));
      ak = ak + t;
      if (j % 2 == 0 && k < N / 2) hk = hk + t;
    }
    if (k < N / 2) deriv = deriv + Cd(double(k)) * to_cd(ak);  // high k are aliased: keep the resolved half
    if (k == 1) a1 = to_cd(ak);
    sum = sum + S(double(k)) * (sqr(ak.r) + sqr(ak.i));
    if (k < N / 2) half = half + S(double(k)) * (sqr(hk.r) + sqr(hk.i));
  }
  const S pi = nearest_pi<S>();
  const S A = pi * sum / S(double(N) * N), A2 = pi * half / S(double(N / 2) * (N / 2));
  res.ok = true;
  res.center = to_cd(center);
  res.area = E2(A);
  res.conv = double((A - A2) / A);
  res.cusp = cabs(deriv) / cabs(a1);
  return res;
}

}  // namespace

LavaursResult lavaurs_area(const int r, const int n, const Complex<double> guess, const int N) {
  return area_t<E2>(r, n, guess, N);
}

namespace {

const LavaursModel<double>& model_d() { static const LavaursModel<double> L; return L; }

// ζ0 = Φ_a(v) for the critical value v = -1/4 and the petal it enters
void zeta0(Cd& z0, int& petal0) {
  static Cd z; static int p = 0;
  if (!p) {
    Cd d, dd;
    slow_assert(model_d().phi_a(Cd(-0.25, 0), z, d, dd, p), "zeta0 failed");
  }
  z0 = z; petal0 = p;
}

inline char side(const Cd w) { return w.r - 0.5 < 0 ? 'L' : 'R'; }
inline Cd Fd(const Cd w) { return sqr(w) - w; }

}  // namespace

std::string lavaurs_address(const int n, const Complex<double> sigma) {
  Cd z0, x, d, dd; int p0;
  zeta0(z0, p0);
  if (!model_d().psi(z0 + sigma, p0, x, d, dd)) return "?";
  std::string a;
  for (int i = 0; i < n; i++) { a += side(x); x = Fd(x); }
  return a;
}

bool lavaurs_theta(const Complex<double> sigma, Complex<double>& theta, Complex<double>& dtheta,
                   Complex<double>& d2theta, const int r) {
  // p_1 = ζ0 + σ, p_{i+1} = H(p_i) + σ (H = Φ_a ∘ Ψ with each exit through the petal of the previous entry);
  // Θ_r = p_r - ζ0 with its σ-derivatives
  Cd z0; int petal;
  zeta0(z0, petal);
  Cd p = z0 + sigma, dp(1), ddp(0);
  for (int i = 1; i < r; i++) {
    Cd q0, q1, q2, a0, a1, a2;
    if (!model_d().psi(p, petal, q0, q1, q2)) return false;
    if (!model_d().phi_a(q0, a0, a1, a2, petal)) return false;
    const Cd h1 = a1 * q1, h2 = a2 * sqr(q1) + a1 * q2;  // H', H''
    ddp = h2 * sqr(dp) + h1 * ddp;
    dp = h1 * dp + Cd(1);
    p = a0 + sigma;
  }
  theta = p - z0;
  dtheta = dp;
  d2theta = ddp;
  return true;
}

bool lavaurs_island(const Complex<double> center, const Complex<double> target, const int branch,
                    Complex<double>& sigma, const int r) {
  // Path tracking of Θ(σ) = Θ_c + t (target - Θ_c), t: 0 → 1: Euler predictor, Newton corrector, adaptive steps (a step
  // is redone at half size if Newton needs more than a small correction relative to the predicted move), so the lift
  // stays on its branch inside the island
  Cd t0, d, dd;
  if (!lavaurs_theta(center, t0, d, dd, r)) return false;
  const Cd span = target - t0;
  double t = 0, dt = 1e-3;
  // Start a little way out on the branch: the quadratic Θ ≈ Θ_c + δ + Θ''δ²/2
  {
    const double t1 = 1e-6;
    const Cd b = cdiv(Cd(2) * d, dd), c = cdiv(Cd(-2) * Cd(t1) * span, dd);
    const Cd disc = sqr(b) - Cd(4) * c;
    const auto r = std::sqrt(std::complex<double>(disc.r, disc.i));
    const Cd rt(r.real(), r.imag());
    sigma = center + Cd(0.5) * (-b + (branch > 0 ? rt : -rt));
    t = t1;
  }
  const auto correct = [&](Cd& s, const Cd tgt, double& moved) {
    moved = 0;
    for (int it = 0; it < 12; it++) {
      Cd v, dv, ddv;
      if (!lavaurs_theta(s, v, dv, ddv, r)) return false;
      const Cd st = cdiv(v - tgt, dv);
      s = s - st;
      moved += cabs(st);
      if (cabs(v - tgt) < 1e-11 * std::max(1.0, cabs(tgt)) || cabs(st) < 1e-15) return true;
    }
    return false;
  };
  for (int guard = 0; t < 1 && guard < 200000; guard++) {
    const double step = std::min(dt, 1 - t);
    Cd v, dv, ddv;
    if (!lavaurs_theta(sigma, v, dv, ddv, r)) return false;
    const Cd pred = sigma + cdiv(Cd(step) * span, dv);  // Euler: dσ/dt = span / Θ'
    Cd s = pred;
    double moved;
    const bool ok = correct(s, t0 + Cd(t + step) * span, moved);
    const double predicted = cabs(pred - sigma);
    if (!ok || moved > 0.1 * predicted) {
      dt = step / 2;
      if (dt < 1e-12) return false;
      continue;
    }
    sigma = s;
    t += step;
    dt = std::min(step * 1.5, 0.02);
  }
  return t >= 1;
}

bool lavaurs_island_center(Complex<double>& sigma, const int r) {
  // Newton on H'(ζ0 + σ) = 0, H' = Θ' - 1, H'' = Θ''
  for (int it = 0; it < 60; it++) {
    Cd v, d, dd;
    if (!lavaurs_theta(sigma, v, d, dd, r)) return false;
    const Cd st = cdiv(d - Cd(1), dd);
    sigma = sigma - st;
    if (cabs(st) < 1e-14) return true;
  }
  return false;
}

char lavaurs_island_side(const int n_u, const Complex<double> sigma) {
  Cd z0, x, d, dd; int p0;
  zeta0(z0, p0);
  if (!model_d().psi(z0 + sigma, p0, x, d, dd)) return '?';
  for (int i = 0; i < n_u; i++) x = Fd(x);
  return side(x);
}

template struct LavaursModel<double>;
template struct LavaursModel<E2>;

}  // namespace mandelbrot
