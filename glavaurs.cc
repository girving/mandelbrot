// The Lavaurs model at the root of the cardioid's p/q bulb (double precision)

#include "glavaurs.h"
#include "acb_cc.h"
#include "arb_cc.h"
#include "arf_cc.h"
#include "debug.h"
#include "expansion_math.h"
#include <cmath>
#include <complex>
#include <flint/acb.h>
#include <flint/acb_mat.h>
namespace mandelbrot {

namespace {

typedef Complex<double> Cd;
inline Cd cd(const std::complex<double> z) { return Cd(z.real(), z.imag()); }
inline double cabs(const Cd z) { return std::hypot(z.r, z.i); }

}  // namespace

namespace {

// The Fatou coefficients a_j (j = -q..N, index j + q; a_0 = 0), β and λ in acb at 512 bits.  Unknowns a_j, j = -q..Nu
// (j ≠ 0), and β; equations: [w^m] of Φ(f(w)) - Φ(w) - 1/q, m = -q+1 .. Nu+1, with (λw)^j (1 + w/λ)^j - w^j and
// β log(1 + w/λ) (the change of L under f)
void fatou_acb(const int p, const int q, const int N, std::vector<Acb>& a, Acb& beta, Acb& lamb) {
  const int Nu = N + 2 * q, prec = 512;
  std::vector<int> idx;
  for (int j = -q; j <= Nu; j++) if (j) idx.push_back(j);
  const int n = idx.size() + 1;
  acb_mat_t M, X, B;
  acb_mat_init(M, n, n);
  acb_mat_init(X, n, 1);
  acb_mat_init(B, n, 1);
  Acb ilam, t, c;
  {
    Arb th;
    arb_const_pi(th, prec);
    arb_mul_si(th, th, 2 * p, prec);
    arb_div_si(th, th, q, prec);
    Arb cs, sn;
    arb_sin_cos(sn, cs, th, prec);
    acb_set_arb_arb(lamb, cs, sn);
  }
  acb_inv(ilam, lamb, prec);
  const int m0 = -q + 1;
  for (int col = 0; col < int(idx.size()); col++) {
    const int j = idx[col];
    // λ^j Σ_k binom(j, k) λ^{-k} w^{j+k}
    acb_pow_si(c, lamb, j, prec);
    for (int k = 0; j + k <= Nu + 1; k++) {
      if (k > 0) { acb_mul_si(c, c, j - k + 1, prec); acb_div_si(c, c, k, prec); acb_mul(c, c, ilam, prec); }
      const int m = j + k;
      if (m >= m0) acb_add(acb_mat_entry(M, m - m0, col), acb_mat_entry(M, m - m0, col), c, prec);
    }
    if (j >= m0) acb_sub_si(acb_mat_entry(M, j - m0, col), acb_mat_entry(M, j - m0, col), 1, prec);
  }
  // β: log(1 + w/λ) = Σ_{m≥1} (-1)^{m+1} w^m / (m λ^m)
  acb_one(c);
  for (int m = 1; m <= Nu + 1; m++) {
    acb_mul(c, c, ilam, prec);
    acb_div_si(t, c, m, prec);
    if (m % 2 == 0) acb_neg(t, t);
    acb_set(acb_mat_entry(M, m - m0, n - 1), t);
  }
  acb_set_si(acb_mat_entry(B, 0 - m0, 0), 1);
  acb_div_si(acb_mat_entry(B, 0 - m0, 0), acb_mat_entry(B, 0 - m0, 0), q, prec);
  slow_assert(acb_mat_solve(X, M, B, prec), "GeneralLavaurs: singular Fatou system");
  a.clear();
  for (int j = -q; j <= N; j++) a.emplace_back();
  for (int col = 0; col < int(idx.size()); col++) {
    const int j = idx[col];
    if (j <= N) acb_set(a[j + q], acb_mat_entry(X, col, 0));
  }
  acb_set(beta, acb_mat_entry(X, n - 1, 0));
  acb_mat_clear(M);
  acb_mat_clear(X);
  acb_mat_clear(B);
}

// arb midpoint → S: the nearest double, or an expansion by exact successive remainders
template<class S> struct FromArf;
template<> struct FromArf<double> {
  static double get(const arf_t x) { return arf_get_d(x, ARF_RND_NEAR); }
};
template<int n> struct FromArf<Expansion<n>> {
  static Expansion<n> get(const arf_t x) {
    Arf r, d;
    arf_set(r, x);
    Expansion<n> y;
    for (int i = 0; i < n; i++) {
      y.x[i] = arf_get_d(r, ARF_RND_NEAR);
      arf_set_d(d, y.x[i]);
      arf_sub(r, r, d, ARF_PREC_EXACT, ARF_RND_NEAR);
    }
    return y;
  }
};
template<class S> Complex<S> from_acb(const acb_t z) {
  return Complex<S>(FromArf<S>::get(arb_midref(acb_realref(z))), FromArf<S>::get(arb_midref(acb_imagref(z))));
}

}  // namespace

template<class S> GLCoreT<S> gl_make_core(const int p, const int q, const int side, const int N) {
  typedef Complex<S> C;
  slow_assert(N + q + 1 <= GLCoreT<S>::kMaxA, "gl_make_core: too many coefficients");
  std::vector<Acb> a;
  Acb beta, lamb;
  fatou_acb(p, q, N, a, beta, lamb);
  GLCoreT<S> core;
  core.p = p; core.q = q; core.side = side; core.N = N;
  for (int j = 0; j < N + q + 1; j++) core.a[j] = from_acb<S>(a[j]);
  core.beta = from_acb<S>(beta);
  core.lam = from_acb<S>(lamb);
  core.v = -(sqr(core.lam) * C(S(0.25)));
  core.crit = -(core.lam * C(S(0.5)));
  core.A = -glcore::cinv(C(S(double(q))) * core.a[0]);
  if (std::fabs(double(core.A.i)) < 1e-12 * std::hypot(double(core.A.r), double(core.A.i))) core.A.i = S(0.0);   // canonical arg
  core.argA = gl_atan2(core.A.i, core.A.r);
  core.argmA = gl_atan2(-core.A.i, -core.A.r);
  // the attracting series is used where |a_{-q}| |w|^{-q} ≥ 625 (0.02 at q = 2), as trustworthy as psi's
  core.r0 = std::pow(std::hypot(double(core.a[0].r), double(core.a[0].i)) / 625, 1.0 / q);
  core.tau = C(S(0.0), twice(gl_pi<S>()) / S(double(q))) * core.beta;
  core.kv = 0;
  C s, d, dd;
  int kv;
  slow_assert(core.phi_a(core.v, s, d, dd, kv), "the critical value does not enter a petal");
  core.kv = kv;
  return core;
}
template GLCoreT<double> gl_make_core(int, int, int, int);
template GLCoreT<Expansion<2>> gl_make_core(int, int, int, int);
template GLCoreT<Expansion<3>> gl_make_core(int, int, int, int);

GeneralLavaurs::GeneralLavaurs(const int p_, const int q_, const int side_, const int N_) : p(p_), q(q_), N(N_), side(side_) {
  core = gl_make_core<double>(p, q, side, N);
  lam = core.lam; A = core.A; v = core.v; crit = core.crit; beta = core.beta; r0 = core.r0; kv = core.kv;
  a.assign(core.a, core.a + N + q + 1);
}

Cd GeneralLavaurs::transit_shift(const int k) const { return core.transit_shift(k); }
int GeneralLavaurs::exit_petal(const int k) const { return core.exit_petal(k); }
void GeneralLavaurs::series(const Cd w, const double ax, Cd& s, Cd& d, Cd& dd) const { core.series(w, ax, s, d, dd); }
int GeneralLavaurs::petal(const Cd w, const int kind) const { return core.petal(w, kind); }
double GeneralLavaurs::axis(const int kind, const int k) const { return core.axis(kind, k); }
bool GeneralLavaurs::phi_a(Cd w, Cd& s, Cd& d, Cd& dd, int& pet, const int max_steps) const {
  return core.phi_a(w, s, d, dd, pet, max_steps);
}
bool GeneralLavaurs::psi(const Cd zeta, const int k, Cd& w, Cd& d, Cd& dd) const { return core.psi(zeta, k, w, d, dd); }

Cd GeneralLavaurs::zeta0() const {
  Cd s, d, dd;
  int k;
  slow_assert(phi_a(v, s, d, dd, k), "zeta0: the critical value does not enter a petal");
  return s;
}

namespace {
double newton(const GeneralLavaurs& L, const int r, const int n, Cd& w, Cd& sigma, const Cd mu, const double tol,
              const int iters = 40) {
  return L.core.newton(r, n, w, sigma, mu, tol, iters);
}
}  // namespace

bool GeneralLavaurs::return_map(const int r, const int n, const Cd w, const Cd sigma, Cd& v, Cd& dw, Cd& dww) const {
  GLJet x;
  if (!core.return_map(r, n, w, sigma, x)) return false;
  v = x.v; dw = x.w; dww = x.ww;
  return true;
}

bool GeneralLavaurs::center(const int r, const int n, Cd& sigma) const {
  Cd w = crit;
  return newton(*this, r, n, w, sigma, Cd(0), 1e-14) < 1e-10;
}

bool GeneralLavaurs::multiplier_point(const int r, const int n, const Cd cen, const Cd mu, Cd& sigma) const {
  Cd w = crit, s = cen;
  const int K = 32;
  for (int k = 1; k <= K; k++)
    if (!(newton(*this, r, n, w, s, Cd(mu.r * k / K, mu.i * k / K), 1e-14) < 1e-9)) return false;
  sigma = s;
  return true;
}

bool GeneralLavaurs::area(const int r, const int n, const Cd guess, Cd& cen, double& ar, double& conv, double& cusp,
                          Cd& a1out, const int Nb) const {
  Cd w = crit, s = guess;
  if (!(newton(*this, r, n, w, s, Cd(0), 1e-14) < 1e-10)) return false;
  cen = s;
  const auto tw = [Nb](const int m) { return cd(std::polar(1.0, M_PI * m / Nb)); };
  const int radial = 8, sub = 4;
  // boundary points accepted at the roundoff floor (Ψ runs q m ~ 400 q iterations): ~1e-8
  const double acc = 1e-8;
  for (int i = 1; i <= radial; i++)
    if (!(newton(*this, r, n, w, s, Cd(double(i) / radial) * tw(1), 1e-14) < acc)) return false;
  std::vector<Cd> pts(Nb);
  for (int j = 0; j < Nb; j++) {
    if (j) {
      const Cd x0 = tw(2 * j - 1), x1 = tw(2 * j + 1);
      for (int t = 1; t < sub; t++) {
        Cd m = x0 + Cd(double(t) / sub) * (x1 - x0);
        m = Cd(m.r / cabs(m), m.i / cabs(m));
        if (!(newton(*this, r, n, w, s, m, 1e-9) < 1e-6)) return false;
      }
    }
    if (!(newton(*this, r, n, w, s, tw(2 * j + 1), 1e-14) < acc)) return false;
    pts[j] = s - cen;
  }
  double sum = 0, half = 0;
  Cd deriv(0), a1(0);
  for (int k = 1; k < Nb; k++) {
    Cd ak(0), hk(0);
    for (int j = 0; j < Nb; j++) {
      const Cd t = pts[j] * conj(tw(int((int64_t(k) * (2 * j + 1)) % (2 * Nb))));
      ak = ak + t;
      if (j % 2 == 0 && k < Nb / 2) hk = hk + t;
    }
    if (k < Nb / 2) deriv = deriv + Cd(double(k)) * ak;
    if (k == 1) a1 = ak;
    sum += k * (sqr(ak.r) + sqr(ak.i));
    if (k < Nb / 2) half += k * (sqr(hk.r) + sqr(hk.i));
  }
  const double A1 = M_PI * sum / (double(Nb) * Nb), A2 = M_PI * half / (double(Nb / 2) * (Nb / 2));
  ar = A1;
  conv = (A1 - A2) / A1;
  cusp = cabs(deriv) / cabs(a1);
  a1out = Cd(a1.r / Nb, a1.i / Nb);
  return true;
}

bool GeneralLavaurs::horn(const Cd p, const int pet, Cd& h, Cd& dh, Cd& ddh, int& pet2) const {
  return core.horn(p, pet, h, dh, ddh, pet2);
}

bool GeneralLavaurs::theta(const Cd sigma, const int r, Cd& th, Cd& dth, Cd& d2th, Cd& hprod, int* pet_out) const {
  return core.theta(sigma, r, th, dth, d2th, hprod, pet_out);
}


// The area of a component at precision S: Newton for the center, then boundary points σ(μ) at μ = e^{iπ(2j+1)/Nb} by
// continuation (radially out, then around with sub intermediate steps), area = π Σ_k k |a_k|² from the Fourier
// coefficients of σ(μ) - center; conv = the relative change from the half grid, cusp = |σ'(1)|/|a_1|
template<class S> bool gl_area(const GLCoreT<S>& L, const int r, const int n, const Complex<S> guess, Complex<S>& cen,
                               S& ar, double& conv, double& cusp, const int Nb, const double tol, const double acc) {
  typedef Complex<S> C;
  C w = L.crit, s = guess;
  if (!(L.newton(r, n, w, s, C(0), tol, 60) < acc)) return false;
  cen = s;
  const auto tw = [Nb](const int m) { return glcore::cpolar(S(1.0), gl_pi<S>() * S(double(m)) / S(double(Nb))); };
  const int radial = 8, sub = 4;
  for (int i = 1; i <= radial; i++)
    if (!(L.newton(r, n, w, s, C(S(double(i) / radial)) * tw(1), tol, 60) < acc)) return false;
  std::vector<C> pts(Nb);
  for (int j = 0; j < Nb; j++) {
    if (j) {
      const C x0 = tw(2 * j - 1), x1 = tw(2 * j + 1);
      for (int t = 1; t < sub; t++) {
        C m = x0 + C(S(double(t) / sub)) * (x1 - x0);
        m = C(S(1.0) / glcore::cabs(m)) * m;
        if (!(L.newton(r, n, w, s, m, 1e-9, 60) < 1e-6)) return false;
      }
    }
    if (!(L.newton(r, n, w, s, tw(2 * j + 1), tol, 60) < acc)) return false;
    pts[j] = s - cen;
  }
  S sum(0.0), half_(0.0);
  C deriv(0), a1(0);
  for (int k = 1; k < Nb; k++) {
    C ak(0), hk(0);
    for (int j = 0; j < Nb; j++) {
      const C tt = tw(int((int64_t(k) * (2 * j + 1)) % (2 * Nb)));
      const C t = pts[j] * C(tt.r, -tt.i);
      ak = ak + t;
      if (j % 2 == 0 && k < Nb / 2) hk = hk + t;
    }
    if (k < Nb / 2) deriv = deriv + C(S(double(k))) * ak;
    if (k == 1) a1 = ak;
    sum = sum + S(double(k)) * (ak.r * ak.r + ak.i * ak.i);
    if (k < Nb / 2) half_ = half_ + S(double(k)) * (hk.r * hk.r + hk.i * hk.i);
  }
  const S A1 = gl_pi<S>() * sum / (S(double(Nb)) * S(double(Nb)));
  const S A2 = gl_pi<S>() * half_ / (S(double(Nb / 2)) * S(double(Nb / 2)));
  ar = A1;
  conv = double((A1 - A2) / A1);
  cusp = double(glcore::cabs(deriv) / glcore::cabs(a1));
  return true;
}
template bool gl_area(const GLCoreT<Expansion<2>>&, int, int, Complex<Expansion<2>>, Complex<Expansion<2>>&,
                      Expansion<2>&, double&, double&, int, double, double);
template bool gl_area(const GLCoreT<Expansion<3>>&, int, int, Complex<Expansion<3>>, Complex<Expansion<3>>&,
                      Expansion<3>&, double&, double&, int, double, double);

}  // namespace mandelbrot
