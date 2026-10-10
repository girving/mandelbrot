// The Lavaurs model at the root of the cardioid's p/q bulb (double precision)

#include "glavaurs.h"
#include "acb_cc.h"
#include "arb_cc.h"
#include "debug.h"
#include <cmath>
#include <complex>
#include <flint/acb.h>
#include <flint/acb_mat.h>
namespace mandelbrot {

namespace {

typedef Complex<double> Cd;
typedef std::complex<double> SC;
inline SC sc(const Cd z) { return SC(z.r, z.i); }
inline Cd cd(const SC z) { return Cd(z.real(), z.imag()); }
inline double cabs(const Cd z) { return std::hypot(z.r, z.i); }
inline Cd cdiv(const Cd a, const Cd b) { return cd(sc(a) / sc(b)); }
inline double wrap(const double t) { return std::remainder(t, 2 * M_PI); }   // to [-π, π]

constexpr double kRepel = 400;    // Re ζ' ≤ -max(kRepel, 2|Im ζ'|) for the repelling local inverse

struct Jet {
  Cd v, w, s, ww, ws;  // value, ∂w, ∂σ, ∂ww, ∂wσ
  Jet apply(const Cd f0, const Cd f1, const Cd f2) const {
    return Jet{f0, f1 * w, f1 * s, f2 * sqr(w) + f1 * ww, f2 * w * s + f1 * ws};
  }
};

}  // namespace

GeneralLavaurs::GeneralLavaurs(const int p_, const int q_, const int side_, const int N_) : p(p_), q(q_), N(N_), side(side_) {
  lam = cd(std::polar(1.0, 2 * M_PI * p / q));
  v = -(sqr(lam) * Cd(0.25));
  crit = -(lam * Cd(0.5));
  // Unknowns a_j, j = -q..Nu (j ≠ 0), and β; equations: [w^m] of Φ(f(w)) - Φ(w) - 1/q, m = -q+1 .. Nu+1, with
  // (λw)^j (1 + w/λ)^j - w^j and β log(1 + w/λ) (the change of L under f)
  const int Nu = N + 2 * q, prec = 512;
  std::vector<int> idx;
  for (int j = -q; j <= Nu; j++) if (j) idx.push_back(j);
  const int n = idx.size() + 1;
  acb_mat_t M, X, B;
  acb_mat_init(M, n, n);
  acb_mat_init(X, n, 1);
  acb_mat_init(B, n, 1);
  Acb lamb, ilam, t, c;
  acb_set_d_d(lamb, lam.r, lam.i);
  // exact λ: e^{2πi p/q}
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
  a.assign(N + q + 1, Cd(0));
  for (int col = 0; col < int(idx.size()); col++) {
    const int j = idx[col];
    if (j > N) continue;
    const acb_srcptr e = acb_mat_entry(X, col, 0);
    a[j + q] = Cd(arf_get_d(arb_midref(acb_realref(e)), ARF_RND_NEAR), arf_get_d(arb_midref(acb_imagref(e)), ARF_RND_NEAR));
  }
  {
    const acb_srcptr e = acb_mat_entry(X, n - 1, 0);
    beta = Cd(arf_get_d(arb_midref(acb_realref(e)), ARF_RND_NEAR), arf_get_d(arb_midref(acb_imagref(e)), ARF_RND_NEAR));
  }
  acb_mat_clear(M);
  acb_mat_clear(X);
  acb_mat_clear(B);
  A = -cdiv(Cd(1), Cd(q) * a[0]);
  if (std::fabs(A.i) < 1e-12 * cabs(A)) A.i = 0;   // a canonical arg for the petal axes
  // the attracting series is used where |a_{-q}| |w|^{-q} ≥ 625 (0.02 at q = 2), as trustworthy as psi's
  r0 = std::pow(cabs(a[0]) / 625, 1.0 / q);
  Cd s, d, dd;
  slow_assert(phi_a(v, s, d, dd, kv), "the critical value does not enter a petal");
}

Cd GeneralLavaurs::transit_shift(const int k) const {
  const int j = (exit_petal(k) - k) - (exit_petal(kv) - kv);
  return cd(SC(0, 2 * M_PI * j / q) * sc(beta));
}

int GeneralLavaurs::exit_petal(const int k) const {
  const double att = -std::arg(-sc(A)) / q + 2 * M_PI * k / q + side * M_PI / q;
  const double base = -std::arg(sc(A)) / q;
  const int j = int(std::lround((att - base) / (2 * M_PI / q)));
  return ((j % q) + q) % q;
}

void GeneralLavaurs::series(const Cd w, const double ax, Cd& s, Cd& d, Cd& dd) const {
  const SC W = sc(w);
  // log(w^q)/q continued from the axis: a fixed integer branch of the principal log(w^q) would jump where arg(w^q) = π,
  // which for q = 2 is the repelling axis itself
  const SC L(std::log(std::abs(W)), ax + wrap(std::arg(W) - ax));
  SC S = sc(beta) * L, D = sc(beta) / W, DD = -sc(beta) / (W * W);
  SC pw = std::pow(W, -q);   // w^j for j = -q
  for (int j = -q; j <= N; j++, pw *= W) {
    if (!j) continue;
    const SC c = sc(a[j + q]);
    S += c * pw;
    D += c * double(j) * pw / W;
    DD += c * double(j) * double(j - 1) * pw / (W * W);
  }
  s = cd(S); d = cd(D); dd = cd(DD);
}

int GeneralLavaurs::petal(const Cd w, const int kind) const {
  const double t = std::arg(sc(A) * std::pow(sc(w), q) * double(kind));
  if (std::fabs(t) > 0.6 * M_PI) return -1;
  const double base = -std::arg(sc(A) * double(kind)) / q;
  const int k = int(std::lround((std::arg(sc(w)) - base) / (2 * M_PI / q)));
  return ((k % q) + q) % q;
}

double GeneralLavaurs::axis(const int kind, const int k) const {
  return -std::arg(sc(A) * double(kind)) / q + 2 * M_PI * k / q;
}

bool GeneralLavaurs::phi_a(Cd w, Cd& s, Cd& d, Cd& dd, int& pet, const int max_steps) const {
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

bool GeneralLavaurs::psi(const Cd zeta, const int k, Cd& w, Cd& d, Cd& dd) const {
  if (!(cabs(zeta) < 1e6)) return false;
  const int64_t m = std::max(int64_t(0), int64_t(std::ceil(zeta.r + std::max(kRepel, 2 * std::fabs(zeta.i)))));
  if (m > 200000) return false;
  const Cd zl = zeta - Cd(double(m));
  // leading term a_{-q} w^{-q} ≈ ζ': the q-th root in repelling petal k
  const double base = -std::arg(sc(A)) / q + 2 * M_PI * k / q;
  const SC r = std::pow(sc(a[0]) / sc(zl), 1.0 / q);
  SC best = r;
  double bd = 10;
  for (int t = 0; t < q; t++) {
    const SC u = r * std::polar(1.0, 2 * M_PI * t / q);
    const double dlt = std::fabs(wrap(std::arg(u) - base));
    if (dlt < bd) { bd = dlt; best = u; }
  }
  Cd u = cd(best);
  const double br = axis(1, k);
  Cd s, sd, sdd;
  for (int it = 0; it < 60; it++) {
    series(u, br, s, sd, sdd);
    const Cd step = cdiv(s - zl, sd);
    u = u - step;
    if (cabs(step) < 1e-16 * cabs(u)) break;
  }
  series(u, br, s, sd, sdd);
  Cd d1 = cdiv(Cd(1), sd);
  Cd d2 = -(sdd * d1 * sqr(d1));
  for (int64_t i = 0; i < q * m; i++) {
    const Cd df = lam + twice(u);
    d2 = Cd(2) * sqr(d1) + df * d2;
    d1 = df * d1;
    u = lam * u + sqr(u);
    if (cabs(u) > 10) return false;
  }
  w = u; d = d1; dd = d2;
  return true;
}

Cd GeneralLavaurs::zeta0() const {
  Cd s, d, dd;
  int k;
  slow_assert(phi_a(v, s, d, dd, k), "zeta0: the critical value does not enter a petal");
  return s;
}

namespace {

bool return_map(const GeneralLavaurs& L, const int r, const int n, const Cd w, const Cd sigma, Jet& x) {
  x = Jet{w, Cd(1), Cd(0), Cd(0), Cd(0)};
  const auto F = [&]() { x = x.apply(L.lam * x.v + sqr(x.v), L.lam + twice(x.v), Cd(2)); };
  F();
  for (int t = 0; t < r; t++) {
    Cd p0, p1, p2, q0, q1, q2;
    int pet;
    if (!L.phi_a(x.v, p0, p1, p2, pet)) return false;
    x = x.apply(p0, p1, p2);
    x.v = x.v + sigma;
    x.s = x.s + Cd(1);
    if (!L.psi(x.v + L.transit_shift(pet), L.exit_petal(pet), q0, q1, q2)) return false;
    x = x.apply(q0, q1, q2);
  }
  for (int i = 0; i < n; i++) F();
  return true;
}

double newton(const GeneralLavaurs& L, const int r, const int n, Cd& w, Cd& sigma, const Cd mu, const double tol,
              const int iters = 40) {
  double last = INFINITY;
  for (int it = 0; it < iters; it++) {
    Jet x;
    if (!return_map(L, r, n, w, sigma, x)) return INFINITY;
    const Cd F1 = x.v - w, F2 = x.w - mu, a = x.w - Cd(1), b = x.s, d = x.ww, e = x.ws, det = a * e - b * d;
    const Cd dw = cdiv(F1 * e - b * F2, det), ds = cdiv(a * F2 - d * F1, det);
    w = w - dw;
    sigma = sigma - ds;
    last = cabs(dw) + cabs(ds);
    if (!(last < 1)) return INFINITY;
    if (last < tol) break;
  }
  return last;
}

}  // namespace

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

bool GeneralLavaurs::theta(const Cd sigma, const int r, Cd& th, Cd& dth, Cd& d2th, Cd& hprod) const {
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
  return true;
}

}  // namespace mandelbrot
