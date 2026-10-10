// Tests of the expansion elementary functions (expansion_math.h) against arb

#include "expansion_math.h"
#include "arb_cc.h"
#include "tests.h"
#include <flint/arb.h>
#include <random>
namespace mandelbrot {
namespace {

using std::mt19937;
using std::uniform_real_distribution;

const int prec = 400;

// A random expansion near d whose lower parts are filled
template<int n> Expansion<n> filled(mt19937& mt, const double d) {
  uniform_real_distribution<double> u(-1, 1);
  Expansion<n> x(d);
  for (int i = 1; i < n; i++) x = x + Expansion<n>(std::ldexp(d * u(mt), -53 * i - 1));
  return x;
}

// Relative error of an expansion against an arb value, in units of 2^{-52 n}
template<int n> double rel_err(const Expansion<n> y, const Arb& truth) {
  Arb e, a;
  arb_sub(e, exact_arb(y), truth, prec);
  arb_abs(e, e);
  arb_abs(a, truth);
  arb_div(e, e, a, prec);
  arb_mul_2exp_si(e, e, 52 * n);
  return arf_get_d(arb_midref(e.x), ARF_RND_NEAR);
}

template<int n> void unary(const char* name, const double lo, const double hi, Expansion<n> (*f)(Expansion<n>),
                           void (*g)(arb_t, const arb_t, slong), const double tol) {
  mt19937 mt(7);
  uniform_real_distribution<double> u(lo, hi);
  double worst = 0;
  for (int i = 0; i < 400; i++) {
    const auto x = filled<n>(mt, u(mt));
    Arb t;
    g(t, exact_arb(x), prec);
    const double e = rel_err(f(x), t);
    worst = std::max(worst, e);
    ASSERT_LE(e, tol) << tfm::format("%s<%d>(%.17g): error %g units", name, n, double(x), e);
  }
}

template<int n> Expansion<n> sqrt_(Expansion<n> x) { return gl_sqrt(x); }
template<int n> Expansion<n> exp_(Expansion<n> x) { return gl_exp(x); }
template<int n> Expansion<n> log_(Expansion<n> x) { return gl_log(x); }
template<int n> Expansion<n> sin_(Expansion<n> x) { Expansion<n> s, c; gl_sincos(x, s, c); return s; }
template<int n> Expansion<n> cos_(Expansion<n> x) { Expansion<n> s, c; gl_sincos(x, s, c); return c; }

template<int n> void all() {
  unary<n>("sqrt", 1e-3, 1e3, sqrt_<n>, arb_sqrt, 16);
  unary<n>("exp", -30, 30, exp_<n>, arb_exp, 64);
  unary<n>("log", 1e-3, 1e3, log_<n>, arb_log, 64);
  unary<n>("sin", -10, 10, sin_<n>, arb_sin, 64);
  unary<n>("cos", -10, 10, cos_<n>, arb_cos, 64);
  // atan2 over all quadrants
  mt19937 mt(11);
  uniform_real_distribution<double> u(-5, 5);
  for (int i = 0; i < 400; i++) {
    const auto y = filled<n>(mt, u(mt)), x = filled<n>(mt, u(mt));
    Arb t;
    arb_atan2(t, exact_arb(y), exact_arb(x), prec);
    const double e = rel_err(gl_atan2(y, x), t);
    ASSERT_LE(e, 64) << tfm::format("atan2<%d>(%.17g, %.17g): error %g units", n, double(y), double(x), e);
  }
  // π itself
  Arb p;
  arb_const_pi(p, prec);
  ASSERT_LE(rel_err(gl_pi<Expansion<n>>(), p), 2);
}

TEST(math2) { all<2>(); }
TEST(math3) { all<3>(); }

}  // namespace
}  // namespace mandelbrot
