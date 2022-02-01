// Expansion tests

#include "expansion.h"
#include "arb_cc.h"
#include "arith.h"
#include "mag_cc.h"
#include "nearest.h"
#include "noncopyable.h"
#include "print.h"
#include "rand.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <fenv.h>
#include <functional>
#include <random>
namespace mandelbrot {
namespace {

using std::bernoulli_distribution;
using std::function;
using std::get;
using std::max;
using std::mt19937;
using std::mt19937_64;
using std::nextafter;
using std::numeric_limits;
using std::tuple;
using std::uniform_int_distribution;
using std::uniform_real_distribution;
using std::vector;

int exponent(const double x) {
  int e;
  frexp(x, &e);
  return e;
}

// Goldberg's definition, according to Collange et al.
double ulp(const double x) {
  constexpr int p = numeric_limits<double>::digits;
  static_assert(p == 52 + 1);
  return ldexp(double(1), exponent(x) - p);
}

// Verify Figure 2.6 of Muller et al., Handbook of floating point arithmetic
//   https://doc.lagout.org/science/0_Computer Science/3_Theory/Handbook of Floating Point Arithmetic.pdf
TEST(expansion, ulp) {
  const int p = 52 + 1;
  const double u = pow(2, -p);
  ASSERT_EQ(ulp(0.7), u);
  ASSERT_EQ(ulp(1), 2*u);
  ASSERT_EQ(ulp(1.2), 2*u);
}

// Is ulp(x) >= |y|?
bool ulp_valid(const double x, const double y) {
  if (!x || !y) return true;
  return ulp(x) >= abs(y);
}

template<int n> bool ulp_valid(const Expansion<n>& e) {
  for (int i = 0; i+1 < n; i++)
    if (!ulp_valid(e.x[i], e.x[i+1]))
      return false;
  return true;
}

template<int n> int ulp_slop(const Expansion<n>& e) {
  int slop = 0;
  for (int i = 0; i+1 < n; i++) {
    const double x = e.x[i];
    const double y = e.x[i+1];
    if (x && y)
      slop = max(slop, exponent(y) - exponent(ulp(x)));
  }
  return slop;
}

template<int n> Expansion<n> random_expansion_with_exponent(mt19937& mt, int e) {
  Expansion<n> a;
  for (int i = 0; i < n; i++) {
    a.x[i] = ldexp(uniform_real_distribution<double>(-1, 1)(mt), e);
    e = exponent(a.x[i]) - 52 - 1;
  }
  return a;
}

template<int n> Expansion<n> random_expansion(mt19937& mt) {
  const int e = uniform_int_distribution<int>(-20, 20)(mt);
  return random_expansion_with_exponent<n>(mt, e);
}

template<int n> Expansion<n> random_expansion_near(mt19937& mt, const Expansion<n> x) {
  const int e = exponent(x.x[0]) - uniform_int_distribution<int>(1, 52*n)(mt);
  return (bernoulli_distribution(0.5)(mt) ? x : -x) + random_expansion_with_exponent<n>(mt, e);
}

// Returns an upper bound for |x|
template<int n> double bound(const Expansion<n> x) {
  const auto inf = numeric_limits<double>::infinity();
  double b = 0;
  for (const auto a : x.x)
    b = nextafter(b + abs(a), inf);  // Make sure we round up without bothering with rounding mode
  return b;
}

template<int n> void convert_test() {
  typedef Expansion<n> E;
  const bool verbose = false;
  mt19937 mt(7);
  Rand rand;

  Arb x, e;
  for (int i = 0; i < 1024; i++) {
    const int prec = 2*100*n;
    const int rand_prec = uniform_int_distribution<int>(2, 2*60*n)(mt);
    arb_randtest_exact(x, rand, rand_prec, 7);
    const double ax = bound(x);
    const E y = round_near<E>(arb_midref(x.x));
    if (verbose)
      print("i %d, rp %d, x %g, y %s", i, rand_prec, x, y.span());
    ASSERT_TRUE(ulp_valid(y));
    const Arb z = y.arb(prec);
    arb_sub(e, x, z, prec);
    const double ae = bound(e);
    ASSERT_LE(ae, ldexp(ax, -52*n)) << format("n %d, i %d, ax %g", n, i, ax);
  }
}

template<int n,class I> void convert_int_test() {
  typedef Expansion<n> E;
  mt19937_64 mt(7);

  Arb x, e;
  for (int i = 0; i < 1024; i++) {
    const I a = I(mt());
    const E x(a);
    ASSERT_TRUE(ulp_valid(x))
        << format("n %d, i %d, a %d", n, i, a)
        << format("\nx %.30", x)
        << format("\nx %s (ulp %g)", x.span(), ulp(x.x[0]));
    const Arf y = exact_arf(x);
    ASSERT_TRUE(arf_equal_si(y, a));
  }
}

template<int n> void random_expansion_test() {
  const bool verbose = false;
  typedef Expansion<n> E;
  mt19937 mt(7);
  vector<tuple<string,function<bool(E)>>> anys = {
    {"hi_pos", [](E x) { return bound(x) > 1e3 && x.x[0] > 0; }},
    {"hi_neg", [](E x) { return bound(x) > 1e3 && x.x[0] < 0; }},
    {"lo_pos", [](E x) { return bound(x) < 1e-3 && x.x[0] > 0; }},
    {"lo_pos", [](E x) { return bound(x) < 1e-3 && x.x[0] < 0; }},
  };
  vector<bool> found(anys.size());
  for (int i = 0; i < 64; i++) {
    const E x = random_expansion<n>(mt);
    ASSERT_TRUE(ulp_valid(x)) << x.span();
    for (int i = 0; i < int(anys.size()); i++)
      if (get<1>(anys[i])(x))
        found[i] = true;
    if (verbose)
      print(x.span());
  }
  for (int i = 0; i < int(anys.size()); i++)
    ASSERT_TRUE(found[i]) << format("%s (i %d)", get<0>(anys[i]), i);
}

template<int n,class Op,class ArbOp> void exact_unary_test(const Op& op, const ArbOp& arb_op) {
  typedef Expansion<n> E;
  mt19937 mt(7);

  Arb correct;
  for (int i = 0; i < 1024; i++) {
    const int prec = 400*n;
    const E x = random_expansion<n>(mt);
    const E z = op(x);
    arb_op(correct, x.arb(prec));
    const auto az = z.arb(prec);
    ASSERT_TRUE(arb_is_exact(correct));
    ASSERT_TRUE(arb_equal(correct, az));
  }
}

template<int n> void neg_test() { exact_unary_test<n>(std::negate<void>(), arb_neg); }
template<int n> void abs_test() { exact_unary_test<n>([](auto x) { return abs(x); }, arb_abs); }
template<int n> void ldexp_test() {
  for (const int e : {-10, -5, -1, 0, 1, 5, 10})
    exact_unary_test<n>(
        [e](auto x) { return ldexp(x, e); },
        [e](Arb& y, const Arb& x) { arb_mul_2exp_si(y, x, e); });
}

template<int n> void inv_test() {
  typedef Expansion<n> E;
  mt19937 mt(7);
  const int slop = 2;

  Arb correct, error;
  for (int i = 0; i < 1024; i++) {
    const int prec = 400*n;
    const E x = random_expansion<n>(mt);
    const E z = inv(x);
    arb_inv(correct, x.arb(prec), prec);
    arb_sub(error, correct, z.arb(prec), prec);
    const auto want = ldexp(bound(correct), -52*n + slop);
    ASSERT_LE(bound(error), want);
  }
}

template<int n,class Op,class ArbOp> void
binary_test(const string& name, const Op& op, const ArbOp& arb_op, const bool cancellation,
            const int uslop, const int eslop) {
  typedef Expansion<n> E;
  mt19937 mt(7);

  Arb correct, error;
  for (int i = 0; i < (1<<15); i++) {
    const int prec = 400*n;
    const bool near = cancellation && bernoulli_distribution(.25)(mt);
    const E x = random_expansion<n>(mt);
    const E y = near ? random_expansion_near(mt, x) : random_expansion<n>(mt);
    const E z = op(x, y);
    arb_op(correct, x.arb(prec), y.arb(prec), prec);
    arb_sub(error, z.arb(prec), correct, prec);
    const auto want = ldexp(bound(correct), -52*n + eslop);
    ASSERT_LE(ulp_slop(z), uslop)
        << format("op %s, n %d, i %d, near %d, error %g > %g", name, n, i, near, bound(error), want)
        << format("\n\nx %.100g\ny %.100g\nz %.100g", x, y, z);
    ASSERT_LE(bound(error), want)
        << format("op %s, n %d, i %d, near %d, error %g > %g", name, n, i, near, bound(error), want)
        << format("\n\nx %.100g\ny %.100g\nz %.100g", x, y, z)
        << format("\n\nx %.17g (valid %d)\ny %.17g (valid %d)\nz %.17g (valid %d)",
                  x.span(), ulp_valid(x), y.span(), ulp_valid(y), z.span(), ulp_valid(z))
        << format("\n\nax %s\nay %s\naz %s", safe(x), safe(y), safe(z));
  }
}

template<int n> void add_test() {
  binary_test<n>("+", std::plus<void>(), arb_add, true, n == 2 ? 0 : 11, 0);
}
template<int n> void sub_test() {
  binary_test<n>("-", std::minus<void>(), arb_sub, true, n == 2 ? 0 : 11, 0);
}
template<int n> void mul_test() {
  binary_test<n>("*", std::multiplies<void>(), arb_mul, false, n == 2 ? 1 : 14, n == 2 ? 2 : 1);
}
template<int n> void div_test() {
  binary_test<n>("/", std::divides<void>(), arb_div, false, 0, 2);
}

template<int n> void equal_test() {
  typedef Expansion<n> E;
  mt19937 mt(7);

  for (int i = 0; i < (1<<15); i++) {
    double pieces[2*n];
    for (int i = 0; i < 2*n; i++)
      pieces[i] = ldexp(double(uniform_int_distribution<int>(-32,32)(mt)), -6*i);

    // Equality should hold between first + second halves and odds + evens
    E x, y;
    for (int i = 0; i < n; i++) {
      x.x[0] += pieces[i];
      x.x[1] += pieces[n+i];
      y.x[0] += pieces[2*i];
      y.x[1] += pieces[2*i+1];
    }
    ASSERT_EQ(x, y);

    // Equality should not hold if we tweak anything a little bit
    const double inf = numeric_limits<double>::infinity();
    for (int i = 0; i < n; i++) {
      for (const double to : {inf, -inf}) {
        E tx = x, ty = y;
        tx.x[i] = nextafter(x.x[i], to);
        ty.x[i] = nextafter(y.x[i], to);
        ASSERT_NE(tx, x);
        ASSERT_NE(tx, y);
        ASSERT_NE(ty, x);
        ASSERT_NE(ty, y);
      }
    }
  }
}

#define TESTS(n) \
  TEST(expansion, random_expansion##n) { random_expansion_test<n>(); } \
  TEST(expansion, convert##n) { convert_test<n>(); } \
  TEST(expansion, convert_int##n##_32) { convert_int_test<n,int32_t>(); } \
  TEST(expansion, convert_int##n##_64) { convert_int_test<n,int64_t>(); } \
  TEST(expansion, neg##n) { neg_test<n>(); } \
  TEST(expansion, abs##n) { abs_test<n>(); } \
  TEST(expansion, ldexp##n) { ldexp_test<n>(); } \
  TEST(expansion, add##n) { add_test<n>(); } \
  TEST(expansion, sub##n) { sub_test<n>(); } \
  TEST(expansion, mul##n) { mul_test<n>(); } \
  TEST(expansion, equal##n) { equal_test<n>(); } \
  TEST(expansion, inv##n) { inv_test<n>(); } \
  TEST(expansion, div##n) { div_test<n>(); }
TESTS(2)
TESTS(3)

}  // namespace
}  // namespace mandelbrot
