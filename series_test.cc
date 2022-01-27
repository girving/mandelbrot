// Series tests

#include "series.h"
#include "poly.h"
#include "rand.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
namespace mandelbrot {
namespace {

using std::numeric_limits;
using std::runtime_error;
using testing::ElementsAre;

Series<double> approx(const Poly& x, const int64_t n) {
  slow_assert(n >= x.length());
  Series<double> y(n);
  y.set_terms(n);
  for (int64_t i = 0; i < n; i++)
    y[i] = arf_get_d(arb_midref(x[i]), ARF_RND_NEAR);
  return y;
}

double error(const Series<const double>& x, const Series<const double>& y, const bool relative = false) {
  if (x.terms() != y.terms())
    return numeric_limits<double>::infinity();
  double e = 0;
  for (int64_t i = 0; i < x.terms(); i++) {
    double d = abs(x[i] - y[i]);
    if (relative)
      d /= max(1., abs(y[i]));
    e = max(e, d);
  }
  return e;
}
double error(const Series<const double>& x, initializer_list<double>&& ys, const bool relative = false) {
  return error(x, Series<const double>(move(ys)), relative);
}
double error(const Series<const double>& x, const Poly& y, const bool relative = false) {
  if (x.terms() < y.length())
    return numeric_limits<double>::infinity();
  return error(x, approx(y, x.terms()), relative);
}

const int prec = 100, mag_bits = 0;

void poly_rand(Poly& x, Rand& random, const int64_t n) {
  arb_poly_randtest(x, random, n, prec, mag_bits);
  poly_mid(x, x);
}

#define ASSERT_EXACT(x, ...) ASSERT_THAT(x, ElementsAre(__VA_ARGS__))
#define ASSERT_CLOSE(x, ...) { \
  const auto e = error(x, {__VA_ARGS__}); \
  ASSERT_LT(e, 3e-14) << format("e %g, x %g", e, x); }

TEST(series, construct) {
  Series<double> x(5);
  ASSERT_EQ(x.terms(), 0);
  ASSERT_EQ(x.limit(), 5);
}

TEST(series, initializer_list) {
  Series<double> x(3, {5, 7});
  ASSERT_EQ(x.limit(), 3);
  ASSERT_EQ(x.terms(), 2);
  ASSERT_EQ(x[0], 5);
  ASSERT_EQ(x[1], 7);

  ASSERT_THROW(Series<double> y(2, {1, 2, 3}), runtime_error);
}

TEST(series, move_construct) {
  Series<double> x(5);
  x.extend(1);
  x = 7;
  Series<double> y(move(x));
  ASSERT_EQ(x.terms(), 0);
  ASSERT_EQ(y.terms(), 1);
  ASSERT_EQ(y[0], 7);
}

TEST(series, clear) {
  Series<double> x(5);
  x.extend(1);
  x = 7;
  x.clear();
  ASSERT_EQ(x.terms(), 0);
  ASSERT_EQ(x.limit(), 0);
}

TEST(series, assign_int) {
  Series<double> x(2);
  x.extend(1);
  x = 7;
  ASSERT_EXACT(x, 7);
  x.extend(2);
  x[1] = 7;
  x = 3;
  ASSERT_EXACT(x, 3, 0);
}

TEST(series, assign_double) {
  Series<double> x(2);
  x.extend(1);
  x = 7.5;
  ASSERT_EQ(x.terms(), 1);
  ASSERT_EQ(x[0], 7.5);
  x.extend(2);
  x[1] = 7;
  x = 3;
  ASSERT_EXACT(x, 3, 0);
}

TEST(series, assign_series) {
  Series<double> x(2, {7, 13});
  Series<double> y(2);
  ASSERT_EQ(y.terms(), 0); 
  y = x;
  ASSERT_EXACT(y, 7, 13);
}

TEST(series, alias) {
  Series<double> x(2), y(2);
  x.extend(2);
  y.extend(2);
  Series<double> t(x);
  const auto lo = x.low(1);
  const auto hi = x.high(1);
  ASSERT_FALSE(x.alias(y));
  ASSERT_TRUE(x.alias(x));
  ASSERT_TRUE(x.alias(t));
  ASSERT_TRUE(x.alias(lo));
  ASSERT_TRUE(x.alias(hi));
}

TEST(series, assert_low_near_zero) {
  Series<double> x(3); 
  x.extend(3);
  x[2] = 1;
  for (int n = 0; n <= 2; n++)
    x.assert_low_near_zero(n);
  ASSERT_THROW(x.assert_low_near_zero(3), runtime_error);
  x.truncate(1);
  ASSERT_THROW(x.assert_low_near_zero(2), runtime_error);
}

TEST(series, print) {
  Series<double> x(2);
  ASSERT_EQ(format("%g", x), "[]");
  x.extend(1); x[0] = 7;
  ASSERT_EQ(format("%g", x), "[7]");
  x.extend(2); x[1] = 4.5;
  ASSERT_EQ(format("%g", x), "[7, 4.5]");
}

TEST(series, truncate_extend) {
  Series<double> x(2); 

  // Truncate
  ASSERT_EQ(x.terms(), 0);
  x.truncate(-1);
  ASSERT_EQ(x.terms(), 0);
  x.truncate(3);
  ASSERT_EQ(x.terms(), 0);
  
  // Extend
  ASSERT_THROW(x.extend(3), runtime_error);
  x.extend(2);
  ASSERT_EQ(x.terms(), 2);
  x[0] = 7;
  x[1] = 13;
  ASSERT_EXACT(x, 7, 13);

  // Truncate and extend should fill in zeros
  x.truncate(1);
  x.extend(2);
  ASSERT_EXACT(x, 7, 0);
}

TEST(series, low_high) {
  Series<double> x(3, {3, 7, 13});
  ASSERT_EXACT(x, 3, 7, 13);
  ASSERT_EXACT(x.low(-1));
  ASSERT_EXACT(x.low(0));
  ASSERT_EXACT(x.low(1), 3);
  ASSERT_EXACT(x.low(2), 3, 7);
  ASSERT_EXACT(x.low(3), 3, 7, 13);
  ASSERT_EXACT(x.low(4), 3, 7, 13);
  ASSERT_EXACT(x.high(-1), 3, 7, 13);
  ASSERT_EXACT(x.high(0), 3, 7, 13);
  ASSERT_EXACT(x.high(1), 7, 13);
  ASSERT_EXACT(x.high(2), 13);
  ASSERT_EXACT(x.high(3));
  ASSERT_EXACT(x.high(4));

  // Verify that they alias
  const auto lo = x.low(1);
  const auto hi = x.high(1);
  x[0] = 23;
  x[2] = 37;
  ASSERT_EQ(lo[0], 23);
  ASSERT_EQ(hi[1], 37);
}

TEST(series, add_scalar) {
  Series<double> x(2);
  x.extend(2);
  x[1] = 3;
  x += 7;
  ASSERT_EXACT(x, 7, 3);
  x -= 1;
  ASSERT_EXACT(x, 6, 3);
  x += 1.5;
  ASSERT_EXACT(x, 7.5, 3);
  x -= 0.25;
  ASSERT_EXACT(x, 7.25, 3);
}

TEST(series, add_series) {
  Series<double> x(2), y(2, {3, 5});
  x.extend(2);
  ASSERT_EXACT(x, 0, 0);
  x += y;
  ASSERT_EXACT(x, 3, 5);
  x += y;
  ASSERT_EXACT(x, 6, 10);
  x -= y;
  ASSERT_EXACT(x, 3, 5);

  // Adding into high part
  x.high(1) += x.high(1);
  ASSERT_EXACT(x, 3, 10);
  x.high(1) -= x.low(1);
  ASSERT_EXACT(x, 3, 7);
}

TEST(series, mul) {
  // Small n
  Series<double> x(3, {3, 5, 7}), y(3, {2, 3, 4});
  Series<double> small, z(3);
  { Series<double> small; ASSERT_THROW(small = mul(x, y), runtime_error); }
  { Series<double> z(3); z = mul(x, y); ASSERT_CLOSE(z, 6, 19, 41); }
  { Series<double> z(2); z = mul(x, y.low(2)); ASSERT_CLOSE(z, 6, 19); }
  { Series<double> z(2); z = mul(x.low(2), y); ASSERT_CLOSE(z, 6, 19); }
  { Series<double> z(2); z = mul(x.low(2), y.low(2)); ASSERT_CLOSE(z, 6, 19); }
  { Series<double> z(1); z = mul(x.low(1), y); ASSERT_EXACT(z, 6); }
  { Series<double> z(1); z = mul(x, y.low(1)); ASSERT_EXACT(z, 6); }
  { Series<double> z; z = mul(x.low(0), y); ASSERT_EXACT(z); }
  { Series<double> z; z = mul(x, y.low(0)); ASSERT_EXACT(z); }

  // Aliasing should work
  { Series<double> z(3); z = x; z = mul(z, z); ASSERT_CLOSE(z, 9, 30, 67); }

  // Large n
  Rand random;
  Poly az, ax, ay;
  for (const int n : {5, 11, 16, 23}) {
    poly_rand(ax, random, n);
    poly_rand(ay, random, n);
    arb_poly_mullow(az, ax, ay, n, prec);
    Series<double> z(n);
    const auto x = approx(ax, n);
    const auto y = approx(ay, n);
    {
      z = mul(x, y);
      const auto e = error(z, az);
      ASSERT_LT(e, 1e-14) << format("e = %g\nz = %.3g\naz = %.3g", e, z, az);
    }
    // Aliased versions
    {
      z = x;
      z = mul(z, y);
      const auto e = error(z, az);
      ASSERT_LT(e, 1e-14) << format("e = %g\nz = %.3g\naz = %.3g", e, z, az);
    } {
      z = y;
      z = mul(x, z);
      const auto e = error(z, az);
      ASSERT_LT(e, 1e-14) << format("e = %g\nz = %.3g\naz = %.3g", e, z, az);
    }
  }
}

TEST(series, sqr) {
  // Small n
  Series<double> x(3, {2, 3, 4});
  { Series<double> small; ASSERT_THROW(small = sqr(x), runtime_error); }
  { Series<double> y; y = sqr(x.low(0)); ASSERT_CLOSE(y); }
  { Series<double> y(1); y = sqr(x.low(1)); ASSERT_CLOSE(y, 4); }
  { Series<double> y(2); y = sqr(x.low(2)); ASSERT_CLOSE(y, 4, 12); }
  { Series<double> y(3); y = sqr(x); ASSERT_CLOSE(y, 4, 12, 25); }

  // Aliasing should work
  { Series<double> z(3); z = x; z = sqr(z); ASSERT_CLOSE(z, 4, 12, 25); }

  // Large n
  Rand random;
  Poly ay, ax;
  for (const int n : {5, 11, 16, 23}) {
    poly_rand(ax, random, n);
    arb_poly_mullow(ay, ax, ax, n, prec);
    Series<double> y(n);
    const auto x = approx(ax, n);
    {
      y = sqr(x);
      const auto e = error(y, ay);
      ASSERT_LT(e, 1e-14) << format("e = %g\ny = %.3g\nay = %.3g", e, y, ay);
    }
    // Aliased version
    {
      y = x;
      y = sqr(y);
      const auto e = error(y, ay);
      ASSERT_LT(e, 1e-14) << format("e = %g\ny = %.3g\nay = %.3g", e, y, ay);
    }
  }
}

TEST(series, newton_steps) {
  const auto slow_steps = [](const int64_t n0, int64_t n) {
    int steps = 0;
    while (n0 < n) {
      n = (n+1)/2;
      steps++;
    }
    return steps;
  };
  for (int n0 = 1; n0 < 20; n0++) {
    for (int n = n0; n < 100; n++) {
      const auto slow = slow_steps(n0, n);
      const auto fast = newton_steps(n0, n);
      ASSERT_EQ(fast, slow) << format("n0 %d, n %d, slow %d, fast %d", n0, n, slow, fast);
    }
  }
}

TEST(series, inv) {
  // Small n
  Series<double> x(3, {2, 3, 4});
  ASSERT_THROW(x = inv(x), runtime_error);
  { Series<double> small; ASSERT_THROW(small = inv(x), runtime_error); }
  { Series<double> y; y = inv(x.low(0)); ASSERT_CLOSE(y); }
  { Series<double> y(1); y = inv(x.low(1)); ASSERT_CLOSE(y, .5); }
  { Series<double> y(2); y = inv(x.low(2)); ASSERT_CLOSE(y, .5, -.75); }
  { Series<double> y(3); y = inv(x); ASSERT_CLOSE(y, .5, -.75, .125); }

  // Large n
  Rand random;
  Poly ay, ax;
  for (const int n : {5, 11, 16, 23}) {
    poly_rand(ax, random, n);
    arb_poly_inv_series(ay, ax, n, prec);
    Series<double> y(n);
    y = inv(approx(ax, n));
    const auto e = error(y, ay);
    ASSERT_LT(e, 1e-6) << format("e = %g\ny = %.3g\nay = %.3g", e, y, ay);
  }
}

TEST(series, div) {
  Rand random;
  Poly az, ax, ay;
  for (const int n : {0, 1, 2, 3, 5, 11, 16, 23}) {
    poly_rand(ax, random, n);
    do {
      poly_rand(ay, random, n);
    } while (n && !arb_is_nonzero(ay[0]));
    arb_poly_div_series(az, ax, ay, n, prec);
    Series<double> z(n);
    z = div(approx(ax, n), approx(ay, n));
    const auto e = error(z, az);
    ASSERT_LT(e, 4e-6) << format("n = %d, e = %g\n\nx = %.3g\n\ny = %.3g\n\nz = %.3g\n\naz = %.3g",
                                 n, e, approx(ax, n), approx(ay, n), z, az);
  }
}

TEST(series, mul1p) {
  Rand random;
  Poly az, ax, ay, t;
  for (const int s : {1, 2, 3}) {
    for (const int n : {0, 1, 2, 3, 5, 11, 16, 23}) {
      poly_rand(ax, random, n);
      poly_rand(ay, random, n);
      arb_poly_shift_left(t, ay, s);
      arb_poly_add_si(t, t, 1, prec);
      arb_poly_mullow(az, ax, t, n, prec);
      const auto x = approx(ax, n);
      const auto y = approx(ay, n);
      Series<double> z(n);
      z = mul1p(x, y.low(n-s), s);
      const auto e = error(z, az);
      ASSERT_LT(e, 1e-14) << format("\ne = %g\nx = %.3g\n\ny = %.3g\n\nz = %.3g\n\naz = %.3g",
                                    e, approx(ax, n), approx(ay, n), z, az);
      // Aliasing works with the second argument works
      Series<double> w(n);
      if (n) {
        w = x;
        ASSERT_THROW(w = mul1p(w, y.low(n-s), s), runtime_error);
      }
      w = y;
      w = mul1p(x, w.low(n-s), s);
      ASSERT_LT(error(z, w), 1e-14);
    }
  }
}

TEST(series, inv1p) {
  Rand random;
  Poly ay, ax, t;
  for (const int s : {1, 2, 3}) {
    for (const int n : {0, 1, 2, 3, 5, 11, 16, 23}) {
      poly_rand(ax, random, n);
      arb_poly_shift_left(t, ax, s);
      arb_poly_add_si(t, t, 1, prec);
      arb_poly_inv_series(ay, t, n+s, prec);
      ay >>= s;
      Series<double> y(n);
      y = inv1p(approx(ax, n), s);
      const auto e = error(y, ay);
      ASSERT_LT(e, 1e-8) << format("\ns = %d, e = %g\n\nx = %.3g\n\ny = %.3g\n\nay = %.3g",
                                   s, e, approx(ax, n), y, approx(ay, n));
    }
  }
}

TEST(series, div1p) {
  Rand random;
  Poly az, ax, ay, t;
  for (const int s : {1, 2, 3}) {
    for (const int n : {0, 1, 2, 3, 5, 11, 16, 23}) {
      poly_rand(ax, random, n);
      poly_rand(ay, random, n);
      arb_poly_shift_left(t, ay, s);
      arb_poly_add_si(t, t, 1, prec);
      arb_poly_div_series(az, ax, t, n, prec);
      Series<double> z(n);
      z = div1p(approx(ax, n), approx(ay, n).low(n-s), s);
      const auto e = error(z, az);
      ASSERT_LT(e, 2e-10) << format("\nn %d, e = %g\n\nx = %.3g\n\ny = %.3g\n\nz = %.3g\n\naz = %.3g",
                                    n, e, approx(ax, n), approx(ay, n), z, az);
    }
  }
}

TEST(series, log) {
  Rand random;
  Poly ay, ax;
  for (const int n : {0, 1, 2, 3, 4, 5, 11, 16, 23}) {
    poly_rand(ax, random, n);
    if (n) arb_poly_set_coeff_si(ax, 0, 1);
    arb_poly_log_series(ay, ax, n, prec);
    Series<double> y(n);
    y = log(approx(ax, n));
    const auto e = error(y, ay);
    ASSERT_LT(e, 1e-12) << format("e = %g\ny = %.3g\nay = %.3g", e, y, ay);
  }
}

TEST(series, log1p) {
  Rand random;
  Poly ay, ax, t;
  for (const int s : {1, 2, 3}) {
    for (const int n : {0, 1, 2, 3, 5, 11, 16, 23}) {
      poly_rand(ax, random, n);
      arb_poly_shift_left(t, ax, s);
      arb_poly_log1p_series(ay, t, n+s, prec);
      ay >>= s;
      Series<double> y(n);
      y = log1p(approx(ax, n), s);
      const auto e = error(y, ay);
      ASSERT_LT(e, 1e-9) << format("\ns = %d, e = %g\n\nx = %.3g\n\ny = %.3g\n\nay = %.3g",
                                   s, e, approx(ax, n), y, approx(ay, n));
    }
  }
}

TEST(series, derivative_shift) {
  Rand random;
  Poly ay, ax;
  for (const int s : {0, 1, 2, 3}) {
    for (const int n : {0, 1, 2, 3, 5, 11, 16, 23}) {
      poly_rand(ax, random, n);
      arb_poly_shift_left(ay, ax, s);
      arb_poly_derivative(ay, ay, prec);
      if (s)
        arb_poly_shift_right(ay, ay, s-1);
      else
        arb_poly_shift_left(ay, ay, 1);
      Series<double> y(n);
      y = derivative_shift(approx(ax, n), s);
      const auto e = error(y, ay);
      ASSERT_LT(e, 1e-14) << format("\ne = %g\nx = %.3g\n\ny = %.3g\n\nay = %.3g",
                                    e, approx(ax, n), y, approx(ay, n));
    }
  }
}

TEST(series, integral_shift) {
  Rand random;
  Poly ay, ax;
  for (const int s : {0, 1, 2, 3}) {
    for (const int n : {0, 1, 2, 3, 5, 11, 16, 23}) {
      poly_rand(ax, random, n);
      if (s)
        arb_poly_shift_left(ay, ax, s-1);
      else
        arb_poly_shift_right(ay, ax, 1);
      arb_poly_integral(ay, ay, prec);
      arb_poly_shift_right(ay, ay, s);
      Series<double> y(n);
      y = integral_shift(approx(ax, n), s);
      const auto e = error(y, ay);
      ASSERT_LT(e, 1e-14) << format("\ne = %g\nx = %.3g\n\ny = %.3g\n\nay = %.3g",
                                    e, approx(ax, n), y, approx(ay, n));
    }
  }
}

TEST(series, exp) {
  Rand random;
  Poly ay, ax;
  for (const int n : {0, 1, 2, 3, 4, 5, 11, 16, 23}) {
    poly_rand(ax, random, n);
    arb_poly_set_coeff_si(ax, 0, 0);
    arb_poly_exp_series(ay, ax, n, prec);
    Series<double> y(n);
    y = exp(approx(ax, n));
    const auto e = error(y, ay);
    ASSERT_LT(e, 1e-13) << format("e = %g\ny = %.3g\nay = %.3g", e, y, ay);
  }
}

TEST(series, expm1) {
  Rand random;
  Poly ay, ax, t;
  for (const int s : {1, 2, 3}) {
    for (const int a : {-1, 1}) {
      for (const int n : {0, 1, 2, 3, 5, 11, 16, 23}) {
        poly_rand(ax, random, n);
        arb_poly_shift_left(t, ax, s);
        if (a < 0) arb_poly_neg(t, t);
        arb_poly_exp_series(ay, t, n+s, prec);
        ay >>= s;
        Series<double> y(n);
        y = expm1(approx(ax, n), a, s);
        const auto e = error(y, ay);
        ASSERT_LT(e, 1.2e-13) << format("\na %d, s %d, e %g\n\nx = %.3g\n\ny = %.3g\n\nay = %.3g",
                                        a, s, e, approx(ax, n), y, approx(ay, n));
      }
    }
  }
}

TEST(series, log1p_exp) {
  Rand random;
  Poly ay, ax, t;
  for (const int s : {1, 2, 3}) {
    for (const int n : {0, 1, 2, 3, 5, 11, 16, 23}) {
      poly_rand(ax, random, n);
      if (n) arb_poly_set_coeff_si(ax, 0, 0);
      arb_poly_exp_series(t, ax, n+s, prec);
      arb_poly_shift_left(t, t, s);
      arb_poly_log1p_series(ay, t, n+s, prec);
      ay >>= s;
      Series<double> y(n);
      y = log1p_exp(approx(ax, n), s);
      const auto e = error(y, ay, true);
      ASSERT_LT(e, 5.6e-5) << format("\nn %d, s %d, e %g\n\nx = %.3g\n\ny = %.3g\n\nay = %.3g",
                                     n, s, e, approx(ax, n), y, approx(ay, n));
    }
  }
}

TEST(series, ldexp) {
  Series<double> x(3, {3, 5, 7}), y(3);
  y = ldexp(x, 3); ASSERT_EXACT(y, 3<<3, 5<<3, 7<<3);
  y = ldexp(x, -1); ASSERT_EXACT(y, 1.5, 2.5, 3.5);
}

}  // namespace
}  // namespace mandelbrot
