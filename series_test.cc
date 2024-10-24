// Series tests

#include "series.h"
#include "expansion.h"
#include "poly.h"
#include "rand.h"
#include "tests.h"
#include <random>
namespace mandelbrot {
namespace {

using std::mt19937;
using std::runtime_error;
using std::uniform_real_distribution;

const int prec = 100, mag_bits = 0;

void poly_rand(Poly& x, Rand& random, const int64_t n) {
  arb_poly_randtest(x, random, n, prec, mag_bits);
  poly_mid(x, x);
}

TEST(construct) {
  Series<double> x(5);
  ASSERT_EQ(x.known(), 0);
  ASSERT_EQ(x.nonzero(), 0);
  ASSERT_EQ(x.limit(), 5);
}

TEST(initializer_list) {
  Series<double> x(3, {5, 7});
  ASSERT_EQ(x.limit(), 3);
  ASSERT_EQ(x.known(), 2);
  ASSERT_EQ(x.nonzero(), 2);
  ASSERT_EQ(x[0], 5);
  ASSERT_EQ(x[1], 7);

  ASSERT_THROW(Series<double> y(2, {1, 2, 3}), runtime_error);
}

TEST(move_construct) {
  Series<double> x(5);
  x.set_scalar(1, 7);
  Series<double> y(move(x));
  ASSERT_EQ(x.known(), 0);
  ASSERT_EQ(x.nonzero(), 0);
  ASSERT_EQ(y.known(), 1);
  ASSERT_EQ(y.nonzero(), 1);
  ASSERT_EQ(y[0], 7);
}

TEST(clear) {
  Series<double> x(5);
  x.set_scalar(1, 7);
  x.clear();
  ASSERT_EQ(x.known(), 0);
  ASSERT_EQ(x.nonzero(), 0);
  ASSERT_EQ(x.limit(), 0);
}

TEST(assign_int) {
  Series<double> x(2);
  x.set_scalar(1, 7);
  ASSERT_EQ(x.known(), 1);
  ASSERT_EQ(x.nonzero(), 1);
  ASSERT_EXACT(x, 7);
  x.set_scalar(2, 3);
  ASSERT_EQ(x.known(), 2);
  ASSERT_EQ(x.nonzero(), 1);
  ASSERT_EXACT(x, 3, 0);
}

TEST(assign_double) {
  Series<double> x(2);
  x.set_scalar(1, 7.5);
  ASSERT_EQ(x.known(), 1);
  ASSERT_EQ(x.nonzero(), 1);
  ASSERT_EXACT(x, 7.5);
  x.set_scalar(2, 3);
  ASSERT_EQ(x.known(), 2);
  ASSERT_EQ(x.nonzero(), 1);
  ASSERT_EXACT(x, 3, 0);
}

TEST(assign_series) {
  Series<double> x(2, {7, 13});
  Series<double> y(2);
  ASSERT_EQ(y.known(), 0);
  ASSERT_EQ(y.nonzero(), 0);
  y = x;
  ASSERT_EXACT(y, 7, 13);
}

TEST(alias) {
  Series<double> x(2), y(2);
  x.set_counts(2, 2);
  y.set_counts(2, 2);
  SeriesView<double> t(x);
  const auto lo = x.low(1);
  const auto hi = x.high(1);
  ASSERT_FALSE(x.alias(y));
  ASSERT_TRUE(x.alias(x));
  ASSERT_TRUE(x.alias(t));
  ASSERT_TRUE(x.alias(lo));
  ASSERT_TRUE(x.alias(hi));
}

TEST(assert_low_near_zero) {
  Series<double> x(3);
  x.set_counts(3, 3);
  x[0] = x[1] = 0;
  x[2] = 1;
  for (int n = 0; n <= 2; n++)
    x.assert_low_near_zero(n);
  ASSERT_THROW(x.assert_low_near_zero(3), runtime_error);
  x.set_counts(1, 1);
  ASSERT_THROW(x.assert_low_near_zero(2), runtime_error);
}

TEST(print) {
  Series<double> x(2);
  ASSERT_EQ(tfm::format("%g", x), "[]");
  x.set_scalar(1, 7);
  ASSERT_EQ(tfm::format("%g", x), "[7]");
  x.set_counts(2, 2); x[1] = 4.5;
  ASSERT_EQ(tfm::format("%g", x), "[7, 4.5]");
}

TEST(set_known) {
  Series<double> x(2);

  // Extending past limit is fine, since we fill in symbolic zeros
  x.set_known(7);
  ASSERT_EQ(x.known(), 7);
  ASSERT_EQ(x.nonzero(), 0);

  // Extend
  x.set_counts(2, 2);
  ASSERT_EQ(x.known(), 2);
  ASSERT_EQ(x.nonzero(), 2);
  x[0] = 7;
  x[1] = 13;
  ASSERT_EXACT(x, 7, 13);

  // Extending more should just change the counts
  x.set_counts(1, 1);
  x.set_known(7);
  ASSERT_EQ(x.known(), 7);
  ASSERT_EQ(x.nonzero(), 1);
  ASSERT_EXACT(x, 7, 0, 0, 0, 0, 0, 0);
  x.set_known(4);
  ASSERT_EQ(x.known(), 4);
  ASSERT_EQ(x.nonzero(), 1);
  ASSERT_EXACT(x, 7, 0, 0, 0);
}

TEST(low_high) {
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

TEST(add_scalar) {
  Series<double> x(2);
  x.set_counts(2, 2);
  x[0] = 0;
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

TEST(add_series) {
  Series<double> x(2), y(2, {3, 5});
  x.set_known(2);
  ASSERT_EXACT(x, 0, 0);
  x += y;
  ASSERT_EXACT(x, 3, 5);
  x += y;
  ASSERT_EXACT(x, 6, 10);
  x -= y;
  ASSERT_EXACT(x, 3, 5);

  // Adding into high part
  x.high_add(1, y.high(1));
  ASSERT_EXACT(x, 3, 10);
  x.high_sub(1, y.low(1));
  ASSERT_EXACT(x, 3, 7);

  // Add/sub with x extension
  for (const int sign : {1, -1}) {
    Series<double> x(2), y(2, {3, 5});
    x.set_counts(2, 1); x[0] = 4; x.data()[1] = 99;
    if (sign > 0) x += y;
    else x -= y;
    ASSERT_EXACT(x, 4 + sign*3, sign*5);
  }

  // Add/sub with y extension
  for (const int sign : {1, -1}) {
    Series<double> x(2), y(2, {3});
    y.set_known(2); y.data()[1] = 99;
    x.set_counts(2, 2); x[0] = 4; x[1] = 2;
    if (sign > 0) x += y;
    else x -= y;
    ASSERT_EXACT(x, 4 + sign*3, 2);
  }

  // High add/sub with x extension
  for (const int sign : {1, -1}) {
    Series<double> x(3), y(2, {3, 5});
    x.set_counts(3, 0);
    x.data()[0] = x.data()[1] = x.data()[2] = 99;
    if (sign > 0) x.high_add(1, y);
    else x.high_sub(1, y);
    ASSERT_EXACT(x, 0, sign*3, sign*5);
  }
}

TEST(neg) {
  Series<double> x(3, {3, 5, 7}), y(3);
  y = neg(x);
  ASSERT_EXACT(y, -3, -5, -7);
}

TEST(mul) {
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
      ASSERT_LT(e, 1e-14) << tfm::format("e = %g\nz = %.3g\naz = %.3g", e, z, az);
    }
    // Aliased versions
    {
      z = x;
      z = mul(z, y);
      const auto e = error(z, az);
      ASSERT_LT(e, 1e-14) << tfm::format("e = %g\nz = %.3g\naz = %.3g", e, z, az);
    } {
      z = y;
      z = mul(x, z);
      const auto e = error(z, az);
      ASSERT_LT(e, 1e-14) << tfm::format("e = %g\nz = %.3g\naz = %.3g", e, z, az);
    }
  }
}

TEST(sqr) {
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
      ASSERT_LT(e, 1e-14) << tfm::format("e = %g\ny = %.3g\nay = %.3g", e, y, ay);
    }
    // Aliased version
    {
      y = x;
      y = sqr(y);
      const auto e = error(y, ay);
      ASSERT_LT(e, 1e-14) << tfm::format("e = %g\ny = %.3g\nay = %.3g", e, y, ay);
    }
  }
}

TEST(newton_steps) {
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
      ASSERT_EQ(fast, slow) << tfm::format("n0 %d, n %d, slow %d, fast %d", n0, n, slow, fast);
    }
  }
}

TEST(inv) {
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
    ASSERT_LT(e, 1e-6) << tfm::format("e = %g\ny = %.3g\nay = %.3g", e, y, ay);
  }
}

TEST(div) {
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
    ASSERT_LT(e, 4e-6) << tfm::format("n = %d, e = %g\n\nx = %.3g\n\ny = %.3g\n\nz = %.3g\n\naz = %.3g",
                                      n, e, approx(ax, n), approx(ay, n), z, az);
  }
}

TEST(mul1p) {
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
      ASSERT_LT(e, 1e-14) << tfm::format("\ne = %g\nx = %.3g\n\ny = %.3g\n\nz = %.3g\n\naz = %.3g",
                                         e, approx(ax, n), approx(ay, n), z, az);
      // Aliasing works
      Series<double> w(n);
      w = x;
      w = mul1p(w, y.low(n-s), s);
      ASSERT_LT(error(z, w), 1e-14);
      w = y;
      w = mul1p(x, w.low(n-s), s);
      ASSERT_LT(error(z, w), 1e-14);
    }
  }
}

TEST(inv1p) {
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
      ASSERT_LT(e, 1e-8) << tfm::format("\ns = %d, e = %g\n\nx = %.3g\n\ny = %.3g\n\nay = %.3g",
                                        s, e, approx(ax, n), y, approx(ay, n));
    }
  }
}

TEST(div1p) {
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
      ASSERT_LT(e, 2e-10) << tfm::format("\nn %d, e = %g\n\nx = %.3g\n\ny = %.3g\n\nz = %.3g\n\naz = %.3g",
                                         n, e, approx(ax, n), approx(ay, n), z, az);
    }
  }
}

TEST(log) {
  Rand random;
  Poly ay, ax;
  for (const int n : {0, 1, 2, 3, 4, 5, 11, 16, 23}) {
    poly_rand(ax, random, n);
    if (n) arb_poly_set_coeff_si(ax, 0, 1);
    arb_poly_log_series(ay, ax, n, prec);
    Series<double> y(n);
    y = log(approx(ax, n));
    const auto e = error(y, ay);
    ASSERT_LT(e, 1e-12) << tfm::format("e = %g\ny = %.3g\nay = %.3g", e, y, ay);
  }
}

TEST(log1p) {
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
      ASSERT_LT(e, 2.3e-10) << tfm::format("\ns = %d, e = %g\n\nx = %.3g\n\ny = %.3g\n\nay = %.3g",
                                           s, e, approx(ax, n), y, approx(ay, n));
    }
  }
}

TEST(derivative_shift) {
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
      ASSERT_LT(e, 1e-14) << tfm::format("\ne = %g\nx = %.3g\n\ny = %.3g\n\nay = %.3g",
                                         e, approx(ax, n), y, approx(ay, n));
    }
  }
}

TEST(integral_shift) {
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
      ASSERT_LT(e, 1e-14) << tfm::format("\ne = %g\nx = %.3g\n\ny = %.3g\n\nay = %.3g",
                                         e, approx(ax, n), y, approx(ay, n));
    }
  }
}

TEST(exp) {
  Rand random;
  Poly ay, ax;
  for (const int n : {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 23}) {
    poly_rand(ax, random, n);
    arb_poly_set_coeff_si(ax, 0, 0);
    arb_poly_exp_series(ay, ax, n, prec);
    Series<double> y(n);
    y = exp(approx(ax, n));
    const auto e = error(y, ay);
    ASSERT_LT(e, 1e-13) << tfm::format("n %d, e = %g\nx = %.3g\ny = %.3g\nay = %.3g", n, e, approx(ax, n), y, ay);
  }
}

TEST(expm1) {
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
        ASSERT_LT(e, 1.2e-13) << tfm::format("\na %d, s %d, e %g\n\nx = %.3g\n\ny = %.3g\n\nay = %.3g",
                                             a, s, e, approx(ax, n), y, approx(ay, n));
      }
    }
  }
}

TEST(log1p_exp) {
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
      ASSERT_LT(e, 4.5e-5) << tfm::format("\nn %d, s %d, e %g\n\nx = %.3g\n\ny = %.3g\n\nay = %.3g",
                                          n, s, e, approx(ax, n), y, approx(ay, n));
    }
  }
}

TEST(ldexp) {
  Series<double> x(3, {3, 5, 7}), y(3);
  y = ldexp(x, 3); ASSERT_EXACT(y, 3<<3, 5<<3, 7<<3);
  y = ldexp(x, -1); ASSERT_EXACT(y, 1.5, 2.5, 3.5);
}

template<class S> void write_series_test() {
  mt19937 mt;
  const int batch_size = 17;
  for (const int n : {0, 1, 2, 5, 32, 1024}) {
    Series<S> x(n);
    x.set_counts(n + 17, n);
    for (int i = 0; i < n; i++) {
      if constexpr (is_same_v<S,double>)
        x[i] = uniform_real_distribution<double>(-1, 1)(mt);
      else
        x[i] = random_expansion<S::n>(mt);
    }
    vector<string> comments{"one", "two"};
    Tmpfile tmp("series-");
    write_series(tmp.path, comments, x.view(), batch_size);
    const auto [cs, y] = read_series<S>(tmp.path);
    for (const auto& c : comments)
      ASSERT_TRUE(std::find(cs.begin(), cs.end(), c) != cs.end());
    ASSERT_EQ(x.known(), y.known());
    ASSERT_EQ(x.nonzero(), y.nonzero());
    for (int64_t i = 0; i < n; i++)
      ASSERT_EQ(x[i], y[i]);
  }
}
TEST(write_series_double) { write_series_test<double>(); }
TEST(write_series_expansion) { write_series_test<Expansion<2>>(); }

}  // namespace
}  // namespace mandelbrot
