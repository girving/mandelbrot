// Poly tests

#include "poly.h"
#include "rand.h"
#include "tests.h"
namespace mandelbrot {
namespace {

using std::max;

TEST(log_refine) {
  const int prec = 200, mag_bits = 5;
  const bool verbose = false;
  Rand random;
  Poly x, fast, slow;
  for (int n = 0; n < 32; n++) {
    do {
      arb_poly_randtest(x, random, n, prec, mag_bits);
    } while (n && !arb_is_positive(x[0]));
    // Fast
    poly_log_refine(fast, x, n, prec);
    // Slow
    arb_poly_log_series(slow, x, n, prec);
    // Compare
    if (verbose)
      print("\nn %d\n  x %.3g\n  fast %.3g\n  slow %.3g", n, x, fast, slow);
    ASSERT_TRUE(overlaps(slow, fast));
  }
}

TEST(log1p_shift_refine) {
  const int prec = 200, mag_bits = 5;
  const bool verbose = false;
  Rand random;
  Poly x, t, fast, slow;
  for (int n = 0; n < 32; n++) {
    for (int s = 1; s <= 3; s++) {
      arb_poly_randtest(x, random, n, prec, mag_bits);
      // Fast
      poly_log1p_shift_refine(fast, x, s, n, prec);
      // Slow
      arb_poly_shift_left(t, x, s);
      arb_poly_add_si(t, t, 1, prec);
      arb_poly_log_series(slow, t, n+s, prec);
      slow.assert_low_zero(s);
      slow >>= s;
      // Compare
      if (verbose)
        print("\nn %d, s %d\n  x %.3g\n  fast %.3g\n  slow %.3g", n, s, x, fast, slow);
      ASSERT_TRUE(overlaps(slow, fast));
    }
  }
}

TEST(expm1_shift_refine) {
  const int prec = 200, mag_bits = 5;
  const bool verbose = false;
  Rand random;
  Poly x, fast, slow;
  for (int n = 0; n < 32; n++) {
    for (int s = 1; s <= 3; s++) {
      for (int a = -1; a <= 1; a += 2) {
        arb_poly_randtest(x, random, n, prec, mag_bits);
        // Slow
        arb_poly_shift_left(slow, x, s);
        if (a < 0)
          arb_poly_neg(slow, slow);
        arb_poly_exp_series(slow, slow, n+s, prec);
        slow >>= s;
        if (verbose)
          print("\nn %d, s %d, a %d\n  x %.3g\n  slow %.3g", n, s, a, x, slow);
        // Fast
        poly_expm1_shift_refine(fast, x, a, s, n, prec);
        if (verbose)
          print("  fast %.3g", fast);
        // Compare
        ASSERT_TRUE(overlaps(slow, fast));
      }
    }
  }
}

TEST(log1p_exp_shift_refine) {
  const int prec = 200, mag_bits = 5;
  const bool verbose = false;
  Rand random;
  Poly x, fast, slow;
  for (int n = 0; n < 32; n++) {
    for (int s = 1; s <= 3; s++) {
      arb_poly_randtest(x, random, n, prec, mag_bits);
      // Slow
      arb_poly_exp_series(slow, x, n, prec);
      poly_log1p_shift_refine(slow, slow, s, n, prec);
      if (verbose)
        print("\nn %d, s %d\n  x %.3g\n  slow %.3g", n, s, x, slow);
      // Fast
      poly_log1p_exp_shift_refine(fast, x, s, n, prec);
      if (verbose)
        print("  fast %.3g", fast);
      // Compare
      ASSERT_TRUE(overlaps(slow, fast));
    }
  }
}

TEST(addsub_shift_series) {
  const int prec = 200, mag_bits = 5;
  Rand random;
  Poly f, g, fast0, fast1, fast2, fast3, slow, slow3;
  for (int n = 0; n < 32; n++) {
    for (int s = 0; s <= 3; s++) {
      for (int a = -1; a <= 1; a += 2) {
        arb_poly_randtest(f, random, n, prec, mag_bits);
        arb_poly_randtest(g, random, n, prec, mag_bits);
        // Slow
        const auto S = a > 0 ? arb_poly_add_series : arb_poly_sub_series;
        arb_poly_shift_left(slow, g, s);
        S(slow, f, slow, n, prec);
        arb_poly_shift_left(slow3, f, s);
        S(slow3, f, slow3, n, prec);
        // Fast
        const auto F = a > 0 ? poly_add_shift_series : poly_sub_shift_series;
        F(fast0, f, g, s, n, prec);
        fast1 = f;
        F(fast1, fast1, g, s, n, prec);
        fast2 = g;
        F(fast2, f, fast2, s, n, prec);
        fast3 = f;
        F(fast3, fast3, fast3, s, n, prec);
        // Compare
        ASSERT_TRUE(overlaps(slow, fast0));
        ASSERT_TRUE(overlaps(slow, fast1));
        ASSERT_TRUE(overlaps(slow, fast2));
        ASSERT_TRUE(overlaps(slow3, fast3));
      }
    }
  }
}

}  // namespace
}  // namespace mandelbrot
