// Tests for the Lavaurs model at the 1/2 root

#include "lavaurs.h"
#include "expansion_arith.h"
#include "tests.h"
#include <cmath>
namespace mandelbrot {
namespace {

typedef Expansion<2> E2;
typedef Complex<E2> Ce;

double err(const Ce a, const Ce b) { const Ce d = a - b; return std::sqrt(double(d.r) * double(d.r) + double(d.i) * double(d.i)); }

TEST(series) {
  const LavaursModel<double> L;
  ASSERT_EQ(L.c_log, 11.0 / 8);
  ASSERT_EQ(L.a[0], -5.0 / 16);
  ASSERT_EQ(L.a[1], 75.0 / 64);
  ASSERT_LT(std::fabs(L.a[2] + 149.0 / 192), 1e-16);
}

TEST(fatou) {
  const LavaursModel<E2> L;
  // Φ(F(w)) = Φ(w) + 1/2 on the attracting series region (both points in it: a test of the series itself)
  for (const Ce w : {Ce(E2(0.03), E2(0.004)), Ce(E2(-0.035), E2(0.002))}) {
    const Ce fw = sqr(w) - w;
    Ce s0, d0, dd0, s1, d1, dd1;
    L.series(w, true, s0, d0, dd0);
    L.series(fw, true, s1, d1, dd1);
    ASSERT_LT(err(s1 - s0, Ce(E2(0.5))), 1e-26);
  }
  // Ψ₊(ζ + 1) = F²(Ψ₊(ζ)) with different local points, and Φ(Ψ₊(ζ)) = ζ deep in the petal
  for (const Ce z : {Ce(E2(-3.25), E2(1.5)), Ce(E2(0.75), E2(-2.0)), Ce(E2(-400.5), E2(7.0))}) {
    Ce w0, d0, dd0, w1, d1, dd1;
    ASSERT_TRUE(L.psi(z, 1, w0, d0, dd0));
    ASSERT_TRUE(L.psi(z + Ce(1), 1, w1, d1, dd1));
    const Ce f2 = sqr(sqr(w0) - w0) - (sqr(w0) - w0);
    ASSERT_LT(err(w1, f2), 1e-24 * std::max(1.0, err(f2, Ce(0))));
  }
  Ce w, d, dd, s, sd, sdd;
  const Ce z(E2(-900.0), E2(50.0));
  ASSERT_TRUE(L.psi(z, 1, w, d, dd));
  L.series(w, false, s, sd, sdd);
  ASSERT_LT(err(s, z), 1e-26);
}

TEST(family_constants) {
  // Family constants C = (π²/4) area from the M data: S1 (n = 1), S2 (n = 0), as lim k^4 area by Neville in 1/k
  const double pi2 = M_PI * M_PI / 4;
  const auto S1 = lavaurs_area(1, 1, Complex<double>(0.17958, 0.63460));
  ASSERT_TRUE(S1.ok);
  ASSERT_LT(std::fabs(pi2 * double(S1.area) / 1.724974063498964762e-3 - 1), 1e-12);
  ASSERT_LT(std::fabs(S1.conv), 1e-18);
  const auto S2 = lavaurs_area(1, 0, Complex<double>(-0.28524, 0.95846));
  ASSERT_TRUE(S2.ok);
  ASSERT_LT(std::fabs(pi2 * double(S2.area) / 2.70222138515e-5 - 1), 1e-10);
}

}  // namespace
}  // namespace mandelbrot
