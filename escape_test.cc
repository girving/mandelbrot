// Escape-time tests

#include "escape.h"
#include "tests.h"
#include <cmath>
#include <complex>
#include <random>
namespace mandelbrot {
namespace {

using std::mt19937;
using std::uniform_real_distribution;

// Reference log2 of the Green's function in long double, escaping at a much larger radius
long double slow_log2g(const long double x, const long double y) {
  std::complex<long double> c(x, y), z = c;
  for (int n = 1; n < 100000; n++) {
    if (std::abs(z) > 1e100L) return std::log2(std::log(std::abs(z))) - (n - 1);
    z = z * z + c;
  }
  return -INFINITY;
}

TEST(cardioid_or_disk) {
  mt19937 rand(7);
  uniform_real_distribution<double> ux(-2, 0.5), uy(-1.2, 1.2);
  int inside = 0;
  for (int i = 0; i < 20000; i++) {
    const double x = ux(rand), y = uy(rand);
    if (!in_cardioid_or_disk(x, y)) continue;
    inside++;
    // Points in the cardioid or disk never escape
    std::complex<double> c(x, y), z = c;
    for (int n = 0; n < 2000 && std::abs(z) <= 2; n++) z = z * z + c;
    ASSERT_LE(std::abs(z), 2) << tfm::format("c = %g + %gi", x, y);
  }
  ASSERT_LT(4000, inside);  // The two components are most of M's area
  // Boundary points: slightly outside the cardioid, and slightly inside
  for (int i = 0; i < 100; i++) {
    const double t = 2 * M_PI * (i + 0.5) / 100;
    const std::complex<double> w(std::cos(t), std::sin(t));
    const auto b = w / 2.0 - w * w / 4.0;
    const auto n = std::complex<double>(0, 1) * (w / 2.0 - w * w / 2.0);  // derivative in t
    const auto out = std::complex<double>(0, -1) * n / std::abs(n);        // outward normal
    const auto o = b + 1e-6 * out, in = b - 1e-6 * out;
    if (std::abs(b + 1.0) < 0.26) continue;  // Near the period 2 disk
    ASSERT_FALSE(in_cardioid_or_disk(o.real(), o.imag())) << t;
    ASSERT_TRUE(in_cardioid_or_disk(in.real(), in.imag())) << t;
  }
}

TEST(green) {
  // Escaping points: g matches a long double reference with a larger escape radius.  Rounding along the
  // orbit acts like perturbing c by ~1e-16, which moves g by ~1e-16 / dist(c, M) relatively, so the
  // tolerance is loose for slowly escaping points.
  mt19937 rand(11);
  uniform_real_distribution<double> ux(-2.5, 1), uy(-1.5, 1.5);
  int checked = 0;
  for (int i = 0; i < 3000; i++) {
    const double x = ux(rand), y = uy(rand);
    const auto e = escape(x, y, 10000);
    // Very slow orbits are chaotic enough that double and long double rounding diverge; see below_precision
    if (e.steps < 0 || e.steps > 300) continue;
    const long double l = slow_log2g(x, y);
    // Relative error in g is |Δ log2 g| · ln 2
    ASSERT_LE(std::abs(e.log2g - double(l)) * M_LN2, 1e-6) << tfm::format("c = %g + %gi: %g vs %g", x, y, e.log2g, double(l));
    checked++;
  }
  ASSERT_LT(1000, checked);
  // Known values: g(-2) = 0 (in M); large c has g ≈ log|c|
  ASSERT_TRUE(escape(-2, 0, 1000).steps < 0 || escape(-2, 0, 1000).g() < 1e-200);
  const auto big = escape(100, 0, 100);
  ASSERT_LE(std::abs(big.g() - std::log(100.0)), 1e-2);  // g(c) = log|c| + ½ log(1 + 1/c) + ...
}

TEST(below) {
  // Slowly escaping real parameters near the cusp: g ≈ exp(-π ln 2 / √ε) roughly, and `below` agrees with g
  for (const double eps : {1e-2, 1e-3, 1e-4}) {
    const auto e = escape(0.25 + eps, 0, 1 << 24);
    ASSERT_LT(0, e.steps);
    ASSERT_LT(0.5 * M_PI / std::sqrt(eps), double(e.steps));  // Passage takes ≈ π / √ε steps
    const int k = int(-e.log2g);
    ASSERT_TRUE(below(e, k - 1) && !below(e, k + 1));  // g < 2^-(k-1) but not < 2^-(k+1)
  }
  // Interior points (a period 3 component, not caught by the cardioid/disk test) are below every threshold
  const auto e = escape(-0.1225611668766536, 0.7448617666197442, 1 << 20);
  ASSERT_TRUE(e.steps < 0 && below(e, 1000000));
}

TEST(period) {
  // Centers of known components report their period
  const struct { double x, y; int p; } cs[] = {
      {0, 0, 1}, {-1, 0, 2}, {-0.1225611668766536, 0.7448617666197442, 3}, {-1.754877666246693, 0, 3},
      {-1.310702641336833, 0, 4}, {0.2822713907669139, 0.5300606175785253, 4}, {-1.985424253054205, 0, 5}};
  for (const auto& c : cs) {
    const auto e = escape(c.x + 1e-8, c.y + 1e-8, 1 << 16);  // Near, not at, the center
    ASSERT_EQ(e.steps, -1);
    ASSERT_EQ(e.period, c.p) << tfm::format("c = %g + %gi", c.x, c.y);
  }
}

// Reference classification without the Newton shortcut: plain iteration with Brent cycle detection
Escape slow_escape(const double x, const double y, const int64_t max_iter) {
  double zx = x, zy = y, cx = x, cy = y;
  int64_t next_check = 16;
  for (int64_t n = 1; n <= max_iter; n++) {
    const double r2 = zx * zx + zy * zy;
    if (r2 > std::ldexp(1.0, 64)) return {n, std::log2(0.5 * std::log(r2)) - double(n - 1), 0, n};
    const double t = zx * zx - zy * zy + x;
    zy = 2 * zx * zy + y;
    zx = t;
    const double dx = zx - cx, dy = zy - cy;
    if (dx * dx + dy * dy < 1e-26) return {-1, -INFINITY, 0, n};
    if (n == next_check) { cx = zx; cy = zy; next_check *= 2; }
  }
  return {-1, -INFINITY, 0, max_iter};
}

TEST(newton_interior) {
  // The Newton certificate never contradicts plain iteration: every escaping point is still classified the
  // same way, and points certified interior do not escape within the slow run
  mt19937 rand(3);
  uniform_real_distribution<double> ux(-2, 0.5), uy(0, 1.2);
  int faster = 0;
  for (int i = 0; i < 200000; i++) {
    const double x = ux(rand), y = uy(rand);
    const auto e = escape(x, y, 1 << 16), s = slow_escape(x, y, 1 << 16);
    ASSERT_EQ(e.steps, s.steps) << tfm::format("c = %.17g + %.17gi", x, y);
    if (e.steps > 0) ASSERT_EQ(e.log2g, s.log2g);
    faster += e.iters < s.iters;
  }
  ASSERT_LT(10000, faster);
}

}  // namespace
}  // namespace mandelbrot
