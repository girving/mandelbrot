// Octave tests

#include "octaves.h"
#include "tests.h"
#include <cmath>
#include <random>
namespace mandelbrot {
namespace {

using std::mt19937;
using std::uniform_real_distribution;

TEST(octave_energy) {
  mt19937 rand(7);
  uniform_real_distribution<double> uniform(-1, 1);
  for (int j = 1; j <= 7; j++) {
    const int64_t lo = int64_t(1) << j, n = 2*lo;
    vector<double> f(2*n);
    for (auto& v : f) v = uniform(rand);
    const auto e = octave_energy(f, j);
    ASSERT_EQ(int64_t(e.size()), lo + 1);

    // Direct DFT of a_m, placed at offsets m - 2^j, then folded
    double total = 0;
    for (int64_t m = lo; m < n; m++) total += (m - 1) * std::pow(f[2*m] + f[2*m+1], 2);
    double sum = 0;
    for (int64_t i = 0; i <= lo; i++) {
      double re = 0, im = 0;
      for (int64_t m = lo; m < n; m++) {
        const double a = std::sqrt(double(m - 1)) * (f[2*m] + f[2*m+1]), t = -2 * M_PI * i * (m - lo) / n;
        re += a * std::cos(t);
        im += a * std::sin(t);
      }
      const double want = (i == 0 || i == lo ? 1 : 2) * (re*re + im*im) / n;
      ASSERT_LE(std::abs(e[i] - want), 1e-12 * total) << tfm::format("j %d, i %d: %g != %g", j, i, e[i], want);
      sum += e[i];
    }
    ASSERT_LE(std::abs(sum - total), 1e-13 * total);
  }
}

TEST(fit_line) {
  const vector<double> x = {1, 2, 3, 4, 5}, y = {3.5, 5.5, 7.5, 9.5, 11.5};
  const auto L = fit_line(x, y);
  ASSERT_LE(std::abs(L.a - 1.5), 1e-14);
  ASSERT_LE(std::abs(L.b - 2), 1e-14);
  // Perturbing y by (0, .1, 0, -.1, 0) changes the slope by sum (x - 3) dy / sum (x - 3)^2 = -.2/10
  const vector<double> z = {3.5, 5.6, 7.5, 9.4, 11.5};
  ASSERT_LE(std::abs(fit_line(x, z).b - 1.98), 1e-14);
}

}  // namespace
}  // namespace mandelbrot
