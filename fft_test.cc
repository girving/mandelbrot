// FFT tests

#include "fft.h"
#include "acb_dft.h"
#include "debug.h"
#include "print.h"
#include "gtest/gtest.h"
#include <cmath>
#include <complex>
#include <random>

namespace std {
template<class T> ostream& operator<<(ostream& out, const vector<T>& x) {
  out << '[';
  for (size_t i = 0; i < x.size(); i++) {
    if (i) out << ", ";
    out << x[i];
  }
  return out << ']';
}
}  // namespace std

namespace mandelbrot {
namespace {

using std::abs;
using std::hypot;
using std::max;
using std::mt19937;
using std::uniform_int_distribution;
using std::uniform_real_distribution;
using std::vector;

double abs(const Complex<double> z) {
  return hypot(z.r, z.i);
}

TEST(fft, fft) {
  mt19937 rand(7);
  uniform_real_distribution<double> uniform(-1, 1);
  for (int p = -1; p <= 15; p++) {
    for (const bool full : {true, false}) {
      const double tol = 2.1e-14 * max(p, 0);
      const int64_t n = p < 0 ? 0 : int64_t(1) << p;
      const int64_t xn = full ? n : uniform_int_distribution<int64_t>(0, n)(rand);
      vector<double> x(xn);
      for (int64_t i = 0; i < xn; i++)
        x[i] = uniform(rand);

      // Accurate forward FFT using arb
      const int prec = 200;
      acb_ptr v = _acb_vec_init(n);
      for (int64_t i = 0; i < xn; i++)
        acb_set_d(v + i, x[i]);
      acb_dft(v, v, n, prec);
      vector<Complex<double>> sy(n);
      const auto mid = [](arb_t x) { return arf_get_d(arb_midref(x), ARF_RND_NEAR); };
      for (int64_t i = 0; i < n; i++)
        sy[i] = Complex<double>(mid(acb_realref(v + i)), mid(acb_imagref(v + i)));
      _acb_vec_clear(v, n);

      // Forward FFT
      vector<Complex<double>> y(n);
      fft(y.data(), x.data(), n, xn);
      for (int64_t i = 0; i < n; i++) {
        const auto e = abs(y[i]-sy[i]);
        ASSERT_LE(e, tol)
            << format("p %d, n %d, i %d, e %g%s", p, n, i, e,
                      n > 4 ? "" : format(":\nx = %g\ny = %g\nsy = %g", x, y, sy));
      }

      // Inverse FFT
      vector<double> z(xn);
      ifft(z.data(), y.data(), n, xn);
      for (int64_t i = 0; i < xn; i++)
        z[i] /= n;
      for (int64_t i = 0; i < xn; i++) {
        const auto e = abs(x[i]-z[i]);
        ASSERT_LE(e, tol) << format("p %d, n %d, i %d, e %g:\nx = %g\nz = %g\ny = %g",
                                    p, n, i, e, x, z, sy);
      }
    }
  }
}

}  // namespace
}  // namespace mandelbrot
