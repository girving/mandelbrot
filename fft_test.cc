// FFT tests

#include "fft.h"
#include "acb_dft.h"
#include "arb_cc.h"
#include "debug.h"
#include "print.h"
#include "tests.h"
#include <cmath>
#include <complex>
#include <random>
namespace mandelbrot {
namespace {

using std::hypot;
using std::max;
using std::mt19937;
using std::uniform_int_distribution;
using std::uniform_real_distribution;
using std::vector;

TEST(fft) {
  mt19937 rand(7);
  uniform_real_distribution<double> uniform(-1, 1);
  for (int p = -1; p <= 15; p++) {
    for (const bool full : {true, false}) {
      const double tol = 3e-14 * max(p, 0);
      const int64_t n = p < 0 ? 0 : int64_t(1) << p;
      const int64_t xn = full ? 2*n : uniform_int_distribution<int64_t>(0, 2*n)(rand);
      vector<double> x(xn);
      for (int64_t i = 0; i < xn; i++)
        x[i] = uniform(rand);
      const auto safe_x = [&x,xn](const int i) { return i < xn ? x[i] : 0; };

      // Accurate forward FFT using arb
      const int prec = 200;
      acb_ptr v = _acb_vec_init(n);
      for (int64_t i = 0; i < n; i++)
        acb_set_d_d(v + i, safe_x(2*i), safe_x(2*i+1));
      acb_dft(v, v, n, prec);
      vector<Complex<double>> sy(n);
      const auto mid = [](arb_t x) { return arf_get_d(arb_midref(x), ARF_RND_NEAR); };
      for (int64_t i = 0; i < n; i++)
        sy[i] = Complex<double>(mid(acb_realref(v + i)), mid(acb_imagref(v + i)));
      _acb_vec_clear(v, n);

      // Forward FFT
      vector<Complex<double>> y(n);
      fft<double>(y, x);
      for (int64_t i = 0; i < n; i++) {
        const auto e = abs(y[i]-sy[i]);
        ASSERT_LE(e, tol)
            << format("p %d, n %d, i %d, e %g%s", p, n, i, e,
                      n > 4 ? "" : format(":\nx = %g\ny = %g\nsy = %g", x, y, sy));
      }

      // Inverse FFT
      vector<double> z(xn);
      ifft<double>(z, y);
      for (int64_t i = 0; i < xn; i++)
        z[i] /= n;
      for (int64_t i = 0; i < xn; i++) {
        const auto e = abs(x[i]-z[i]);
        ASSERT_LE(e, tol) << format("p %d, n %d, xn %d, i %d, e %g:\nx = %g\nz = %g\ny = %g",
                                    p, n, xn, i, e, x, z, sy);
      }
    }
  }
}

TEST(rfft) {
  mt19937 rand(7);
  uniform_real_distribution<double> uniform(-1, 1);
  for (int p = -1; p <= 15; p++) {
    for (const bool full : {true, false}) {
      const double tol = 2.1e-14 * max(p, 0);
      const int64_t n = p < 0 ? 0 : int64_t(1) << p;
      if (n == 1) continue;
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
      vector<Complex<double>> y(n/2);
      rfft<double>(y, x);
      for (int64_t i = 0; i < n; i++) {
        const auto yi = i == 0 ? Complex<double>(y[0].r, 0)
                      : 2*i == n ? Complex<double>(y[0].i, 0)
                      : 2*i < n ? y[i] : conj(y[n-i]);
        const auto e = abs(yi-sy[i]);
        ASSERT_LE(e, tol)
            << format("p %d, n %d, xn %d, i %d, e %g%s", p, n, xn, i, e,
                      n > 8 ? "" : format(":\nx = %g\ny = %g\nsy = %g", x, y, sy));
      }

      // Inverse FFT
      vector<double> z(xn);
      irfft<double>(z, y);
      for (int64_t i = 0; i < xn; i++)
        z[i] /= n;
      for (int64_t i = 0; i < xn; i++) {
        const auto e = abs(x[i]-z[i]);
        ASSERT_LE(e, tol) << format("p %d, n %d, xn %d, i %d, e %g:\nx = %g\nz = %g\ny = %g",
                                    p, n, xn, i, e, x, z, sy);
      }
    }
  }
}

TEST(srfft) {
  mt19937 rand(7);
  uniform_real_distribution<double> uniform(-1, 1);
  for (int p = -1; p <= 15; p++) {
    for (const bool full : {true, false}) {
      const double tol = 2.1e-14 * max(p, 0);
      const int64_t n = p < 0 ? 0 : int64_t(1) << p;
      if (n == 1) continue;
      const int64_t xn = full ? n : uniform_int_distribution<int64_t>(0, n)(rand);
      vector<double> x(xn);
      for (int64_t i = 0; i < xn; i++)
        x[i] = uniform(rand);

      // Accurate forward FFT using arb.
      // The shifted DFT is
      //   y_k = sum_j w_n^(j(k+1/2)) x_j
      //       = sum_j w_n^(jk) w_(2n)^j x_j
      // which we can compute using some premultiplication of x
      const int prec = 200;
      acb_ptr v = _acb_vec_init(n);
      Arb jn, xj;
      for (int64_t j = 0; j < xn; j++) {
        // v[j] = e^(2𝜋i(-j)/(2n)) x[j] = e^(𝜋i(-j/n)) = cos(𝜋(-j/n)) + i sin(𝜋(-j/n))
        const auto vj = v + j;
        arb_set_si(jn, -j);
        arb_div_ui(jn, jn, n, prec);
        arb_sin_cos_pi(acb_imagref(vj), acb_realref(vj), jn, prec);
        arb_set_d(xj, x[j]);
        acb_mul_arb(vj, vj, xj, prec);
      }
      acb_dft(v, v, n, prec);
      vector<Complex<double>> sy(n/2);
      const auto mid = [](arb_t x) { return arf_get_d(arb_midref(x), ARF_RND_NEAR); };
      for (int64_t i = 0; i < n/2; i++) {
        const auto u = v + i;
        sy[i] = Complex<double>(mid(acb_realref(u)), mid(acb_imagref(u)));
      }
      _acb_vec_clear(v, n);

      // Forward FFT
      vector<Complex<double>> y(n/2);
      srfft<double>(y, x);
      for (int64_t i = 0; i < n/2; i++) {
        const auto e = abs(y[i]-sy[i]);
        ASSERT_LE(e, tol)
            << format("srfft: p %d, n %d, xn %d, i %d, e %g%s", p, n, xn, i, e,
                      n > 16 ? "" : format(":\nx = %g\ny = %g\nsy = %g", x, y, sy));
      }

      // Inverse FFT
      vector<double> z(xn);
      isrfft<double>(z, y);
      for (int64_t i = 0; i < xn; i++)
        z[i] /= n/2;
      for (int64_t i = 0; i < xn; i++) {
        const auto e = abs(x[i]-z[i]);
        ASSERT_LE(e, tol) << format("isrfft: p %d, n %d, xn %d, i %d, e %g:\nx = %g\nz = %g\ny = %g",
                                    p, n, xn, i, e, x, z, sy);
      }
    }
  }
}

}  // namespace
}  // namespace mandelbrot
