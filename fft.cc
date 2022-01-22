// Fast Fourier Transforms

#include "fft.h"
#include "debug.h"
#include <bit>
#include <cmath>
#include <cstdint>
#include <memory>
namespace mandelbrot {

using std::bit_ceil;
using std::cos;
using std::countr_zero;
using std::has_single_bit;
using std::sin;
using std::unique_ptr;

static Complex<double> twiddle(const int64_t a, const int64_t b) {
  const double t = 2 * M_PI / b * a;
  return Complex<double>(cos(t), sin(t));
}

template<class S> void fft(Complex<S>* y, const S* x, const int64_t n, const int64_t xn) {
  // Find our power of two
  if (!n) return;
  const int p = countr_zero(uint64_t(n));
  slow_assert(n == int64_t(1) << p && xn <= n);

  // Bit-reverse from x to y, and expand from real to complex
  const auto np = 64 - p;
  for (int64_t i = 0; i < n; i++)
    y[__builtin_bitreverse64(i) >> np] = Complex<S>(i < xn ? x[i] : 0);

  // Cooley-Tukey FFT
  //   Details: https://en.wikipedia.org/wiki/Cooley-Tukey_FFT_algorithm
  for (int s = 0; s < p; s++) {
    const auto m = int64_t(1) << s;
    for (int64_t k = 0; k < n; k += 2*m) {
      for (int64_t j = 0; j < m; j++) {
        auto& y0 = y[k + j];
        auto& y1 = y[k + j + m];
        const auto u0 = y0;
        const auto u1 = y1 * twiddle(-j, 2*m);
        y0 = u0 + u1;
        y1 = u0 - u1;
      }
    }
  }
}

template<class S> void ifft(S* x, Complex<S>* y, const int64_t n, const int64_t xn) {
  // Find our power of two
  if (!n) return;
  const int p = countr_zero(uint64_t(n));
  slow_assert(n == int64_t(1) << p && xn <= n);

  // Cooley-Tuken inverse FFT
  for (int s = p-1; s >= 0; s--) {
    const auto m = int64_t(1) << s;
    for (int64_t k = 0; k < n; k += 2*m) {
      for (int64_t j = 0; j < m; j++) {
        auto& y0 = y[k + j];
        auto& y1 = y[k + j + m];
        const auto u0 = y0 + y1;
        const auto u1 = y0 - y1;
        y0 = u0;
        y1 = u1 * twiddle(j, 2*m);
      }
    }
  }

  // Bit-reverse from y to x, and take real parts
  const auto np = 64 - p;
  for (int64_t i = 0; i < xn; i++)
    x[i] = y[__builtin_bitreverse64(i) >> np].r;
}

// z[:n] = x[:n] * y[:n]
template<class S> void fft_mul(S* z, const S* x, const S* y, const int64_t n) {
  if (n == 0)
    return;
  else if (n == 1)
    z[0] = x[0] * y[0];
  else {
    // FFT multiplication for large n
    const int64_t fn = bit_ceil(uint64_t(2*n));
    unique_ptr<Complex<S>[]> buffer(new Complex<S>[2*fn]);
    const auto fx = buffer.get();
    const auto fy = buffer.get() + fn;
    fft(fx, x, fn, n);
    fft(fy, y, fn, n);
    const auto a = inv(S(fn));
    for (int64_t i = 0; i < fn; i++)
      fx[i] = a * fx[i] * fy[i];
    ifft(z, fx, fn, n);
  }
}

// y[:n] = x[:n]^2
template<class S> void fft_sqr(S* y, const S* x, const int64_t n) {
  if (n == 0)
    return;
  else if (n == 1)
    y[0] = sqr(x[0]);
  else {
    // FFT squaring for large n
    const int64_t fn = bit_ceil(uint64_t(2*n));
    unique_ptr<Complex<S>[]> buffer(new Complex<S>[fn]);
    const auto fx = buffer.get();
    fft(fx, x, fn, n);
    const auto a = inv(S(fn));
    for (int64_t i = 0; i < fn; i++)
      fx[i] = a * sqr(fx[i]);
    ifft(y, fx, fn, n);
  }
}

#define FFT(S) \
  template void fft(Complex<S>* y, const S* x, const int64_t n, const int64_t xn); \
  template void ifft(S* x, Complex<S>* y, const int64_t n, const int64_t xn); \
  template void fft_mul(S* z, const S* x, const S* y, const int64_t n); \
  template void fft_sqr(S* y, const S* x, const int64_t n);
FFT(double)

}  // namespace mandelbrot
