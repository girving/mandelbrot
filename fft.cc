// Fast Fourier Transforms

#include "fft.h"
#include "debug.h"
#include "print.h"
#include <bit>
#include <cmath>
#include <cstdint>
#include <memory>
#include <type_traits>
#include <vector>
namespace mandelbrot {

using std::bit_ceil;
using std::cos;
using std::countr_zero;
using std::has_single_bit;
using std::is_trivially_copyable_v;
using std::max;
using std::sin;
using std::swap;
using std::unique_ptr;
using std::vector;

static Complex<double> slow_twiddle(const int64_t a, const int64_t b) {
  const double t = 2 * M_PI / b * a;
  return Complex<double>(cos(t), sin(t));
}

namespace {
template<class S> class FullTwiddle {
  // twiddles = concat(
  //   [sentinel],  // For better cache alignment
  //   twiddle(:1, 2),
  //   twiddle(:2, 4),
  //   twiddle(:4, 8),
  //   ...,
  //   twiddle(:1<<(p-1), 2<<(p-1)))
  static vector<Complex<S>> twiddles;
  static int p;
public:

  // Ensure that we know up to twiddle(:1<<s, 2<<s)
  void ensure(const int s) {
    slow_assert(0 <= s && s <= 22);
    if (s < p) return;
    twiddles.reserve(size_t(1) << max(10, s+1));
    twiddles.resize(size_t(1) << (s+1));
    while (p <= s) {
      const auto m = int64_t(1) << p;
      for (int64_t j = 0; j < m; j++)
        twiddles[m + j] = slow_twiddle(j, 2*m);
      p++;
    }
  }

  // slow_twiddle(j, 2<<s)
  Complex<S> operator()(const int64_t j, const int s) {
    // i = j + sum_{a < s} 2^a = j + 2^s - 1
    const auto m = int64_t(1) << s;
    const auto i = m + j;
    assert(size_t(j) < size_t(m) && size_t(i) < twiddles.size());
    return twiddles[i];
  }
};

template<class S> vector<Complex<S>> FullTwiddle<S>::twiddles;
template<class S> int FullTwiddle<S>::p;
}  // namespace

// Bit-reverse y in place.
// Don't worry about speed since we don't use bit reversal in fft_mul.
template<class T> static void bitrev(span<T> y) {
  const int64_t n = y.size();
  const int p = countr_zero(uint64_t(n));
  const auto np = 64 - p;
  for (int64_t i = 0; i < n; i++) {
    const int64_t j = __builtin_bitreverse64(i) >> np;
    if (i < j)
      swap(y[i], y[j]);
  }
}

// Shifted-permute y in place, very slowly.
// Ideally I'd understand why this is the transformation.
template<class T> static void reverses(span<T> y) {
  const int64_t n = y.size();
  if (n <= 2) return;
  const auto p = y.data();
  reverses(span<T>(p, p + n/2));
  reverses(span<T>(p + n/2, p + n));
  std::reverse(p + n/2, p + n);
}
template<class T> static void unreverses(span<T> y) {
  const int64_t n = y.size();
  if (n <= 2) return;
  const auto p = y.data();
  std::reverse(p + n/2, p + n);
  unreverses(span<T>(p, p + n/2));
  unreverses(span<T>(p + n/2, p + n));
}
template<class T> static void shifted_permute(span<T> y) {
  bitrev(y);
  reverses(y);
}
template<class T> static void shifted_unpermute(span<T> y) {
  unreverses(y);
  bitrev(y);
}

// fft, but skip the bitrev at the end
template<class S> static void fft_bitrev(span<Complex<S>> y, span<const S> x) {
  // Find our power of two
  const int64_t n = y.size(), xn = x.size();
  if (!n) return;
  const int p = countr_zero(uint64_t(n));
  slow_assert(n == int64_t(1) << p && xn <= 2*n);

  // Copy from x to y, without bit reversing
  static_assert(is_trivially_copyable_v<S>);
  memcpy(y.data(), x.data(), xn*sizeof(S));
  memset(reinterpret_cast<S*>(y.data()) + xn, 0, (2*n-xn)*sizeof(S));

  // Cooley-Tukey FFT
  for (int s = p-1; s >= 0; s--) {
    const auto m = int64_t(1) << s;
    for (int64_t k = 0; k < n; k += 2*m) {
      for (int64_t j = 0; j < m; j++) {
        auto& y0 = y[k + j];
        auto& y1 = y[k + j + m];
        const auto u0 = y0 + y1;
        const auto u1 = y0 - y1;
        y0 = u0;
        y1 = u1 * slow_twiddle(-j, 2*m);
      }
    }
  }
}

// ifft, but skip the bitrev at the beginning
template<class S> static void ifft_bitrev(span<S> x, span<Complex<S>> y) {
  // Find our power of two
  const int64_t n = y.size(), xn = x.size();
  if (!n) return;
  const int p = countr_zero(uint64_t(n));
  slow_assert(n == int64_t(1) << p && xn <= 2*n);

  // Cooley-Tukey FFT
  for (int s = 0; s < p; s++) {
    const auto m = int64_t(1) << s;
    for (int64_t k = 0; k < n; k += 2*m) {
      for (int64_t j = 0; j < m; j++) {
        auto& y0 = y[k + j];
        auto& y1 = y[k + j + m];
        const auto u0 = y0;
        const auto u1 = y1 * slow_twiddle(j, 2*m);
        y0 = u0 + u1;
        y1 = u0 - u1;
      }
    }
  }

  // Copy to output
  static_assert(is_trivially_copyable_v<S>);
  memcpy(x.data(), y.data(), xn*sizeof(S));
}

template<class S> void fft(span<Complex<S>> y, span<const S> x) {
  fft_bitrev(y, x);
  bitrev(y);
}

template<class S> void ifft(span<S> x, span<Complex<S>> y) {
  bitrev(y);
  ifft_bitrev(x, y);
}

template<class S> void rfft(span<Complex<S>> y, span<const S> x) {
  const int64_t n = 2*y.size();
  if (!n) return;

  // Half-size complex FFT
  fft(y, x);

  // Postprocess into real FFT
  const auto c = y[0];
  y[0].r = c.r + c.i;
  y[0].i = c.r - c.i;
  for (int64_t i = 1; i <= n/4; i++) {
    const auto e = left(slow_twiddle(-i, n));
    const auto a = y[i];
    const auto b = y[n/2-i];
    const auto u = a + conj(b);
    const auto v = e * (a - conj(b));
    const auto s = 0.5 * (u - v);
    const auto t = 0.5 * conj(u + v);
    y[i] = s;
    y[n/2-i] = t;
  }
}

template<class S> void irfft(span<S> x, span<Complex<S>> y) {
  const int64_t n = 2*y.size();
  if (!n) return;

  // Preprocess into half-size complex FFT
  const auto c = y[0];
  y[0].r = c.r + c.i;
  y[0].i = c.r - c.i;
  y[0] = y[0];
  for (int64_t i = 1; i <= n/4; i++) {
    const auto e = right(slow_twiddle(i, n));
    const auto m = n >= 8 || 4*i == n ? 1 : 0.5;
    const auto s = y[i];
    const auto t = y[n/2-i];
    const auto u = m * (conj(t) + s);
    const auto v = m * (conj(t) - s);
    const auto a = u + e * v;
    const auto b = conj(u - e * v);
    y[i] = a;
    y[n/2-i] = b;
  }

  // Half-size complex inverse FFT
  ifft(x, y);
}

// srfft with permuted output
template<class S> static void srfft_scramble(span<Complex<S>> y, span<const S> x) {
  // Find our power of two
  const int64_t n = 2*y.size(), xn = x.size();
  if (!n) return;
  const int p = countr_zero(uint64_t(n));
  slow_assert(n == int64_t(1) << p && xn <= n);

  // Precompute twiddle factors
  FullTwiddle<S> T;
  T.ensure(p);

  // Decimination-in-frequency shifted real-to-complex FFT, using the commutator notation:
  //   t j2 j1 j0 -> t k2/2 j1 j0
  //              -> k2 t k1/2 j0
  //              -> k2 k1 t k0/2
  //              == k2 k1 k0

  // First butterfly, copying real to complex, without twiddle factors:
  //   t j2 j1 j0 -> t k2/2 j1 j0
  for (int64_t j = 0; j < n/2; j++) {
    const auto x0 = j < xn ? x[j] : 0;
    const auto x1 = j + n/2 < xn ? x[j + n/2] : 0;
    y[j] = Complex<S>(x0, -x1);
  }

  // Second butterfly, in place, shifted twiddling on input:
  //   t k(p-1)/2 j(p-2) ...j... -> k(p-1) t k(p-2)/2 ...j...
  if (p > 1) {
    for (int64_t j = 0; j < n/4; j++) {
      auto& y0 = y[j];
      auto& y1 = y[j + n/4];
      const auto z1 = diag<-1>(sqrt(S(0.5)), y1);
      const auto u0 = conj(T(j, p)) * (y0 + z1);
      const auto u1 = conj(T(3*j, p) * (y0 - z1));
      y0 = u0;
      y1 = u1;
    }
  }

  // Remaining butterflies, in place, unshifted twiddling on input:
  //   ...k... t k(s+1)/2 js ...j... -> ...k... k(s+1) t ks/2 ...j...
  for (int s = p-3; s >= 0; s--) {
    const auto m = int64_t(1) << s;
    for (int64_t k = 0; k < n/2; k += 2*m) {
      for (int64_t j = 0; j < m; j++) {
        auto& y0 = y[k + j];
        auto& y1 = y[k + j + m];
        const auto u0 = y0 + y1;
        const auto u1 = conj(T(j, s) * (y0 - y1));
        y0 = u0;
        y1 = u1;
      }
    }
  }
}

// isrfft, but taking permuted input
template<class S> static void isrfft_scramble(span<S> x, span<Complex<S>> y) {
  // Find our power of two
  const int64_t n = 2*y.size(), xn = x.size();
  if (!n) return;
  const int p = countr_zero(uint64_t(n));
  slow_assert(n == int64_t(1) << p && xn <= 2*n);

  // Precompute twiddle factors
  FullTwiddle<S> T;
  T.ensure(p);

  // First butterflies, in place, unshifted twiddling on output:
  //   ...k... t k(s+1)/2 js ...j... <- ...k... k(s+1) t ks/2 ...j...
  for (int s = 0; s < p-2; s++) {
    const auto m = int64_t(1) << s;
    for (int64_t k = 0; k < n/2; k += 2*m) {
      for (int64_t j = 0; j < m; j++) {
        auto& y0 = y[k + j];
        auto& y1 = y[k + j + m];
        const auto u0 = y0;
        const auto u1 = conj(T(j, s) * y1);
        y0 = u0 + u1;
        y1 = u0 - u1;
      }
    }
  }

  // Second to last butterfly, in place, shifted twiddling on output:
  //   t k(p-1)/2 j(p-2) ...j... -> k(p-1) t k(p-2)/2 ...j...
  if (p > 1) {
    for (int64_t j = 0; j < n/4; j++) {
      auto& y0 = y[j];
      auto& y1 = y[j + n/4];
      const auto z0 = T(j, p) * y0;
      const auto z1 = T(3*j, p) * y1;
      const auto u0 = z0 + conj(z1);
      const auto u1 = diag<1>(sqrt(S(0.5)), z0 - conj(z1));
      y0 = u0;
      y1 = u1;
    }
  }

  // Final butterfly, copying complex to real, without twiddle factors:
  //   t j2 j1 j0 <- t k2/2 j1 j0
  for (int64_t j = 0; j < n/2; j++) {
    const auto yj = y[j];
    if (j < xn) x[j] = yj.r;
    if (j + n/2 < xn) x[j + n/2] = -yj.i;
  }
}

template<class S> void srfft(span<Complex<S>> y, span<const S> x) {
  srfft_scramble(y, x);
  shifted_permute(y);
}

template<class S> void isrfft(span<S> x, span<Complex<S>> y) {
  shifted_unpermute(y);
  isrfft_scramble(x, y);
}

template<class S> void fft_mul(span<S> z, span<const S> x, span<const S> y) {
  const int64_t n = z.size();
  slow_assert(n == int64_t(x.size()) && n == int64_t(y.size()));
  if (!n)
    return;
  else if (n == 1)
    z[0] = x[0] * y[0];
  else {
    // FFT multiplication for large n
    const int64_t fn = bit_ceil(uint64_t(2*n));
    unique_ptr<Complex<S>[]> buffer(new Complex<S>[fn]);
    const span<Complex<S>> fx(buffer.get(), fn/2);
    const span<Complex<S>> fy(buffer.get() + fn/2, fn/2);
    srfft_scramble(fx, x);
    srfft_scramble(fy, y);
    const auto a = inv(S(fn/2));
    for (int64_t i = 0; i < fn/2; i++)
      fx[i] = a * fx[i] * fy[i];
    isrfft_scramble(z, fx);
  }
}

template<class S> void fft_sqr(span<S> y, span<const S> x) {
  const int64_t n = y.size();
  slow_assert(n == int64_t(x.size()));
  if (!n)
    return;
  else if (n == 1)
    y[0] = sqr(x[0]);
  else {
    // FFT squaring for large n
    const int64_t fn = bit_ceil(uint64_t(2*n));
    unique_ptr<Complex<S>[]> buffer(new Complex<S>[fn/2]);
    const span<Complex<S>> fx(buffer.get(), fn/2);
    srfft_scramble(fx, x);
    const auto a = inv(S(fn/2));
    for (int64_t i = 0; i < fn/2; i++)
      fx[i] = a * sqr(fx[i]);
    isrfft_scramble(y, fx);
  }
}

#define FFT(S) \
  template void fft(span<Complex<S>> y, span<const S> x); \
  template void ifft(span<S> x, span<Complex<S>> y); \
  template void rfft(span<Complex<S>> y, span<const S> x); \
  template void irfft(span<S> x, span<Complex<S>> y); \
  template void srfft(span<Complex<S>> y, span<const S> x); \
  template void isrfft(span<S> x, span<Complex<S>> y); \
  template void fft_mul(span<S> z, span<const S> x, span<const S> y); \
  template void fft_sqr(span<S> y, span<const S> x);
FFT(double)

}  // namespace mandelbrot
