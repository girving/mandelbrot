// Fast Fourier Transforms

#include "fft.h"
#include "arith.h"
#include "array.h"
#include "bit.h"
#include "cutil.h"
#include "device.h"
#include "debug.h"
#include "expansion.h"
#include "gen-mul-bases.h"
#include "loops.h"
#include "nearest.h"
#include "print.h"
#include "singleton.h"
#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <type_traits>
#include <vector>
namespace mandelbrot {

using std::cos;
using std::is_trivially_copyable_v;
using std::max;
using std::sin;
using std::swap;
using std::unique_ptr;
using std::vector;

namespace {
template<class S> struct FullTwiddleView;

template<class T> class FullTwiddle : public Singleton<FullTwiddle<T>> {
  static_assert(!is_device<T>);
  typedef Undevice<T> S;
  template<class U> friend class FullTwiddle;

  // twiddles = concat(
  //   [sentinel],  // For better cache alignment
  //   twiddle(:1, 2),
  //   twiddle(:2, 4),
  //   twiddle(:4, 8),
  //   ...,
  //   twiddle(:1<<(p-1), 2<<(p-1)))
  vector<Complex<S>> twiddles;
public:
  FullTwiddle(Single s) : Singleton<FullTwiddle<T>>(s) {}

  // Ensure that we know up to twiddle(:1<<s, 2<<s)
  void ensure(const int s) {
    const int64_t size = int64_t(1) << (s+1);
    const int64_t prev = twiddles.size();
    if (size <= prev) return;
    twiddles.reserve(max(size, int64_t(1024)));
    twiddles.resize(size);
    for (int p = 0; p <= s; p++) {
      const auto m = int64_t(1) << p;
      if (prev < 2*m)
        nearest_twiddles(span<Complex<S>>(twiddles.data() + m, m), 2*m, 200);
    }
  }

  void clear() { twiddles.clear(); }
  auto view() const { return FullTwiddleView<S>{twiddles.data(), twiddles.size()}; }

  // twiddle(j, 2<<s)
  Complex<S> operator()(const int64_t j, const int s) const {
    // i = j + sum_{a < s} 2^a = j + 2^s - 1
    const auto m = int64_t(1) << s;
    const auto i = m + j;
    assert(size_t(j) < size_t(m) && size_t(i) < twiddles.size());
    return twiddles[i];
  }
};

// The GPU version computes on the CPU, then moves across
template<class S> class FullTwiddle<Device<S>> : public Singleton<FullTwiddle<Device<S>>> {
  Array<Device<Complex<S>>> gpu;
public:
  FullTwiddle(Single s) : Singleton<FullTwiddle<Device<S>>>(s) {}

  void ensure(const int s) {
    auto& cpu = FullTwiddle<S>::single();
    cpu.ensure(s);
    if (gpu.size() < int64_t(cpu.twiddles.size())) {
      gpu.clear();
      Array<Device<Complex<S>>>(cpu.twiddles.size()).swap(gpu);
      host_to_device<Complex<S>>(gpu, cpu.twiddles);
    }
  }

  void clear() { gpu.clear(); }
  auto view() const { return FullTwiddleView<Device<S>>{gpu.data(), gpu.size()}; }
};

// Lightweight view of already computed twiddles
template<class S> struct FullTwiddleView {
  const AddComplex<S>* twiddles;
  const int size;
  __host__ __device__ Complex<S> operator()(const int j, const int s) const {
    const int m = 1 << s;
    const int i = m + j;
    assert(unsigned(j) < unsigned(m) && unsigned(i) < unsigned(size));
    return twiddles[i];
  }
};
}  // namespace

template<class S> static FullTwiddleView<S> undevice(FullTwiddleView<Device<S>> v) {
  return FullTwiddleView<S>{undevice(v.twiddles), v.size};
}

// Bit-reverse y in place.
// Don't worry about speed since we don't use bit reversal in fft_mul.
template<class T> static void bitrev(span<T> y) {
  const int64_t n = y.size();
  const int p = countr_zero(uint64_t(n));
  const auto np = 64 - p;
  for (int64_t i = 0; i < n; i++) {
    const int64_t j = bitreverse(uint64_t(i)) >> np;
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
  reverses(span<T>(p, n/2));
  reverses(span<T>(p + n/2, n/2));
  std::reverse(p + n/2, p + n);
}
template<class T> static void unreverses(span<T> y) {
  const int64_t n = y.size();
  if (n <= 2) return;
  const auto p = y.data();
  std::reverse(p + n/2, p + n);
  unreverses(span<T>(p, n/2));
  unreverses(span<T>(p + n/2, n/2));
}
template<class T> static void shifted_permute(span<T> y) {
  bitrev(y);
  reverses(y);
}
template<class T> static void shifted_unpermute(span<T> y) {
  unreverses(y);
  bitrev(y);
}

DEF_LOOP(fft_butterfly_0, two_n, i, (S* ys, const S* x, const int xn),
  ys[i] = i < xn ? x[i] : 0;)

DEF_LOOP(fft_butterfly_s, n2, jk, (Complex<S>* y, FullTwiddleView<S> twiddle, const int s),
  const auto m = 1 << s;
  const auto j = jk & (m-1);
  const auto k = (jk - j) << 1;
  auto& y0 = y[k + j];
  auto& y1 = y[k + j + m];
  const auto u0 = y0 + y1;
  const auto u1 = y0 - y1;
  y0 = u0;
  y1 = u1 * conj(twiddle(j, s));)

// fft, but skip the bitrev at the end
template<class S> static void fft_bitrev(span<AddComplex<S>> y, span<const S> x) {
  // Find our power of two
  const int64_t n = y.size(), xn = x.size();
  if (!n) return;
  const int p = countr_zero(uint64_t(n));
  slow_assert(n == int64_t(1) << p && xn <= 2*n);

  // Precompute twiddle factors
  auto& twiddle = FullTwiddle<S>::single();
  twiddle.ensure(p-1);

  // Copy from x to y, without bit reversing
  fft_butterfly_0(2*n, reinterpret_cast<S*>(y.data()), x.data(), xn);

  // Cooley-Tukey FFT
  for (int s = p-1; s >= 0; s--)
    fft_butterfly_s(n/2, y.data(), twiddle.view(), s);
}

DEF_LOOP(ifft_butterfly_0, xn, i, (S* x, const S* ys),
  x[i] = ys[i];)

DEF_LOOP(ifft_butterfly_s, n2, jk, (Complex<S>* y, FullTwiddleView<S> twiddle, const int s),
  const auto m = 1 << s;
  const auto j = jk & (m-1);
  const auto k = (jk - j) << 1;
  auto& y0 = y[k + j];
  auto& y1 = y[k + j + m];
  const auto u0 = y0;
  const auto u1 = y1 * twiddle(j, s);
  y0 = u0 + u1;
  y1 = u0 - u1;)

// ifft, but skip the bitrev at the beginning
template<class S> static void ifft_bitrev(span<S> x, span<AddComplex<S>> y) {
  // Find our power of two
  const int64_t n = y.size(), xn = x.size();
  if (!n) return;
  const int p = countr_zero(uint64_t(n));
  slow_assert(n == int64_t(1) << p && xn <= 2*n);

  // Precompute twiddle factors
  auto& twiddle = FullTwiddle<S>::single();
  twiddle.ensure(p-1);

  // Cooley-Tukey FFT
  for (int s = 0; s < p; s++)
    ifft_butterfly_s(n/2, y.data(), twiddle.view(), s);

  // Copy to output
  ifft_butterfly_0(xn, x.data(), reinterpret_cast<S*>(y.data()));
}

template<class S> void fft(span<Complex<S>> y, span<const S> x) {
  fft_bitrev(y, x);
  bitrev(y);
}

template<class S> void ifft(span<S> x, span<Complex<S>> y) {
  bitrev(y);
  ifft_bitrev(x, y);
}

// bitreverse(-bitreverse(i))
__host__ __device__ static inline int bitrev_neg(const int i) {
  return i ? i ^ ((1 << (31 - countl_zero(uint32_t(i)))) - 1) : 0;
}

DEF_SERIAL(rfft_post_n2, (Complex<S>* y),
  const auto a = y[0];
  y[0] = Complex<S>(a.r + a.i, a.r - a.i);)

DEF_LOOP(rfft_post, n4, j, (Complex<S>* y, FullTwiddleView<S> twiddle, const int p),
  const int j1 = j ? bitrev_neg(2*j) : 1;
  const auto a = y[2*j];
  const auto b = y[j1];
  Complex<S> s;
  Complex<S> t;
  if (j) {
    const int i = bitreverse(uint32_t(j)) >> (34 - p);
    const auto e = left(conj(twiddle(i, p-1)));
    const auto u = a + conj(b);
    const auto v = e * (a - conj(b));
    s = half(u - v);
    t = half(conj(u + v));
  } else {
    s = Complex<S>(a.r + a.i, a.r - a.i);
    t = conj(b);
  }
  y[2*j] = s;
  y[j1] = t;)

DEF_LOOP(irfft_pre, n4, j, (Complex<S>* y, FullTwiddleView<S> twiddle, const int p),
  const int j1 = j ? bitrev_neg(2*j) : 1;
  const auto s = y[2*j];
  const auto t = y[j1];
  Complex<S> a;
  Complex<S> b;
  if (j) {
    const int i = bitreverse(uint32_t(j)) >> (34 - p);
    const auto e = right(twiddle(i, p-1));
    const auto u = conj(t) + s;
    const auto v = e * (conj(t) - s);
    a = u + v;
    b = conj(u - v);
  } else {
    a.r = s.r + s.i;
    a.i = s.r - s.i;
    b = twice(conj(t));
  }
  y[2*j] = a;
  y[j1] = b;)

template<class S> static void rfft_bitrev(span<AddComplex<S>> y, span<const S> x) {
  const int64_t n = 2*y.size();
  if (!n) return;

  // Precompute twiddle factors
  const int p = countr_zero(uint64_t(n));
  auto& twiddle = FullTwiddle<S>::single();
  twiddle.ensure(p-1);

  // Half-size complex FFT
  fft_bitrev(y, x);

  // Postprocess into real FFT
  if (p == 1)
    rfft_post_n2(y.data());
  else
    rfft_post(n/4, y.data(), twiddle.view(), p);
}

template<class S> static void irfft_bitrev(span<S> x, span<AddComplex<S>> y) {
  const int64_t n = 2*y.size();
  if (!n) return;

  // Precompute twiddle factors
  const int p = countr_zero(uint64_t(n));
  auto& twiddle = FullTwiddle<S>::single();
  twiddle.ensure(p-1);

  // Preprocess into half-size complex FFT
  if (p == 1)
    rfft_post_n2(y.data());
  else
    irfft_pre(n/4, y.data(), twiddle.view(), p);

  // Half-size complex inverse FFT
  ifft_bitrev(x, y);
}

template<class S> void rfft(span<Complex<S>> y, span<const S> x) {
  rfft_bitrev(y, x);
  bitrev(y);
}

template<class S> void irfft(span<S> x, span<Complex<S>> y) {
  bitrev(y);
  irfft_bitrev(x, y);
}

// First butterfly, copying real to complex, without twiddle factors:
//   t j2 j1 j0 -> t k2/2 j1 j0
DEF_LOOP(srfft_butterfly_0, n2, j, (Complex<S>* y, const S* x, const int64_t xn),
  const auto x0 = j < xn ? x[j] : 0;
  const auto x1 = j + n2 < xn ? x[j + n2] : 0;
  y[j] = Complex<S>(x0, -x1);)

// Second butterfly, in place, shifted twiddling on input:
//   t k(p-1)/2 j(p-2) ...j... -> k(p-1) t k(p-2)/2 ...j...
DEF_LOOP(srfft_butterfly_1, n4, j, (Complex<S>* y, FullTwiddleView<S> twiddle, const S sqrt_half, const int p),
  auto& y0 = y[j];
  auto& y1 = y[j + n4];
  const auto z1 = diag<-1>(sqrt_half, y1);
  const auto u0 = conj(twiddle(j, p)) * (y0 + z1);
  const auto u1 = conj(twiddle(3*j, p) * (y0 - z1));
  y0 = u0;
  y1 = u1;)

// Remaining butterflies, in place, unshifted twiddling on input:
//   ...k... t k(s+1)/2 js ...j... -> ...k... k(s+1) t ks/2 ...j...
DEF_LOOP(srfft_butterfly_s, n4, jk, (Complex<S>* y, FullTwiddleView<S> twiddle, const int s),
  const auto m = int64_t(1) << s;
  const auto j = jk & (m-1);
  const auto k = (jk - j) << 1;
  auto& y0 = y[k + j];
  auto& y1 = y[k + j + m];
  const auto u0 = y0 + y1;
  const auto u1 = conj(twiddle(j, s) * (y0 - y1));
  y0 = u0;
  y1 = u1;)

// srfft with permuted output
template<class T> static void srfft_scramble(span<AddComplex<T>> y, span<const T> x) {
  typedef Undevice<T> S;

  // Find our power of two
  const int64_t n = 2*y.size(), xn = x.size();
  if (!n) return;
  const int p = countr_zero(uint64_t(n));
  slow_assert(n == int64_t(1) << p && xn <= n);

  // Precompute twiddle factors
  auto& twiddle = FullTwiddle<T>::single();
  twiddle.ensure(p);

  // Decimination-in-frequency shifted real-to-complex FFT, using the commutator notation:
  //   t j2 j1 j0 -> t k2/2 j1 j0
  //              -> k2 t k1/2 j0
  //              -> k2 k1 t k0/2
  //              == k2 k1 k0

  // First butterfly, copying real to complex, without twiddle factors:
  //   t j2 j1 j0 -> t k2/2 j1 j0
  srfft_butterfly_0(n/2, y.data(), x.data(), xn);

  // Second butterfly, in place, shifted twiddling on input:
  //   t k(p-1)/2 j(p-2) ...j... -> k(p-1) t k(p-2)/2 ...j...
  if (p > 1) {
    static const auto sqrt_half = nearest_sqrt<S>(1, 2);
    srfft_butterfly_1(n/4, y.data(), twiddle.view(), sqrt_half, p);
  }

  // Remaining butterflies, in place, unshifted twiddling on input:
  //   ...k... t k(s+1)/2 js ...j... -> ...k... k(s+1) t ks/2 ...j...
  for (int s = p-3; s >= 0; s--)
    srfft_butterfly_s(n/4, y.data(), twiddle.view(), s);
}

// Final butterfly, copying complex to real, without twiddle factors:
//   t j2 j1 j0 <- t k2/2 j1 j0
DEF_LOOP(isrfft_butterfly_0, n2, j, (S* x, const Complex<S>* y, const int xn),
  const auto yj = y[j];
  if (j < xn) x[j] = yj.r;
  if (j + n2 < xn) x[j + n2] = -yj.i;)

// Second to last butterfly, in place, shifted twiddling on output:
//   t k(p-1)/2 j(p-2) ...j... -> k(p-1) t k(p-2)/2 ...j...
DEF_LOOP(isrfft_butterfly_1, n4, j, (Complex<S>* y, FullTwiddleView<S> twiddle, const S sqrt_half, const int p),
  auto& y0 = y[j];
  auto& y1 = y[j + n4];
  const auto z0 = twiddle(j, p) * y0;
  const auto z1 = twiddle(3*j, p) * y1;
  const auto u0 = z0 + conj(z1);
  const auto u1 = diag<1>(sqrt_half, z0 - conj(z1));
  y0 = u0;
  y1 = u1;)

// First butterflies, in place, unshifted twiddling on output:
//   ...k... t k(s+1)/2 js ...j... <- ...k... k(s+1) t ks/2 ...j...
DEF_LOOP(isrfft_butterfly_s, n4, jk, (Complex<S>* y, FullTwiddleView<S> twiddle, const int s),
  const auto m = int64_t(1) << s;
  const auto j = jk & (m-1);
  const auto k = (jk - j) << 1;
  auto& y0 = y[k + j];
  auto& y1 = y[k + j + m];
  const auto u0 = y0;
  const auto u1 = conj(twiddle(j, s) * y1);
  y0 = u0 + u1;
  y1 = u0 - u1;)

// isrfft, but taking permuted input
template<class T> static void isrfft_scramble(span<T> x, span<AddComplex<T>> y) {
  typedef Undevice<T> S;

  // Find our power of two
  const int64_t n = 2*y.size(), xn = x.size();
  if (!n) return;
  const int p = countr_zero(uint64_t(n));
  slow_assert(n == int64_t(1) << p && xn <= 2*n);

  // Precompute twiddle factors
  auto& twiddle = FullTwiddle<T>::single();
  twiddle.ensure(p);

  // First butterflies, in place, unshifted twiddling on output:
  //   ...k... t k(s+1)/2 js ...j... <- ...k... k(s+1) t ks/2 ...j...
  for (int s = 0; s < p-2; s++)
    isrfft_butterfly_s(n/4, y.data(), twiddle.view(), s);

  // Second to last butterfly, in place, shifted twiddling on output:
  //   t k(p-1)/2 j(p-2) ...j... -> k(p-1) t k(p-2)/2 ...j...
  if (p > 1) {
    static const auto sqrt_half = nearest_sqrt<S>(1, 2);
    isrfft_butterfly_1(n/4, y.data(), twiddle.view(), sqrt_half, p);
  }

  // Final butterfly, copying complex to real, without twiddle factors:
  //   t j2 j1 j0 <- t k2/2 j1 j0
  isrfft_butterfly_0(n/2, x.data(), y.data(), xn);
}

template<class S> void srfft(span<Complex<S>> y, span<const S> x) {
  srfft_scramble(y, x);
  shifted_permute(y);
}

template<class S> void isrfft(span<S> x, span<Complex<S>> y) {
  shifted_unpermute(y);
  isrfft_scramble(x, y);
}

DEF_LOOP(mul_cwise_loop, fn2, i, (Complex<S>* fx, const Complex<S>* fy, const S a),
  fx[i] = a * fx[i] * fy[i];)

DEF_LOOP(sqr_cwise_loop, fn2, i, (Complex<S>* fx, const S a),
  fx[i] = a * sqr(fx[i]);)

DEF_LOOP(rmul_cwise_loop, fn2, i, (Complex<S>* fx, const Complex<S>* fy, const S a),
  const auto x = fx[i];
  const auto y = fy[i];
  fx[i] = a * (i ? x * y : hadamard(x, y));)

DEF_LOOP(rsqr_cwise_loop, fn2, i, (Complex<S>* fx, const S a),
  const auto x = fx[i];
  fx[i] = a * (i ? sqr(x) : hadamard_sqr(x));)

// Whether to use rfft or srfft for multiplication
static const bool srfft_mul = true;

template<class T> void fft_mul(span<T> z, span<add_const_t<T>> x, span<add_const_t<T>> y) {
  typedef Undevice<T> S;
  const int64_t nz = z.size(), nx = x.size(), ny = y.size();
  slow_assert(nz <= relu(nx + ny - 1));
  if (nz <= mul_base_n)
    mul_base(z.data(), nz, x.data(), nx, y.data(), ny);
  else {
    // FFT multiplication for large n
    const int64_t fn = bit_ceil(uint64_t(2*nz));
    const Array<AddComplex<T>> buffer(fn);
    const auto fx = buffer.span().first(fn/2);
    const auto fy = buffer.span().last(fn/2);
    if (srfft_mul) {
      srfft_scramble(fx, x);
      srfft_scramble(fy, y);
      mul_cwise_loop(fn/2, fx.data(), fy.data(), inv(S(fn/2)));
      isrfft_scramble(z, fx);
    } else {
      rfft_bitrev(fx, x);
      rfft_bitrev(fy, y);
      rmul_cwise_loop(fn/2, fx.data(), fy.data(), inv(S(fn)));
      irfft_bitrev(z, fx);
    }
  }
}

template<class T> void fft_sqr(span<T> y, span<add_const_t<T>> x) {
  typedef Undevice<T> S;
  const int64_t ny = y.size(), nx = x.size();
  slow_assert(ny <= relu(2*nx - 1));
  if (ny <= sqr_base_n)
    sqr_base(y.data(), ny, x.data(), nx);
  else {
    // FFT squaring for large n
    const int64_t fn = bit_ceil(uint64_t(2*ny));
    const Array<AddComplex<T>> fx(fn/2);
    if (srfft_mul) {
      srfft_scramble(fx, x);
      sqr_cwise_loop(fn/2, fx.data(), inv(S(fn/2)));
      isrfft_scramble(y, fx);
    } else {
      rfft_bitrev(fx, x);
      rsqr_cwise_loop(fn/2, fx.data(), inv(S(fn)));
      irfft_bitrev(y, fx);
    }
  }
}

#define MUL(T) \
  template void fft_mul(span<T> z, span<const T> x, span<const T> y); \
  template void fft_sqr(span<T> y, span<const T> x);
#define REST(S) \
  template void fft(span<Complex<S>> y, span<const S> x); \
  template void ifft(span<S> x, span<Complex<S>> y); \
  template void rfft(span<Complex<S>> y, span<const S> x); \
  template void irfft(span<S> x, span<Complex<S>> y); \
  template void srfft(span<Complex<S>> y, span<const S> x); \
  template void isrfft(span<S> x, span<Complex<S>> y);
REST(double)
MUL(double)
MUL(Expansion<2>)
IF_CUDA(
  MUL(Device<double>)
  MUL(Device<Expansion<2>>)
)

}  // namespace mandelbrot
