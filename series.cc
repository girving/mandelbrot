// Series arithmetic

#include "series.h"
#include "poly.h"
namespace mandelbrot {

Series<double> approx(const Poly& x, const int64_t n) {
  slow_assert(n >= x.length());
  Series<double> y(n);
  y.set_counts(n, n);
  for (int64_t i = 0; i < n; i++)
    y[i] = arf_get_d(arb_midref(x[i]), ARF_RND_NEAR);
  return y;
}

double error(SeriesView<const double> x, SeriesView<const double> y, const bool relative) {
  if (x.known() != y.known())
    return numeric_limits<double>::infinity();
  double e = 0;
  const auto both = min(x.nonzero(), y.nonzero());
  for (int64_t i = 0; i < both; i++) {
    double d = abs(x[i] - y[i]);
    if (relative)
      d /= max(1., abs(y[i]));
    e = max(e, d);
  }
  for (int64_t i = both; i < x.nonzero(); i++)
    e = max(e, abs(x[i]));
  for (int64_t i = both; i < y.nonzero(); i++) {
    double d = abs(y[i]);
    if (relative)
      d = min(d, 1.);
    e = max(e, d);
  }
  return e;
}

double error(SeriesView<const double> x, initializer_list<double>&& ys, const bool relative) {
  return error(x, Series<double>(move(ys)), relative);
}

double error(SeriesView<const double> x, const Poly& y, const bool relative) {
  if (x.known() < y.length())
    return numeric_limits<double>::infinity();
  return error(x, approx(y, x.known()), relative);
}

template<class S> __global__ static void add_scalar_kernel(S* ys, const S a) { ys[0] += a; }

template<class T, bool view> void Series<T,view>::operator+=(const S a) const {
  slow_assert(nonzero_);
  if constexpr (is_device<T>) add_scalar_kernel<<<1,1,0,stream()>>>(device_get(data()), a);
  else x[0] += a;
}

template<class S> __global__
void high_addsub_kernel(const int n, S* y, const S* x, const int ynz, const int xnz, const int sign, const int s) {
  GRID_STRIDE_LOOP(n, i) {
    const auto yi = i < ynz ? y[i] : S(0);
    auto xi = uint32_t(i-s) < uint32_t(xnz) ? x[i-s] : S(0);
    if (sign < 0) xi = -xi;
    y[i] = yi + xi;
  }
}

template<class T, bool view> void Series<T,view>::high_addsub(const int sign, const int64_t s, SeriesView<CT> f) {
  static_assert(!is_const_v<T>);
  const auto fnz = f.nonzero_;
  const auto nk = min(known_, f.known_ + s);
  const auto nz = min(nk, max(nonzero_, fnz ? fnz + s : 0));
  slow_assert(abs(sign) == 1 && !alias(f) && nz <= limit());
  if (nz) {
    if constexpr (is_device<T>)
      INVOKE_GRID_STRIDE_LOOP(high_addsub_kernel, nz, device_get(data()), device_get(f.data()), nonzero_, fnz, sign, s);
    else {
      #define LOOP(op) \
        for (int64_t i = 0; i < nz; i++) \
          x[i] = (i < nonzero_ ? x[i] : S(0)) op (uint64_t(i-s) < uint64_t(fnz) ? f.x[i-s] : S(0));
      if (sign > 1) LOOP(+)
      else LOOP(-)
      #undef LOOP
    }
  }
  known_ = nk;
  nonzero_ = nz;
}


}  // namespace mandelbrot
