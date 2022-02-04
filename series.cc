// Series arithmetic

#include "series.h"
#include "expansion.h"
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

DEF_SERIAL(add_scalar_kernel, (S* ys, const S a), ys[0] += a;)

template<class T> void add_scalar(Series<T>& x, const typename Series<T>::Scalar a) {
  slow_assert(x.nonzero());
  add_scalar_kernel(x.data(), a);
}

DEF_LOOP(high_addsub_loop, n, i, (S* y, const S* x, const int ynz, const int xnz, const int sign, const int s),
  const auto yi = i < ynz ? y[i] : S(0);
  auto xi = uint32_t(i-s) < uint32_t(xnz) ? x[i-s] : S(0);
  if (sign < 0) xi = -xi;
  y[i] = yi + xi;)

template<class T> void high_addsub(Series<T>& y, const int sign, const int64_t s, SeriesView<add_const_t<T>> x) {
  const auto ynz = y.nonzero(), xnz = x.nonzero();
  const auto nk = min(y.known(), x.known() + s);
  const auto nz = min(nk, max(ynz, xnz ? xnz + s : 0));
  slow_assert(abs(sign) == 1 && !y.alias(x) && nz <= y.limit());
  high_addsub_loop(nz, y.data(), x.data(), ynz, xnz, sign, s);
  y.set_counts(nk, nz);
}

DEF_LOOP(mul1p_post_loop, post, i, (S* z, const S* x, const int s, const int xnz),
  z[i] = (i < s ? S(0) : z[i]) + (i < xnz ? x[i] : S(0));)

template<class T> void mul1p_post(Series<T>& z, SeriesView<add_const_t<T>> x,
                                  const int64_t post, const int64_t s, const int64_t xnz) {
  mul1p_post_loop(post, z.data(), x.data(), s, xnz);
}

#define Ss(S) \
  template void add_scalar(Series<S>&, const typename Series<S>::Scalar); \
  template void high_addsub(Series<S>&, const int, const int64_t, SeriesView<const S>); \
  template void mul1p_post(Series<S>&, SeriesView<const S>, const int64_t, const int64_t, const int64_t);
Ss(double)
Ss(Expansion<2>)
IF_CUDA(
  Ss(Device<double>)
  Ss(Device<Expansion<2>>)
)

}  // namespace mandelbrot
