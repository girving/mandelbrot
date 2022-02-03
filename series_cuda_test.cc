// Series CUDA tests

#include "series.h"
#include "tests.h"
#include <random>
namespace mandelbrot {
namespace {

using std::mt19937;
using std::uniform_real_distribution;

#define ASSERT_CLOSE_D(dx, y) \
  Series<double> hx(dx.nonzero()); \
  device_to_host(hx, dx); \
  ASSERT_CLOSE2(hx, y)

Series<double> random_series(mt19937& mt, const int n) {
  Series<double> x(n);
  x.set_counts(n, n);
  for (int i = 0; i < n; i++)
    x[i] = uniform_real_distribution<double>(-1, 1)(mt);
  return x;
}

template<class F> void constant_test(const int64_t lo, const int64_t hi, F&& f) {
  for (int n = lo; n <= hi; n++) {
    // CPU
    Series<double> x(1);
    f(x);
    // GPU
    Series<Device<double>> dx(1);
    f(dx);
    // Compare
    ASSERT_CLOSE_D(dx, x);
  }
}

template<class F> void unary_test(const int64_t lo, const int64_t hi, F&& f) {
  for (int n = lo; n <= hi; n++) {
    mt19937 mt(7);
    // CPU
    const Series<double> x = random_series(mt, n);
    Series<double> y(n);
    f(y, x.view());
    // GPU
    Series<Device<double>> dx(n), dy(n);
    host_to_device(dx, x);
    f(dy, dx.view());
    // Compare
    ASSERT_CLOSE_D(dy, y);
  }
}

template<class F> void binary_test(const int64_t lo, const int64_t hi, F&& f) {
  mt19937 mt(7);
  for (int n = lo; n <= hi; n++) {
    // CPU
    const Series<double> x = random_series(mt, n);
    const Series<double> y = random_series(mt, n);
    Series<double> z(n);
    f(z, x.view(), y.view());
    // GPU
    Series<Device<double>> dx(n), dy(n), dz(n);
    host_to_device(dx, x);
    host_to_device(dy, y);
    f(dz, dx.view(), dy.view());
    // Compare
    ASSERT_CLOSE_D(dz, z);
  }
}

#define SLOOP(lo, hi) for (int s = lo; s <= hi; s++)
#define CONSTANT(lo, hi, exp) constant_test(lo, hi, [=](auto& x) { exp; })
#define UNARY(lo, hi, exp) unary_test(lo, hi, [=](auto& y, const auto x) { exp; })
#define BINARY(lo, hi, exp) binary_test(lo, hi, [=](auto& z, const auto x, const auto y) { exp; })

TEST(set_scalar) { CONSTANT(0, 7, x.set_scalar(1, 7.5)); }
TEST(assign_series) { UNARY(0, 7, y = x); }

TEST(add_int) { UNARY(0, 7, y = x; y += 3); }
TEST(add_scalar) { UNARY(0, 7, y = x; y += 3.5); }
TEST(sub_int) { UNARY(0, 7, y = x; y -= 3); }
TEST(sub_scalar) { UNARY(0, 7, y = x; y -= 3.5); }

template<int s> void addsub_test() {
  for (const int k : {0, 2}) {
    BINARY(0, 32, z = x; z.template high_addsub<s>(k, y));
    UNARY(0, 32, y.set_counts(x.known(), 0); y.high_addsub<s>(k, x));
    BINARY(0, 32, z = x; auto w = y.copy(); w.set_counts(y.known(), min(4, y.known())); z.high_addsub<s>(k, w));
  }
}
TEST(add_series) { addsub_test<1>(); }
TEST(sub_series) { addsub_test<-1>(); }

TEST(mul) {
  BINARY(0, 128, z = mul(x, y));
  BINARY(0, 128, z = x; z = mul(z, y));
  BINARY(0, 128, z = y; z = mul(x, z));
}

TEST(sqr) {
  UNARY(0, 128, y = sqr(x));
  UNARY(0, 128, y = x; y = sqr(y));
}

TEST(mul1p) {
  SLOOP(1, 3) {
    BINARY(0, 128, z = mul1p(x, y, s));
    BINARY(0, 128, z = y; z = mul1p(x, z, s));
  }
}

TEST(inv) { UNARY(0, 128, y = inv(x)); }
TEST(div) { BINARY(0, 128, z = div(x, y)); }
TEST(inv1p) { SLOOP(1, 3) UNARY(0, 128, y = inv1p(x, s)); }
TEST(div1p) { SLOOP(1, 3) BINARY(0, 128, z = div1p(z, x, y, s));}
TEST(log) { UNARY(0, 128, y = log(x)); }
TEST(log1p) { SLOOP(1, 3) UNARY(0, 128, y = log1p(x, s)); }
TEST(derivative_shift) { SLOOP(0, 3) UNARY(0, 128, y = derivative_shift(x, s));}
TEST(integral_shift) { SLOOP(0, 3) UNARY(0, 128, y = integral_shift(x, s)); }
TEST(exp) { UNARY(0, 128, y = exp(x)); }
TEST(expm1) { SLOOP(1, 3) for (const int a : {1,-1}) UNARY(0, 128, y = expm1(x, a, s));}
TEST(log1p_exp) { SLOOP(1, 3) UNARY(0, 128, y = log1p_exp(x, s)); }
TEST(ldexp) { SLOOP(-3, 3) UNARY(0, 128, y = ldexp(x, s)); }

}  // namespace
}  // namespace mandelbrot
