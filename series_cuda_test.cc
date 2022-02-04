// Series CUDA tests

#include "series.h"
#include "tests.h"
#include <random>
namespace mandelbrot {
namespace {

using std::min;
using std::mt19937;
using std::uniform_real_distribution;

void assert_rel(SeriesView<const Device<double>> dz, SeriesView<const double> z, const double tol) {
  Series<double> hz(dz.nonzero());
  device_to_host(hz, dz);
  const auto e = error(hz, z, true);
  ASSERT_LE(e, tol) << format("n %d, e %g\ndz %g\nz %g", z.nonzero(), e, hz, z);
}

Series<double> random_series(mt19937& mt, const int n) {
  Series<double> x(n);
  x.set_counts(n, n);
  for (int i = 0; i < n; i++)
    x[i] = uniform_real_distribution<double>(-1, 1)(mt);
  return x;
}

template<class F> void constant_test(const int64_t lo, const int64_t hi, F&& f, const double tol = 1e-14) {
  for (int n = lo; n <= hi; n++) {
    // CPU
    Series<double> x(1);
    f(x);
    // GPU
    Series<Device<double>> dx(1);
    f(dx);
    // Compare
    assert_rel(dx, x, tol);
  }
}

struct Nop { void operator()(auto& x) const {} };

template<class F,class P=Nop> void
unary_test(const int64_t lo, const int64_t hi, F&& f, const double tol = 1e-14, P project = P()) {
  for (int n = lo; n <= hi; n++) {
    mt19937 mt(7);
    // CPU
    Series<double> x = random_series(mt, n);
    project(x);
    Series<double> y(n);
    f(y, x.view());
    // GPU
    Series<Device<double>> dx(n), dy(n);
    host_to_device(dx, x);
    f(dy, dx.view());
    // Compare
    assert_rel(dy, y, tol);
  }
}

template<class F> void binary_test(const int64_t lo, const int64_t hi, F&& f, const double tol = 1e-14) {
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
    assert_rel(dz, z, tol);
  }
}

#define SLOOP(lo, hi) for (int s = lo; s <= hi; s++)
#define CONSTANT(lo, hi, exp) constant_test(lo, hi, [=](auto& x) { exp; })
#define UNARY(lo, hi, exp, ...) unary_test(lo, hi, [=](auto& y, const auto x) { exp; } __VA_OPT__(,) __VA_ARGS__)
#define BINARY(lo, hi, exp, ...) \
  binary_test(lo, hi, [=](auto& z, const auto x, const auto y) { exp; } __VA_OPT__(,) __VA_ARGS__)

TEST(set_scalar) { CONSTANT(0, 7, x.set_scalar(1, 7.5)); }
TEST(assign_series) { UNARY(0, 7, y = x); }

TEST(add_int) { UNARY(1, 7, y = x; y += 3); }
TEST(add_scalar) { UNARY(1, 7, y = x; y += 3.5); }
TEST(sub_int) { UNARY(1, 7, y = x; y -= 3); }
TEST(sub_scalar) { UNARY(2, 7, y = x; y -= 3.5); }

TEST(addsub) {
  for (const int s : {1, -1}) {
    for (const int k : {0, 2}) {
      BINARY(0, 32, z = x; high_addsub(z, s, k, y));
      UNARY(0, 32, y.set_counts(x.known(), 0); high_addsub(y, s, k, x));
      BINARY(0, 32, z = x; auto w = y.copy(y.limit()); w.set_counts(y.known(), min(4l, y.known()));
                    high_addsub(z, s, k, w));
    }
  }
}

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

static auto constant(const int a) {
  return [a](Series<double>& x) { if (x.nonzero()) x[0] = a; };
}

TEST(inv) { UNARY(0, 32, y = inv(x), 1e-9); }
TEST(div) { BINARY(0, 8, z = div(x, y), 1e-7); }
TEST(inv1p) { SLOOP(1, 3) UNARY(0, 32, y = inv1p(x, s), 1e-9); }
TEST(div1p) { SLOOP(1, 3) BINARY(0, 32, z = div1p(x, y, s), 1e-10);}
TEST(log) { UNARY(0, 32, y = log(x), 1e-10, constant(1)); }
TEST(log1p) { SLOOP(1, 3) UNARY(0, 32, y = log1p(x, s), 1e-10); }
TEST(derivative_shift) { SLOOP(0, 3) UNARY(0, 32, y = derivative_shift(x, s));}
TEST(integral_shift) { SLOOP(0, 3) UNARY(0, 32, y = integral_shift(x, s)); }
TEST(exp) { UNARY(0, 32, y = exp(x), 1e-10, constant(0)); }
TEST(expm1) { SLOOP(1, 3) for (const int a : {1,-1}) UNARY(0, 32, y = expm1(x, a, s), 1e-10);}
TEST(log1p_exp) { SLOOP(1, 3) UNARY(0, 16, y = log1p_exp(x, s), 1e-7, constant(0)); }
TEST(ldexp) { SLOOP(-3, 3) UNARY(0, 32, y = ldexp(x, s), 1e-9); }

}  // namespace
}  // namespace mandelbrot
