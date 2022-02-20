// Loops tests

#include "array.h"
#include "loops.h"
#include "tests.h"
#include <random>
namespace mandelbrot {
namespace {

using std::abs;
using std::make_unique;
using std::mt19937;
using std::uniform_real_distribution;

DEF_LOOP(some_loop, n, i, (S* y, const S* x, const int a),
  y[i] = a - sqr(x[i]);)

TEST(loop) {
  typedef float S;
  const int n = 1231;
  const int a = 3;
  mt19937 mt(7);
  const S tol = 1e-5;

  // Host
  Array<S> x(n), y(n);
  for (int i = 0; i < n; i++)
    x[i] = uniform_real_distribution<S>(-1, 1)(mt);
  some_loop(n, y.data(), x.data(), a);
  for (int i = 0; i < n; i++)
    ASSERT_LE(abs(y[i] - (a - sqr(x[i]))), tol);

  // Device
  IF_CUDA(
    Array<Device<S>> dx(n), dy(n);
    host_to_device<S>(dx, x);
    some_loop(n, dy.data(), dx.data(), a);
    Array<S> hy(n);
    device_to_host<S>(hy, dy);
    for (int i = 0; i < n; i++)
      ASSERT_LE(abs(y[i] - hy[i]), tol);
  )
}

DEF_SERIAL(some_serial, (S* x, const S a), *x = sqr(a);)

TEST(serial) {
  typedef float S;
  const S a = 3;

  // Host
  S x = 0;
  some_serial(&x, a);
  ASSERT_EQ(x, sqr(a));

  // Device
  IF_CUDA(
    Array<Device<S>> dx(1);
    some_serial(dx.data(), a);
    Array<S> hx(1);
    device_to_host<S>(hx, dx);
    ASSERT_EQ(hx[0], sqr(a));
  )
}

template<class T> struct Box : public unique_ptr<T> {
  Box(const T& x) : unique_ptr<T>(new T(x)) {}
};
TEST(map_reduce) {
  const auto reduce = [](Box<string>& y, const Box<string>& fx) {
    if (y->size()) *y += ", ";
    *y += *fx;
  };
  const auto map = [](const int i) {
    return Box<string>(str(i));
  };
  vector<int> xs;
  for (int n = 0; n <= 1024; n++) {
    Box<string> s("");
    map_reduce(s, reduce, map, n);
    ASSERT_TRUE(s);
    ASSERT_EQ("[" + *s + "]", str(xs));
    xs.push_back(n);
  }
}

}  // namespace
}  // namespace mandelbrot
