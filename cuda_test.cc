// CUDA tests

#include "array.h"
#include "cutil.h"
#include "loops.h"
#include "tests.h"
#include <cuda.h>
#include <random>
namespace mandelbrot {
namespace {

using std::abs;
using std::mt19937;
using std::uniform_real_distribution;

template<class S> __global__ void add(const int64_t n, S* x, S* y) {
  for (int64_t i = 0; i < n; i++)
    y[i] = x[i] + y[i];
}

TEST(add_serial) {
  const int64_t n = 1024;
  typedef float S;
  vector<S> x(n), y(n);  // CPU memory

  // Initialize x and y on the host
  for (int64_t i = 0; i < n; i++) {
    x[i] = 1;
    y[i] = 2;
  }
  
  // Copy to GPU
  const auto s = stream();
  Array<Device<S>> dx(n);
  Array<Device<S>> dy(n);
  host_to_device<S>(dx, x);
  host_to_device<S>(dy, y);

  // Add on GPU
  add<<<1, 1, 0, s>>>(n, device_get(dx), device_get(dy));
  dx.clear();

  // Copy back, then wait for stream to synchronize
  device_to_host<S>(y, dy);
  dy.clear();
  cuda_sync();
  
  // Verify results
  for (int64_t i = 0; i < n; i++)
    ASSERT_EQ(y[i], 3);
}

__global__ void sqr1a(const int n, float* y, const float* x, const int a) {
  GRID_STRIDE_LOOP(n, i) y[i] = sqr(x[i] + a);
}

TEST(grid_stride_loop) {
  typedef float S;
  const int max_n = 1 << 13;
  const int sentinels = 256;

  // Prepare memory
  vector<S> x(max_n + sentinels);
  vector<S> y(max_n + sentinels, -1);
  for (int i = 0; i < int(x.size()); i++)
    x[i] = i;
  Array<Device<S>> dx(x.size());
  Array<Device<S>> dy(y.size());
  host_to_device<S>(dx, x);
  host_to_device<S>(dy, y);

  // Compute
  for (int n = 0; n <= max_n; n++) {
    const int ns = n + sentinels;
    const int a = n + 1;
    INVOKE_GRID_STRIDE_LOOP(sqr1a, n, device_get(dy), device_get(dx), a);
    device_to_host<S>(span<S>(y).first(ns), dy.span().first(ns));
    for (int i = 0; i < ns; i++)
      ASSERT_EQ(y[i], i < n ? sqr(x[i] + a) : -1) << format("n %d, ns %d, a %d, i %d", n, ns, a, i);
  }
}

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
  Array<Device<S>> dx(n), dy(n);
  host_to_device<S>(dx, x);
  some_loop(n, dy.data(), dx.data(), a);
  Array<S> hy(n);
  device_to_host<S>(hy, dy);
  for (int i = 0; i < n; i++)
    ASSERT_LE(abs(y[i] - hy[i]), tol);
}

}  // namespace
}  // namespace mandelbrot
