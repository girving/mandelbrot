// CUDA tests

#include "cutil.h"
#include "tests.h"
#include <cuda.h>
namespace mandelbrot {
namespace {

template<class S> __global__ void add(const int64_t n, S* x, S* y) {
  for (int64_t i = 0; i < n; i++)
    y[i] = x[i] + y[i];
}

TEST(add) {
  const int64_t n = 1<<20;
  typedef float S;
  vector<S> x(n), y(n);  // CPU memory

  // Initialize x and y on the host
  for (int64_t i = 0; i < n; i++) {
    x[i] = 1;
    y[i] = 2;
  }
  
  // Copy to GPU
  const auto s = stream();
  auto dx = alloc<S>(n);
  auto dy = alloc<S>(n);
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

}  // namespace
}  // namespace mandelbrot
