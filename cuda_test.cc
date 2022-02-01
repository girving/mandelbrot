// CUDA tests

#include "cuda.h"
#include "gtest/gtest.h"
namespace mandelbrot {
namespace {

template<class S> __global__ void add(const int64_t n, S* x, S* y) {
  for (int64_t i = 0; i < n; i++)
    y[i] = x[i] + y[i];
}

TEST(cuda, add) {
  const int64_t n = 1<<20;
  typedef float S;
  vector<S> x(n), y(n);  // CPU memory

  // Initialize x and y on the host
  for (int64_t i = 0; i < n; i++) {
    x[i] = 1;
    y[i] = 2;
  }
  
  // Copy to GPU
  const Stream s;
  auto dx = alloc<S>(n, s);
  auto dy = alloc<S>(n, s);
  cuda_check(cudaMemcpyAsync(dx.get(), x.data(), n*sizeof(S), cudaMemcpyHostToDevice, s));
  cuda_check(cudaMemcpyAsync(dy.get(), x.data(), n*sizeof(S), cudaMemcpyHostToDevice, s));

  // Add on GPU
  add<<<1, 1, 0, s>>>(n, dx.get(), dy.get());
  dx.reset();

  // Copy back, then wait for stream to synchronize
  cuda_check(cudaMemcpyAsync(y.data(), dy.get(), n*sizeof(S), cudaMemcpyDeviceToHost, s));
  dy.reset();
  s.sync();
  
  // Verify results
  for (int64_t i = 0; i < n; i++)
    ASSERT_EQ(y[i], 3);
}

}  // namespace
}  // namespace mandelbrot
