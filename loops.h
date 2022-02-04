// Cuda/host loops
#pragma once

#include "cutil.h"
#include "debug.h"
#include "preprocessor.h"
#include <type_traits>
namespace mandelbrot {

using std::is_signed_v;

// For now, assume we fit in int32_t
template<class I> __device__ static inline int grid_stride_loop_size(const I n) {
  static_assert(is_signed_v<I>);
  return int(n);
}

// Define a grid stride loop
#define GRID_STRIDE_LOOP(n, i) \
  for (int _n = grid_stride_loop_size(n), \
           _stride = blockDim.x * gridDim.x, \
           i = blockIdx.x * blockDim.x + threadIdx.x; \
       i < _n; i += _stride)

// Call a one-dimensional grid-stride loop
#define INVOKE_GRID_STRIDE_LOOP(name, n, ...) CUDA_OR_DIE(({ \
  const int _n = (n);  /* For now, assume we fit in int32_t */ \
  name<<<32*num_sms(), 256>>>(_n, __VA_ARGS__); }))

// Define 1D loop functions on CPU and GPU
#define DEF_LOOP(name, n, i, args, body) \
  IF_CUDA(template<class S> __global__ static void name##_device(const int n, UNPAREN args) { \
    GRID_STRIDE_LOOP(n, i) { body } \
  }) \
  template<class S> static void name##_host(const int n, UNPAREN args) { \
    for (int64_t i = 0; i < n; i++) { body } \
  } \
  template<class... Args> static inline void name(const int64_t n, Args&&... xs) { \
    if (!n) return; \
    if constexpr ((... || is_device<Args>)) \
      INVOKE_GRID_STRIDE_LOOP(name##_device, n, undevice(xs)...); \
    else \
      name##_host(n, std::forward<Args>(xs)...); \
  }

}  // namespace mandelbrot
