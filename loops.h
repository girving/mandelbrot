// Cuda/host loops
#pragma once

#include "cutil.h"
#include "debug.h"
#include "preprocessor.h"
#include "print.h"
#include <omp.h>
#include <optional>
#include <type_traits>
namespace mandelbrot {

using std::is_signed_v;
using std::optional;
using std::tuple;

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
    _Pragma("omp parallel for") \
    for (int64_t i = 0; i < n; i++) { body } \
  } \
  template<class... Args> static inline void name(const int64_t n, Args&&... xs) { \
    if (!n) return; \
    if constexpr ((... || is_device<Args>)) \
      INVOKE_GRID_STRIDE_LOOP(name##_device, n, undevice(xs)...); \
    else \
      name##_host(n, std::forward<Args>(xs)...); \
  }

// Define a serial function on CPU and GPU (not a loop, but meh).
// This is for reducing the number of total kernel invocations in base cases.
#define DEF_SERIAL(name, args, body) \
  IF_CUDA(template<class S> __global__ static void name##_device(UNPAREN args) { body }) \
  template<class S> static void name##_host(UNPAREN args) { body } \
  template<class... Args> static inline void name(Args&&... xs) { \
    if constexpr ((... || is_device<Args>)) \
      CUDA_OR_DIE(name##_device<<<1, 1, 0, stream()>>>(undevice(xs)...)); \
    else \
      name##_host(std::forward<Args>(xs)...); \
  }

// Chop a loop into [start,end) chunks
tuple<int64_t,int64_t> partition_loop(const int64_t steps, const int threads, const int thread);

// OpenMP reductions that assume only associativity.
// Formally, if (reduce(y, a), reduce(y, b)) is equivalent to (reduce(a, b), reduce(y, a)), then
// this routine is equivalent to:
//   for (int64_t i = 0; i < n; i++)
//     reduce(y, map(i));
template<class Y, class R, class M> void map_reduce(Y& y, R&& reduce, M&& map, const int64_t n) {
  vector<optional<Y>> partials;
  #pragma omp parallel
  {
    const int threads = omp_get_num_threads();
    const int thread = omp_get_thread_num();
    const auto [start, end] = partition_loop(n, threads, thread);
    if (start < end) {
      #pragma omp critical
      {
        partials.resize(threads);
      }
      auto& p = partials[thread];
      for (int64_t i = start; i < end; i++) {
        auto fx = map(i);
        if (i == start) p = move(fx);
        else reduce(*p, fx);
      }
    }
  }
  for (const auto& t : partials)
    if (t)
      reduce(y, *t);
}

}  // namespace mandelbrot
