// Basic math functions
#pragma once

#include "cutil.h"
#include <algorithm>
#include <cmath>
namespace mandelbrot {

using std::abs;
using std::max;

// Floating point helpers
#define ARITH(S) \
  __host__ __device__ static inline S sqr(S x) { return x * x; } \
  __host__ __device__ static inline S half(S x) { return S(0.5) * x; } \
  __host__ __device__ static inline S twice(S x) { return x + x; } \
  __host__ __device__ static inline S inv(S x) { return 1 / x; } \
  __host__ __device__ static inline S bound(S x) { return abs(x); }
ARITH(float)
ARITH(double)
#undef ARITH

// relu(x) = max(0, x)
template<class T> static inline T relu(const T& x) {
  return max(x, T(0));
}

}  // namespace mandelbrot
