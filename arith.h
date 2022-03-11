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
  __host__ __device__ static inline S bound(S x) { return abs(x); } \
  __host__ __device__ static inline S fma(const S x, const S y, const S s) { return __builtin_fma(x, y, s); }
ARITH(float)
ARITH(double)
#undef ARITH

// Integer overload to make generated code happy
__host__ __device__ static inline int twice(int x) { return x << 1; }

// relu(x) = max(0, x)
template<class T> static inline T relu(const T& x) {
  return max(x, T(0));
}

static inline int exponent(const double x) {
  int e;
  frexp(x, &e);
  return e;
}

template<class I> static inline I exact_div(const I a, const I b) {
  const I r = a / b;
  slow_assert(r * b == a);
  return r;
}

}  // namespace mandelbrot
