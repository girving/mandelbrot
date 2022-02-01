// Basic math functions
#pragma once

#include <algorithm>
#include <cmath>
namespace mandelbrot {

using std::abs;
using std::max;

// Double helpers
static inline double sqr(double x) { return x * x; }
static inline double half(double x) { return 0.5 * x; }
static inline double twice(double x) { return x + x; }
static inline double inv(double x) { return 1 / x; }
static inline double bound(double x) { return abs(x); }

// relu(x) = max(0, x)
template<class T> static inline T relu(const T& x) {
  return max(x, T(0));
}

}  // namespace mandelbrot
