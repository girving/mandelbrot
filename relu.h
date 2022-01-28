// relu(x) = max(0, x)
#pragma once

#include <algorithm>
namespace mandelbrot {

using std::max;

template<class T> static inline T relu(const T& x) {
  return max(x, T(0));
}

}  // namespace mandelbrot
