// span<T>
#pragma once

#include <cassert>
#include <ostream>
#include <span>
#include <vector>
namespace mandelbrot {

using std::span;
using std::vector;

}  // namespace mandelbrot
namespace std {
template<class T> ostream& operator<<(ostream& out, span<T> x) {
  out << '[';
  for (size_t i = 0; i < x.size(); i++) {
    if (i) out << ", ";
    out << x[i];
  }
  return out << ']';
}
}  // namespace std
