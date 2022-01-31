// A few overloads on top of tinyformat::format
#pragma once

#include "tinyformat.h"
#include <span>
namespace mandelbrot {

using std::string;
using tinyformat::format;

static inline string format() { return string(); }

static inline string safe(const double x) { return format(".17g", x); }

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

template<class T> ostream& operator<<(ostream& out, const vector<T>& x) {
  return out << span<const T>(x);
}

}  // namespace std
