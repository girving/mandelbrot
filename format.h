// A few overloads on top of tinyformat::format
#pragma once

#include "tinyformat.h"
#include "span.h"
#include <vector>
namespace mandelbrot {

using std::string;
using tinyformat::format;

static inline string format() { return string(); }
static inline const string& format(const string& s) { return s; }

static inline string safe(const double x) { return format("%.17g", x); }

}  // namespace mandelbrot
namespace std {
template<class T> ostream& operator<<(ostream& out, const vector<T>& x) {
  return out << span<const T>(x);
}
}  // namespace std
