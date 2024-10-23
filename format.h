// A few overloads on top of tinyformat::format
#pragma once

#include "tinyformat.h"
#include "span.h"
#include <vector>
namespace tinyformat {

static inline std::string format() { return std::string(); }
static inline const std::string& format(const std::string& s) { return s; }

} namespace mandelbrot {

using std::string;

static inline string safe(const double x) { return tfm::format("%.17g", x); }

}  // namespace mandelbrot
namespace std {
template<class T> ostream& operator<<(ostream& out, const vector<T>& x) {
  return out << span<const T>(x);
}
}  // namespace std
