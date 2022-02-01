// Shim if std::span is missing
#pragma once

#include <cassert>
#include <ostream>
#include <vector>

#if __has_include(<span>)

#define SPAN_NAMESPACE std
#include <span>
namespace mandelbrot {
using std::span;
using std::vector;
}

#else  // No <span>, so roll our own

#define SPAN_NAMESPACE mandelbrot
namespace mandelbrot {

using std::vector;

template<class T> class span {
  T* x;
  size_t n;
public:

  span(T* x, size_t n) : x(x), n(n) {}
  template<class S,size_t k> span(S (&y)[k]) : x(y), n(k) {}
  template<class S> span(span<S> y) : x(y.data()), n(y.size()) {}
  template<class S> span(vector<S>& y) : x(y.data()), n(y.size()) {}
  template<class S> span(const vector<S>& y) : x(y.data()), n(y.size()) {}

  size_t size() const { return n; }
  T* data() const { return x; }
  T& operator[](const size_t i) const { assert(i < n); return x[i]; }
};

}  // namespace mandelbrot
#endif  // __has_include(<span>)

namespace std {
template<class T> ostream& operator<<(ostream& out, SPAN_NAMESPACE::span<T> x) {
  out << '[';
  for (size_t i = 0; i < x.size(); i++) {
    if (i) out << ", ";
    out << x[i];
  }
  return out << ']';
}
}  // namespace std
