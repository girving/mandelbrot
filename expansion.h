// Higher precision floating point numbers via expansion arithmetic
#pragma once

#include "arith.h"
#include "span.h"
#include <cmath>
#include <ostream>
#include <random>
namespace mandelbrot {

// References:
//   Shewchuk 1997, Adaptive precision floating-point arithmetic and fast robust geometric predicates
//     https://people.eecs.berkeley.edu/~jrs/papers/robustr.pdf
//   Lu et al. 2020, Supporting extended precision on graphics processors
//     https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.707.4321&rep=rep1&type=pdf
//   Collange et al. 2016, Parallel floating-point expansions for extended-precision GPU computations
//     https://hal.archives-ouvertes.fr/hal-01298206/document

struct Arb;
struct Arf;
using std::is_same_v;
using std::mt19937;
using std::string;
using std::string_view;
using std::ostream;

struct Nonoverlap {};
constexpr Nonoverlap nonoverlap;

template<int n_> struct Expansion {
  struct Unusable {};
  static constexpr int n = n_;
  static_assert(n >= 2);
  double x[n];  // Ulp-nonoverlapping in decreasing order of magnitude

  __host__ __device__ Expansion() : x{0} {}
  __host__ __device__ Expansion(Unusable*) : x{0} {}  // Allow construction from literal 0

  // These promise that the nonoverlapping invariant already holds
  __host__ __device__ Expansion(double x0, double x1, Nonoverlap) : x{x0, x1} { static_assert(n == 2); }
  __host__ __device__ Expansion(double x0, double x1, double x2, Nonoverlap) : x{x0, x1, x2} { static_assert(n == 3); }
  __host__ __device__ Expansion(double x0, double x1, double x2, double x3, Nonoverlap)
      : x{x0, x1, x2, x3} { static_assert(n == 4); }

  // Conversion from smaller types
  __host__ __device__ explicit Expansion(const double a) : x{a, 0} {}
  __host__ __device__ explicit Expansion(const int32_t a) : x{double(a), 0} {}
  __host__ __device__ explicit Expansion(const int64_t a) : x{0} { x[0] = a; x[1] = a - int64_t(x[0]); }
  __host__ __device__ void operator=(const int32_t a) { x[0] = a; for (int i = 1; i < n; i++) x[i] = 0; }

  // In place arithmetic
  __host__ __device__ void operator+=(const int y) { *this = *this + Expansion(y); }
  __host__ __device__ void operator-=(const int y) { *this = *this - Expansion(y); }
  __host__ __device__ void operator+=(const Expansion y) { *this = *this + y; }
  __host__ __device__ void operator-=(const Expansion y) { *this = *this - y; }

  // Multiplication by integers which might be too big to fit in a double
  friend __host__ __device__ Expansion operator*(const int64_t a, const Expansion x) { return Expansion(a) * x; }

  // Division via Newton iteration
  __host__ __device__ Expansion operator/(const Expansion b) const;
  __host__ __device__ Expansion operator/(const int64_t b) const { return *this / Expansion(b); }

  // Componentwise-safe operations
  #define CWISE(op) Expansion y; for (int i = 0; i < n; i++) y.x[i] = op; return y;
  friend __host__ __device__ Expansion half(const Expansion x) { CWISE(mandelbrot::half(x.x[i])) }
  friend __host__ __device__ Expansion twice(const Expansion x) { CWISE(mandelbrot::twice(x.x[i])) }
  friend __host__ __device__ Expansion ldexp(const Expansion x, const int e) { CWISE(std::ldexp(x.x[i], e)) }
  #undef CWISE

  bool operator==(const int a) const {
    if (x[0] != a)
      return false;
    for (int i = 1; i < n; i++)
      if (x[i])
        return false;
    return true;
  }

  explicit operator double() const {
    double s = x[n-1];
    for (int i = n-2; i >= 0; i--)
      s += x[i];
    return s;
  }

  std::span<const double> span() const {
    return std::span<const double>(x, size_t(n));
  }

  explicit operator bool() const {
    for (int i = 0; i < n; i++)
      if (x[i])
        return true;
    return false;
  }

  bool operator==(const Expansion y) const {
    for (int i = 0; i < n; i++)
      if (x[i] != y.x[i])
        return false;
    return true;
  }

  bool operator!=(const Expansion y) const { return !(*this == y); }

  // These are slow
  Arb arb(const int prec) const;
  explicit Expansion(const string& s);
  explicit Expansion(string_view s);
};

template<int n> int sign(const Expansion<n> x);
template<int n> Expansion<n> abs(const Expansion<n> x);

// For now, don't try to optimize sqr further
template<int n> __host__ __device__ Expansion<n> sqr(const Expansion<n> x) { return x * x; }

// Dangerous if the exponents are far apart
template<int n> Arb exact_arb(const Expansion<n> x);
template<int n> Arf exact_arf(const Expansion<n> x);

// Dangerous if the exponents are far apart, since it calls exact_arf()
template<int n> ostream& operator<<(ostream& out, const Expansion<n> e);

// Print enough digits for exact reconstruction, or fall back to list syntax if we run out of precision
template<int n> string safe(const Expansion<n> x);
template<int n> string maybe_nice_safe(const Expansion<n> x);  // Might return "" if decimal printing fails

// Reciprocal via Newton's method
template<int n> __host__ __device__ Expansion<n> inv(const Expansion<n> x) {
  // Base case
  Expansion<n> y;
  y.x[0] = 1 / x.x[0];

  // Newton step:
  //   1/y = x
  //   f(y) = 1/y - x
  //   f'(y) = -1/y^2
  //   N(y) = y - f(y) / f'(y)
  //        = y - (1/y - x) / (-1/y^2)
  //        = y - y(xy - 1)
  //        = y(2 - xy)
  const Expansion<n> two(2);
  for (int i = 0; i < n; i++)
    y = y*(two - x*y);
  return y;
}

// Division via reciprocal multiplication + Newton's method
template<int n> __host__ __device__ Expansion<n> Expansion<n>::operator/(const Expansion b) const {
  // Compute the inverse and multiply
  const auto inv_b = inv(b);
  Expansion y = *this * inv_b;

  // One more step of Newton refinement:
  //   y = a/b
  //   f(y) = by - a
  //   f'(y) = b
  //   N(y) = y - (b*y - a)/b
  y = y - (b*y - *this) * inv_b;
  return y;
}

// Random numbers for unit tests
template<int n> Expansion<n> random_expansion(mt19937& mt);
template<int n> Expansion<n> random_expansion_near(mt19937& mt, const Expansion<n> x);
template<int n> Expansion<n> random_expansion_with_exponent(mt19937& mt, int e);

// Constants of templated size
template<class S> static inline S constant(double x0, double x1);
template<> inline double constant(double x0, double x1) { return x0; }
template<> inline Expansion<2> constant(double x0, double x1) { return Expansion<2>(x0, x1, nonoverlap); }

}  // namespace mandelbrot
