// Complex numbers with custom real types
#pragma once

#include <iostream>
#include <cmath>
namespace mandelbrot {

using std::abs;
using std::ostream;

// Double helpers
static inline double sqr(double x) { return x * x; }
static inline double twice(double x) { return x + x; }
static inline double inv(double x) { return 1 / x; }

// We need our own complex template to handle custom scalars.
// Ours has very few operations.
template<class S> struct Complex {
  typedef S Real;

  S r, i;

  Complex() : r(0), i(0) {}
  explicit Complex(const S& r) : r(r), i(0) {}
  Complex(const S& r, const S& i) : r(r), i(i) {}

  Complex operator-() const { return Complex(-r, -i); }
  Complex operator+(const Complex& z) const { return Complex(r + z.r, i + z.i); }
  Complex operator-(const Complex& z) const { return Complex(r - z.r, i - z.i); }
  void operator+=(const Complex& z) { r += z.r; i += z.i; }
  void operator-=(const Complex& z) { r -= z.r; i -= z.i; }
  Complex operator*(const Complex& z) const { return Complex(r*z.r - i*z.i, r*z.i + i*z.r); }
  friend Complex operator*(const S a, const Complex& z) { return Complex(a*z.r, a*z.i); }
  friend Complex sqr(const Complex& z) { return Complex(sqr(z.r) - sqr(z.i), twice(z.r*z.i)); }
  friend Complex conj(const Complex& z) { return Complex(z.r, -z.i); }
  friend Complex left(const Complex& z) { return Complex(-z.i, z.r); }  // iz
  friend Complex right(const Complex& z) { return Complex(z.i, -z.r); }  // -iz
  friend Complex twice(const Complex& z) { return Complex(twice(z.r), twice(z.i)); }

  friend ostream& operator<<(ostream& out, const Complex& z) {
    out << z.r;
    if (copysign(S(1), z.i) > 0) out << '+';
    return out << z.i << 'j';
  }
};

// Diagonal complex scaling by a(1 +- i).  2 adds, 2 muls.
template<int sign, class S> static inline Complex<S> diag(const S& a, const Complex<S>& z) {
  static_assert(sign == 1 || sign == -1);
  if constexpr (sign == 1) return a * Complex<S>(z.r - z.i, z.i + z.r);
  else return a * Complex<S>(z.r + z.i, z.i - z.r);
}

static inline double abs(const Complex<double> z) {
  return hypot(z.r, z.i);
}

}  // namespace mandelbrot
