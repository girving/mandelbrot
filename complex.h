// Complex numbers with custom real types
#pragma once

#include "arith.h"
#include <iostream>
#include <cmath>
namespace mandelbrot {

using std::ostream;

// We need our own complex template to handle custom scalars.
// Ours has very few operations.
template<class S> struct Complex {
  typedef S Real;

  S r, i;

  __host__ __device__ Complex() : r(0), i(0) {}
  __host__ __device__ explicit Complex(const S& r) : r(r), i(0) {}
  __host__ __device__ Complex(const S& r, const S& i) : r(r), i(i) {}

  __host__ __device__ Complex operator-() const { return Complex(-r, -i); }
  __host__ __device__ Complex operator+(const Complex& z) const { return Complex(r + z.r, i + z.i); }
  __host__ __device__ Complex operator-(const Complex& z) const { return Complex(r - z.r, i - z.i); }
  __host__ __device__ void operator+=(const Complex& z) { r += z.r; i += z.i; }
  __host__ __device__ void operator-=(const Complex& z) { r -= z.r; i -= z.i; }
  __host__ __device__ Complex operator*(const Complex& z) const { return Complex(r*z.r - i*z.i, r*z.i + i*z.r); }
  __host__ __device__ friend Complex operator*(const S a, const Complex& z) { return Complex(a*z.r, a*z.i); }
  __host__ __device__ friend Complex sqr(const Complex& z) { return Complex(sqr(z.r) - sqr(z.i), twice(z.r*z.i)); }
  __host__ __device__ friend Complex conj(const Complex& z) { return Complex(z.r, -z.i); }
  __host__ __device__ friend Complex left(const Complex& z) { return Complex(-z.i, z.r); }  // iz
  __host__ __device__ friend Complex right(const Complex& z) { return Complex(z.i, -z.r); }  // -iz
  __host__ __device__ friend Complex twice(const Complex& z) { return Complex(twice(z.r), twice(z.i)); }
  __host__ __device__ friend Complex half(const Complex& z) { return Complex(half(z.r), half(z.i)); }

  friend ostream& operator<<(ostream& out, const Complex& z) {
    out << z.r;
    if (copysign(S(1), z.i) > 0) out << '+';
    return out << z.i << 'j';
  }
};

// Diagonal complex scaling by a(1 +- i).  2 adds, 2 muls.
template<int sign, class S> __host__ __device__ static inline Complex<S> diag(const S& a, const Complex<S>& z) {
  static_assert(sign == 1 || sign == -1);
  if constexpr (sign == 1) return a * Complex<S>(z.r - z.i, z.i + z.r);
  else return a * Complex<S>(z.r + z.i, z.i - z.r);
}

static inline double abs(const Complex<double> z) {
  return hypot(z.r, z.i);
}

}  // namespace mandelbrot
