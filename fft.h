// Fast Fourier Transforms
#pragma once

#include <cstdint>
#include <iostream>
#include <cmath>
namespace mandelbrot {

using std::ostream;

// Double helpers
static inline double sqr(double x) { return x * x; }
static inline double twice(double x) { return x + x; }
static inline double inv(double x) { return 1 / x; }

// We need our own complex template to handle custom scalars.
// Ours has very few operations.
template<class S> struct Complex {
  S r, i;

  Complex() : r(0), i(0) {}
  explicit Complex(const S& r) : r(r), i(0) {}
  Complex(const S& r, const S& i) : r(r), i(i) {}

  Complex operator+(const Complex& z) const { return Complex(r + z.r, i + z.i); }
  Complex operator-(const Complex& z) const { return Complex(r - z.r, i - z.i); }
  void operator+=(const Complex& z) { r += z.r; i += z.i; }
  void operator-=(const Complex& z) { r -= z.r; i -= z.i; }
  Complex operator*(const Complex& z) const { return Complex(r*z.r - i*z.i, r*z.i + i*z.r); }
  friend Complex operator*(const S a, const Complex& z) { return Complex(a*z.r, a*z.i); }
  friend Complex sqr(const Complex& z) { return Complex(sqr(z.r) - sqr(z.i), twice(z.r*z.i)); }

  friend ostream& operator<<(ostream& out, const Complex& z) {
    out << z.r;
    if (z.i >= 0) out << '+';
    return out << z.i << 'j';
  }
};

// Real to complex size-n fft of x[:xn] zero padded to size n
template<class S> void fft(Complex<S>* y, const S* x, const int64_t n, const int64_t xn);

// Complex to real inverse fft, dropping outputs past x[:xn].
// WARNING: The input array y is scrambled.
template<class S> void ifft(S* x, Complex<S>* y, const int64_t n, const int64_t xn);

// z[:n] = x[:n] * y[:n]
template<class S> void fft_mul(S* z, const S* x, const S* y, const int64_t n);

// y[:n] = x[:n]^2
template<class S> void fft_sqr(S* y, const S* x, const int64_t n);

}  // namespace mandelbrot
