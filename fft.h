// Fast Fourier Transforms
#pragma once

#include "complex.h"
#include <cstdint>
namespace mandelbrot {

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
