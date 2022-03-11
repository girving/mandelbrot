// Fast Fourier Transforms
#pragma once

#include "complex.h"
#include "device.h"
#include "span.h"
#include <cstdint>
namespace mandelbrot {

using std::add_const_t;
using std::function;

// Complex size-n fft of x, with n = y.size() and x zero extended to 2n reals.
template<class S> void fft(span<Complex<S>> y, span<const S> x);

// Complex to real inverse fft.
// WARNING: The input array y is scrambled.
template<class S> void ifft(span<S> x, span<Complex<S>> y);

// Real to complex size-n fft of x, with n/2 = y.size() and x zero padded to size n.
// Only half of the n-output FFT is computed, taking advantage of Hermitian symmetry.
// y[0].r and y[0].i correspond to y[0] and y[n] of the full FFT.
template<class S> void rfft(span<Complex<S>> y, span<const S> x);

// Complex to real inverse fft.
// WARNING: The input array y is scrambled.
template<class S> void irfft(span<S> x, span<Complex<S>> y);

// Real to complex size-n shifted fft of x, with n/2 = y.size() and x zero padded to size n.
// This corresponds to the odd entries of the size 2n FFT.  Since the odd entries all come
// in complex conjugate pairs (unlike y[0] for rfft above), the computation is more uniform.
// Only half of the n-output FFT is computed, taking advantage of Hermitian symmetry.
template<class S> void srfft(span<Complex<S>> y, span<const S> x);

// Complex to real inverse shifted fft.
// WARNING: The input array y is scrambled.
// { srfft(y, x); isrfft(x, y); } is x *= n/2 = y.size()
template<class S> void isrfft(span<S> x, span<Complex<S>> y);

// z[:] = x[:] * y[:] as polynomials.
// Aliasing is allowed.
template<class T> void fft_mul(span<T> z, span<add_const_t<T>> x, span<add_const_t<T>> y);

// z[:] += x[:] * y[:] as polynomials.
// Aliasing is allowed.  middle() is called between when x,y are used and when z is used.
template<class T> void fft_addmul(span<T> z, span<add_const_t<T>> x, span<add_const_t<T>> y,
                                  const function<void()>& middle);

// y[:] = x[:]^2 as polynomials.
// Aliasing is allowed.
template<class T> void fft_sqr(span<T> y, span<add_const_t<T>> x);

}  // namespace mandelbrot
