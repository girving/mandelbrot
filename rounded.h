// Doubles rounded to fewer significant bits, for measuring how escape-time results depend on precision
//
// Rounded<bits> holds a double and rounds the result of every arithmetic operation to `bits` significant bits
// by Veltkamp splitting (c = (2^(53-bits) + 1) a; a ↦ c - (c - a)), which is exact round-to-nearest when FMA
// contraction is off (as in all our builds).  Orbits over Rounded<bits> use double's algorithm and
// tolerances, so any change in classification relative to double comes from rounding alone.
#pragma once

#include "cutil.h"
#include "orbit.h"
namespace mandelbrot {

template<int bits> struct Rounded {
  static_assert(2 <= bits && bits <= 52);
  double v;

  __host__ __device__ static double round(const double a) {
    const double c = double((1ull << (53 - bits)) + 1) * a;
    return c - (c - a);
  }

  Rounded() = default;
  __host__ __device__ Rounded(const double a) : v(round(a)) {}
  __host__ __device__ explicit operator double() const { return v; }

  __host__ __device__ Rounded operator-() const { Rounded r; r.v = -v; return r; }
  __host__ __device__ friend Rounded operator+(const Rounded a, const Rounded b) { return Rounded(a.v + b.v); }
  __host__ __device__ friend Rounded operator-(const Rounded a, const Rounded b) { return Rounded(a.v - b.v); }
  __host__ __device__ friend Rounded operator*(const Rounded a, const Rounded b) { return Rounded(a.v * b.v); }
  __host__ __device__ friend Rounded operator/(const Rounded a, const Rounded b) { return Rounded(a.v / b.v); }
  __host__ __device__ Rounded& operator+=(const Rounded b) { return *this = *this + b; }
  __host__ __device__ Rounded& operator-=(const Rounded b) { return *this = *this - b; }
  __host__ __device__ friend bool operator<(const Rounded a, const Rounded b) { return a.v < b.v; }
  __host__ __device__ friend bool operator>(const Rounded a, const Rounded b) { return a.v > b.v; }
  __host__ __device__ friend bool operator<=(const Rounded a, const Rounded b) { return a.v <= b.v; }
  __host__ __device__ friend bool operator>=(const Rounded a, const Rounded b) { return a.v >= b.v; }
  // Mixed with doubles (literals like 2 * z or 1 + w): the double is rounded first
  __host__ __device__ friend Rounded operator+(const double a, const Rounded b) { return Rounded(a) + b; }
  __host__ __device__ friend Rounded operator-(const double a, const Rounded b) { return Rounded(a) - b; }
  __host__ __device__ friend Rounded operator*(const double a, const Rounded b) { return Rounded(a) * b; }
  __host__ __device__ friend Rounded operator+(const Rounded a, const double b) { return a + Rounded(b); }
  __host__ __device__ friend Rounded operator-(const Rounded a, const double b) { return a - Rounded(b); }
  __host__ __device__ friend Rounded operator*(const Rounded a, const double b) { return a * Rounded(b); }
};

// Same algorithm and tolerances as double
template<int bits> struct OrbitTol<Rounded<bits>> : public OrbitTol<double> {};

}  // namespace mandelbrot
