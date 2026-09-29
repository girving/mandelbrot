// Double-double arithmetic for orbits, as a precision reference for double
//
// DoubleDouble is an unevaluated sum hi + lo with |lo| ≤ ulp(hi) / 2, so about 106 bits.  Addition uses
// two-sum and multiplication an explicit fma two-product (exact, since FMA contraction is off in all our
// builds).  Orbits over it use double's algorithm and tolerances (OrbitTol), so paired differences against
// double measure double's rounding alone.
#pragma once

#include "cutil.h"
#include "orbit.h"
#include <cmath>
namespace mandelbrot {

struct DoubleDouble {
  double hi, lo;

  DoubleDouble() = default;
  __host__ __device__ DoubleDouble(const double a) : hi(a), lo(0) {}
  __host__ __device__ explicit operator double() const { return hi + lo; }

  __host__ __device__ static DoubleDouble two_sum(const double a, const double b) {
    const double s = a + b, bb = s - a;
    DoubleDouble r;
    r.hi = s; r.lo = (a - (s - bb)) + (b - bb);
    return r;
  }
  __host__ __device__ static DoubleDouble quick_two_sum(const double a, const double b) {
    const double s = a + b;
    DoubleDouble r;
    r.hi = s; r.lo = b - (s - a);
    return r;
  }

  __host__ __device__ DoubleDouble operator-() const { DoubleDouble r; r.hi = -hi; r.lo = -lo; return r; }
  __host__ __device__ friend DoubleDouble operator+(const DoubleDouble a, const DoubleDouble b) {
    DoubleDouble s = two_sum(a.hi, b.hi);
    const DoubleDouble t = two_sum(a.lo, b.lo);
    s.lo += t.hi;
    s = quick_two_sum(s.hi, s.lo);
    s.lo += t.lo;
    return quick_two_sum(s.hi, s.lo);
  }
  __host__ __device__ friend DoubleDouble operator-(const DoubleDouble a, const DoubleDouble b) { return a + (-b); }
  __host__ __device__ friend DoubleDouble operator*(const DoubleDouble a, const DoubleDouble b) {
    const double p = a.hi * b.hi, e = std::fma(a.hi, b.hi, -p);
    return quick_two_sum(p, e + (a.hi * b.lo + a.lo * b.hi));
  }
  __host__ __device__ friend DoubleDouble operator/(const DoubleDouble a, const DoubleDouble b) {
    // Two Newton-style correction steps on the double quotient
    const double q1 = a.hi / b.hi;
    const DoubleDouble r = a - b * DoubleDouble(q1);
    const double q2 = r.hi / b.hi;
    const DoubleDouble r2 = r - b * DoubleDouble(q2);
    const double q3 = r2.hi / b.hi;
    return quick_two_sum(q1, q2) + DoubleDouble(q3);
  }
  __host__ __device__ friend DoubleDouble orbit_fma(const DoubleDouble a, const DoubleDouble b, const DoubleDouble c) {
    return a * b + c;
  }
  __host__ __device__ DoubleDouble& operator+=(const DoubleDouble b) { return *this = *this + b; }
  __host__ __device__ DoubleDouble& operator-=(const DoubleDouble b) { return *this = *this - b; }
  __host__ __device__ friend bool operator<(const DoubleDouble a, const DoubleDouble b) {
    return a.hi < b.hi || (a.hi == b.hi && a.lo < b.lo);
  }
  __host__ __device__ friend bool operator>(const DoubleDouble a, const DoubleDouble b) { return b < a; }
  __host__ __device__ friend bool operator<=(const DoubleDouble a, const DoubleDouble b) { return !(b < a); }
  __host__ __device__ friend bool operator>=(const DoubleDouble a, const DoubleDouble b) { return !(a < b); }
  __host__ __device__ friend bool operator==(const DoubleDouble a, const DoubleDouble b) {
    return a.hi == b.hi && a.lo == b.lo;
  }
  __host__ __device__ friend DoubleDouble operator+(const double a, const DoubleDouble b) { return DoubleDouble(a) + b; }
  __host__ __device__ friend DoubleDouble operator-(const double a, const DoubleDouble b) { return DoubleDouble(a) - b; }
  __host__ __device__ friend DoubleDouble operator*(const double a, const DoubleDouble b) { return DoubleDouble(a) * b; }
  __host__ __device__ friend DoubleDouble operator+(const DoubleDouble a, const double b) { return a + DoubleDouble(b); }
  __host__ __device__ friend DoubleDouble operator-(const DoubleDouble a, const double b) { return a - DoubleDouble(b); }
  __host__ __device__ friend DoubleDouble operator*(const DoubleDouble a, const double b) { return a * DoubleDouble(b); }
};

// Same algorithm and tolerances as double
template<> struct OrbitTol<DoubleDouble> : public OrbitTol<double> {};

}  // namespace mandelbrot
