// Orbits over Expansion<2> (double-double, about 106 bits), as a precision reference for double
//
// Orbit<Expansion<2>> runs double's algorithm with double's tolerances, so paired differences against double
// (escape_tree --prec comparedd) measure double's rounding alone.
#pragma once

#include "expansion_arith.h"
#include "orbit.h"
namespace mandelbrot {

// Same algorithm and tolerances as double
template<> struct OrbitTol<Expansion<2>> : public OrbitTol<double> {};

// The step's fused multiply-add: exact products and sums to about 2^-104
__host__ __device__ inline Expansion<2> orbit_fma(const Expansion<2> a, const Expansion<2> b, const Expansion<2> c) {
  return a * b + c;
}

}  // namespace mandelbrot
