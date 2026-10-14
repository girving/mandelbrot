// Böttcher coefficient octaves: energy density over external angle
#pragma once

#include "span.h"
#include <vector>
namespace mandelbrot {

using std::span;
using std::vector;

// f holds the rows of an f-k{k}.npy array flattened: f[2m] + f[2m+1] = f_m = b_{m-1} as an unevaluated sum.
// Octave j is m ∈ [2^j, 2^(j+1)) with a_m = sqrt(m-1) f_m, and its energy sum a_m^2 is its contribution
// to (1 - area/π).  Returns e with 2^j + 1 entries: e[i] is the energy at θ = i/2^(j+1) with θ and 1-θ
// folded together (b_n is real, so the density is symmetric).  Parseval: sum e = sum a_m^2.
vector<double> octave_energy(span<const double> f, const int j);

// Least squares line y = a + b x
struct Line { double a, b; };
Line fit_line(span<const double> x, span<const double> y);

}  // namespace mandelbrot
