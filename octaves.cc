// Böttcher coefficient octaves

#include "octaves.h"
#include "debug.h"
#include "fft.h"
#include <cmath>
namespace mandelbrot {

vector<double> octave_energy(span<const double> f, const int j) {
  slow_assert(j >= 1 && int64_t(f.size()) >= int64_t(4) << j, "octave %d needs %d rows", j, int64_t(2) << j);
  const int64_t lo = int64_t(1) << j, n = 2*lo;
  vector<double> x(lo);
  for (int64_t m = lo; m < n; m++)
    x[m - lo] = std::sqrt(double(m - 1)) * (f[2*m] + f[2*m+1]);
  vector<Complex<double>> y(lo);
  rfft<double>(y, x);
  // y[0] packs the θ = 0 and θ = 1/2 entries; others pair with their conjugates at 1-θ
  vector<double> e(lo + 1);
  e[0] = sqr(y[0].r) / n;
  e[lo] = sqr(y[0].i) / n;
  for (int64_t i = 1; i < lo; i++)
    e[i] = 2 * (sqr(y[i].r) + sqr(y[i].i)) / n;
  return e;
}

Line fit_line(span<const double> x, span<const double> y) {
  slow_assert(x.size() == y.size() && x.size() >= 2);
  const double n = x.size();
  double sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (size_t i = 0; i < x.size(); i++) {
    sx += x[i]; sy += y[i]; sxx += x[i]*x[i]; sxy += x[i]*y[i];
  }
  const double b = (n*sxy - sx*sy) / (n*sxx - sx*sx);
  return Line{(sy - b*sx) / n, b};
}

}  // namespace mandelbrot
