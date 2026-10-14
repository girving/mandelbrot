// Escape-time classification of Mandelbrot parameters

#include "escape.h"
namespace mandelbrot {

EscapeDE escape_de(const double x, const double y, const int64_t max_iter) {
  OrbitDE o;
  if (!o.start(x, y))
    o.run(max_iter, max_iter);
  return o.result();
}

Escape escape(const double x, const double y, const int64_t max_iter) {
  Orbit<double> o;
  if (!o.start(x, y))
    o.finish(max_iter);
  return o.result();
}

}  // namespace mandelbrot
