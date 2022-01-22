// Mandelbrot area via Arb

#include "arb_area.h"
#include "area.h"
#include "debug.h"
#include <exception>

using namespace mandelbrot;

int main() {
  try {
    const int max_k = 14;
    if (0) {
      const int prec = 200;
      arb_areas(max_k, prec);
    } else
      areas<double>(max_k);
    return 0;
  } catch (const std::exception& e) {
    die(e.what());
  }
}
