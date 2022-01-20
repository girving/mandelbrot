// Mandelbrot area via Arb

#include "area.h"
#include "debug.h"
#include <exception>

int main() {
  try {
    const int max_k = 14;
    const int prec = 200;
    mandelbrot::areas(max_k, prec);
    return 0;
  } catch (const std::exception& e) {
    mandelbrot::die(e.what());
  }
}
