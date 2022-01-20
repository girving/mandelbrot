// C++ interface to arb_t

#include "arb-cc.h"
#include <iostream>
namespace mandelbrot {

std::ostream& operator<<(std::ostream& out, const Arb& a) {
  char* buffer;
  size_t size;
  FILE* f = open_memstream(&buffer, &size);
  arb_fprintd(f, a.x, out.precision());
  fflush(f);
  out << buffer;
  fclose(f);
  return out;
}

}  // namespace mandelbrot
