// C++ interface to acb_t

#include "acb_cc.h"
#include <iostream>
namespace mandelbrot {

std::ostream& operator<<(std::ostream& out, const Acb& a) {
  char* buffer;
  size_t size;
  FILE* f = open_memstream(&buffer, &size);
  acb_fprintd(f, a.x, out.precision());
  fflush(f);
  out << buffer;
  fclose(f);
  return out;
}

}  // namespace mandelbrot
