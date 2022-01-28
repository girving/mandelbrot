// C++ interface to arb_t

#include "arb_cc.h"
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

string Arb::safe(const slong n) const {
  char* p = arb_get_str(x, n, 0);
  string s(p);
  free(p);
  return s;
}

}  // namespace mandelbrot
