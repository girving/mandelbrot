// C++ interface to acb_t

#include "acb_cc.h"
#include "format.h"
#include <iostream>
namespace mandelbrot {

using std::numeric_limits;

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

string Acb::safe() const {
  const auto n = numeric_limits<int64_t>::max();
  return format("%.*g", n, *this);
}

}  // namespace mandelbrot
