// C++ interface to arf_t

#include "arf_cc.h"
#include <iostream>
namespace mandelbrot {

string Arf::safe() const {
  char* buffer;
  size_t size;
  FILE* f = open_memstream(&buffer, &size);
  arf_fprint(f, x);
  fflush(f);
  string s(buffer);
  fclose(f);
  return s;
}

std::ostream& operator<<(std::ostream& out, const Arf& a) {
  char* buffer;
  size_t size;
  FILE* f = open_memstream(&buffer, &size);
  arf_fprintd(f, a.x, out.precision());
  fflush(f);
  out << buffer;
  fclose(f);
  return out;
}

}  // namespace mandelbrot
