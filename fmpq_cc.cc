// C++ interface to fmpq_t

#include "fmpq_cc.h"
#include <iostream>
namespace mandelbrot {

std::ostream& operator<<(std::ostream& out, const Fmpq& a) {
  char* buffer;
  size_t size;
  FILE* f = open_memstream(&buffer, &size);
  fmpq_fprint(f, a.x);
  fflush(f);
  out << buffer;
  fclose(f);
  return out;
}

}  // namespace mandelbrot
