// C++ interface to arf_t

#include "arf-cc.h"
#include "tinyformat.h"
namespace mandelbrot {

using tinyformat::format;

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
