// C++ interface to arb_t

#include "arb_cc.h"
#include "arf_cc.h"
#include "mag_cc.h"
#include <iostream>
namespace mandelbrot {

using std::numeric_limits;

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

string Arb::safe() const {
  const auto n = numeric_limits<int64_t>::max();
  char* p = arb_get_str(x, n, 0);
  string s(p);
  free(p);
  return s;
}

double bound(const Arb& x) {
  Mag m;
  arb_get_mag(m, x);
  return mag_get_d(m);
}

}  // namespace mandelbrot
