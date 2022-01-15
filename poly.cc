// C++ interface with arb_poly_t

#include "poly.h"
#include "tinyformat.h"
namespace mandelbrot {

using std::max;
using std::numeric_limits;
using std::ostream;
using tinyformat::format;

arb_srcptr Poly::operator[](const slong n) const {
  static Arb zero;
  const auto p = arb_poly_get_coeff_ptr(x, n);
  return p ? p : zero.x;
}

string Poly::stats() const {
  const int prec = 200;
  Arb lo, hi;
  arb_pos_inf(lo);
  arb_neg_inf(hi);
  double radius = numeric_limits<double>::infinity();
  const auto n = length();
  for (int i = 0; i < n; i++) {
    const auto c = (*this)[i];
    arb_min(lo, lo, c, prec);
    arb_max(hi, hi, c, prec);
    radius = max(radius, mag_get_d(arb_radref(c)));
  }
  return format("degree %d, coeffs [%.3g,%.3g], radius %g", degree(), lo, hi, radius);
}

void safe_poly_mullow(arb_poly_t fg, const arb_poly_t f, const arb_poly_t g, const int n, const int prec) {
  if (n > 0)
    arb_poly_mullow(fg, f, g, n, prec);
  else
    arb_poly_zero(fg);
}

void Poly::assert_low_zero(const slong n) const {
  for (int i = 0; i < n; i++)
    if (!arb_contains_si((*this)[i], 0))
      throw std::runtime_error("low terms are not zero");
}

ostream& operator<<(ostream& out, const Poly& f) {
  char* buffer;
  size_t size; 
  FILE* file = open_memstream(&buffer, &size);
  arb_poly_fprintd(file, f.x, out.precision());
  fflush(file);
  out << buffer;
  fclose(file);
  return out;
}

}  // namespace mandelbrot
