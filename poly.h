// C++ interface to arb_poly_t
#pragma once

#include <arb_poly.h>
#include "arb.h"
#include "print.h"
namespace mandelbrot {

using std::swap;

struct Poly {
  arb_poly_t x;

  Poly() { arb_poly_init(x); }
  Poly(const Poly& g) = delete;
  Poly(Poly&& g) { *x = *g.x; arb_poly_init(g.x); }
  explicit Poly(const slong a) { arb_poly_init(x); arb_poly_set_coeff_si(x, 0, a); }
  explicit Poly(const Arb& a) { arb_poly_init(x); arb_poly_set_coeff_arb(x, 0, a.x); }
  ~Poly() { arb_poly_clear(x); }

  // Assignment
  void operator=(const slong a) { arb_poly_zero(x); arb_poly_set_coeff_si(x, 0, a); }
  void operator=(const Arb& a) { arb_poly_zero(x); arb_poly_set_coeff_arb(x, 0, a); }
  void operator=(const Poly& g) { arb_poly_set(x, g.x); }
  void operator=(Poly&& g) { swap(*x, *g.x); arb_poly_zero(g.x); }

  // Implicit conversion
  operator arb_poly_struct*() { return x; }
  operator const arb_poly_struct*() const { return x; }

  // Information
  slong degree() const { return arb_poly_degree(x); }
  slong length() const { return arb_poly_length(x); }
  arb_srcptr operator[](const slong n) const;
  arb_ptr mutable_at(const slong n) { return arb_poly_get_coeff_ptr(x, n); }
  string stats() const;
  friend std::ostream& operator<<(std::ostream& out, const Poly& f);
  void assert_low_zero(const slong n) const;

  // Shifts
  void operator<<=(const slong n) { arb_poly_shift_left(x, x, n); }
  void operator>>=(const slong n) { arb_poly_shift_right(x, x, n); }
};

template<class... Args> void stats(const Poly& f, const char* fmt, const Args&... args) {
  print("    %s: %s", format(fmt, args...), f.stats());
}

void safe_poly_mullow(arb_poly_t fg, const arb_poly_t f, const arb_poly_t g, const int n, const int prec);

}  // namespace mandelbrot
