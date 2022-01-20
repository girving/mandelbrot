// C++ interface to arb_poly_t
#pragma once

#include <arb_poly.h>
#include "arb-cc.h"
#include "print.h"
namespace mandelbrot {

using std::max;
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
  bool alias(const Poly& f) const { return +x == +f.x; }

  // Shifts
  void operator<<=(const slong n) { arb_poly_shift_left(x, x, n); }
  void operator>>=(const slong n) { arb_poly_shift_right(x, x, n); }
};

template<class... Args> void stats(const Poly& f, const char* fmt, const Args&... args) {
  print("    %s: %s", format(fmt, args...), f.stats());
}

void safe_poly_mullow(arb_poly_t fg, const arb_poly_t f, const arb_poly_t g, const slong n, const slong prec);
void poly_mid(Poly& mid, const Poly& f);
void poly_add_arb(Poly& f, const Arb& a, const slong prec);
bool overlaps(const Poly& f, const Poly& g);

static inline void safe_poly_set_trunc(Poly& f, const Poly& g, const slong n) {
  arb_poly_set_trunc(f, g, max(slong(0), n));
}

// h = f +/- z^s g
void poly_add_shift_series(Poly& h, const Poly& f, const Poly& g, const slong s, const slong n, const slong prec);
void poly_sub_shift_series(Poly& h, const Poly& f, const Poly& g, const slong s, const slong n, const slong prec);

// C = intersection(C, A - B)
void poly_intersect_sub(Poly& C, const Poly& A, const Poly& B, const slong n, const slong prec);

// Strengthened Newton iteration (arb + one more Newton step)
void poly_inv_refine(Poly& y, const Poly& x, const slong n, const slong prec);
void poly_div_refine(Poly& y, const Poly& a, const Poly& b, const slong n, const slong prec);
void poly_log_refine(Poly& y, const Poly& x, const slong n, const slong prec);
void poly_log1p_refine(Poly& y, const Poly& x, const slong n, const slong prec);
void poly_exp_refine(Poly& y, const Poly& x, const slong n, const slong prec);

// Shifted versions of log1p, expm1, log1p_exp
//   log1p_shift: y = z^-s log(1 + z^s x) + O(z^n)
//   expm1_shift: y = z^-s (exp(az^s x) - 1) + O(z^n) with |a| = 1
//   log1p_exp_shift: y = z^-s log (1 + z^s e^x) + O(z^n)
void poly_log1p_shift_refine(Poly& y, const Poly& x, const slong s, const slong n, const slong prec);
void poly_expm1_shift_refine(Poly& y, const Poly& x, const slong a, const slong s, const slong n, const slong prec);
void poly_log1p_exp_shift_refine(Poly& y, const Poly& x, const slong s, const slong n, const slong prec);

}  // namespace mandelbrot
