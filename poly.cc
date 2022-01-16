// C++ interface with arb_poly_t

#include "poly.h"
#include "arf.h"
#include "tinyformat.h"
namespace mandelbrot {

using std::max;
using std::min;
using std::numeric_limits;
using std::ostream;
using std::runtime_error;
using tinyformat::format;

arb_srcptr Poly::operator[](const slong n) const {
  static Arb zero;
  const auto p = arb_poly_get_coeff_ptr(x, n);
  return p ? p : zero.x;
}

string Poly::stats() const {
  const int prec = 200;
  Arf lo, hi, t;
  arf_pos_inf(lo);
  arf_neg_inf(hi);
  double radius = -numeric_limits<double>::infinity();
  const auto n = length();
  for (int i = 0; i < n; i++) {
    const auto c = (*this)[i];
    arb_get_lbound_arf(t, c, prec);
    arf_min(lo, lo, t);
    arb_get_ubound_arf(t, c, prec);
    arf_max(hi, hi, t);
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

void poly_mid(Poly& mid, const Poly& f) {
  const slong n = f.length();
  arb_poly_fit_length(mid, n);
  for (int i = 0; i < n; i++)
    arb_get_mid_arb(mid.x->coeffs + i, f[i]);
  _arb_poly_set_length(mid, n);
}

void poly_intersect_sub(Poly& C, const Poly& A, const Poly& B, slong n, const slong prec) {
  n = min(n, C.length());
  n = min(n, max(A.length(), B.length()));
  _arb_poly_set_length(C, n);
  Arb t;
  for (int i = 0; i < n; i++) {
    arb_sub(t, A[i], B[i], prec);
    auto c = C.mutable_at(i);
    if (!arb_overlaps(c, t))
      throw runtime_error(format("poly_intersect_sub: i %d, c %.10g, t %.10g", i, Arb(c), t));
    arb_intersection(c, c, t, prec);
  }
}

// Skip all refinement (for ablation purposes)
static const bool skip_refine = false;

void poly_inv_refine(Poly& y, const Poly& x, slong n, const slong prec) {
  if (!arb_contains_si(x[0], 1))
    throw runtime_error("poly_inv_refine: expected leading coefficient 1");
  // 1/y = x
  // f(y) = 1/y - x
  // f'(y) = -1/y^2
  // N(y) = y0 - f(y0) / f'(y)
  //      = y0 - (1/y0 - x) / (-1/y^2)
  //      ⊂ y0 + (y - xy^2)
  //      = y0 + (y - xy^2)
  //      = y0 - y(xy - 1)
  // N(y) ⊂ intersect_{y0 ∈ y} (y0 - y(xy-1))
  arb_poly_inv_series(y, x, n, prec);
  if (skip_refine) return;
  n = min(n, y.length());
  Poly t;
  arb_poly_mullow(t, x, y, n, prec);  // t = xy
  arb_poly_add_si(t, t, -1, prec);  // t = xy - 1
  arb_poly_mullow(t, t, y, n, prec);  // t = y(xy - 1)

  // y = intersect_{y0 ∈ y} (y0 - t)
  //   = (ylo - t) ∩ (yhi - t)
  Arb mid, rad, lo, hi;
  for (int i = 0; i < n; i++) {
    const auto ti = t[i];
    if (!arb_contains_si(ti, 0))
      throw runtime_error("poly_inv_refine: ti does not contain zero");
    auto yi = y.mutable_at(i);
    arb_get_mid_arb(mid, yi);
    arb_set_interval_mag(rad, arb_radref(yi), arb_radref(yi), 2*prec);
    arb_sub(lo, mid, rad, 2*prec);
    arb_add(hi, mid, rad, 2*prec);
    arb_sub(lo, lo, ti, 2*prec);
    arb_sub(hi, hi, ti, 2*prec);
    if (!arb_overlaps(lo, hi))
      throw runtime_error("poly_inv_refine: lo and hi do not overlap");
    arb_intersection(lo, lo, hi, 2*prec);
    if (!arb_overlaps(yi, lo))
      throw runtime_error("poly_inv_refine: yi and lo do not overlap");
    arb_intersection(yi, yi, lo, prec);
  }
}

void poly_div_refine(Poly& y, const Poly& a, const Poly& b, const slong n, const slong prec) {
  if (skip_refine) {
    arb_poly_div_series(y, a, b, n, prec);
    return;
  } else if (+y.x == +a.x || +y.x == +b.x) {
    Poly t;
    poly_div_refine(t, a, b, n, prec);
    y = t;
    return;
  }
  // y = a/b
  // f(y) = by - a
  // f'(y) = b
  // N(y) = y0 - (b*y0 - a)/b
  //      = y0 - (b*y0 - a)(1/b)
  Poly inv_b;
  poly_inv_refine(inv_b, b, n, prec);  // inv_b = 1/b
  arb_poly_mullow(y, a, inv_b, n, prec);  // y = a(1/b) = a/b
  Poly y0;
  poly_mid(y0, y);
  Poly dy;
  arb_poly_mullow(dy, b, y0, n, prec);  // dy = b*y0
  arb_poly_sub_series(dy, dy, a, n, prec);  // dy = b*y0 - a
  arb_poly_mullow(dy, dy, inv_b, n, prec);  // dy = (b*y0 - a)(1/b)
  poly_intersect_sub(y, y0, dy, n, prec);
}

void poly_log_refine(Poly& y, const Poly& x, const slong n, const slong prec) {
  if (skip_refine) {
    arb_poly_log_series(y, x, n, prec);
    return;
  } else if (!arb_contains_si(x[0], 1))
    throw runtime_error("poly_log_refine: expected 1 constant term");
  // y' = x'/x
  Poly dx;
  arb_poly_derivative(dx, x, prec);
  poly_div_refine(y, dx, x, n, prec);
  arb_poly_integral(y, y, prec);
}

void poly_log1p_refine(Poly& y, const Poly& x, const slong n, const slong prec) {
  if (skip_refine) {
    arb_poly_log1p_series(y, x, n, prec);
    return;
  } else if (!arb_contains_si(x[0], 0))
    throw runtime_error("poly_log1p_refine: expected 0 constant term");
  Poly x1;
  arb_poly_set_trunc(x1, x, n);
  arb_poly_set_coeff_si(x1, 0, 1);
  poly_log_refine(y, x1, n, prec);
}

void poly_exp_refine(Poly& y, const Poly& x, const slong n, const slong prec) {
  // y = exp(x)
  // f(y) = log(y) - x
  // f'(y) = 1/y
  // N(y) = y0 - f(y0)/f'(y)
  //      = y0 - (log(y0) - x)/(1/y)
  //      = y0 - y*(log(y0) - x)
  arb_poly_exp_series(y, x, n, prec);
  if (skip_refine) return;
  Poly y0;
  poly_mid(y0, y);
  Poly dy;
  poly_log_refine(dy, y0, n, prec);  // dy = log(y0)
  arb_poly_sub_series(dy, dy, x, n, prec);  // dy = log(y0) - x
  arb_poly_mullow(dy, y, dy, n, prec);  // dy = y(log(y0) - x)
  poly_intersect_sub(y, y0, dy, n, prec);
}

}  // namespace mandelbrot
