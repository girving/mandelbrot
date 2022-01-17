// C++ interface with arb_poly_t

#include "poly.h"
#include "arf.h"
#include "tinyformat.h"
#include <vector>
namespace mandelbrot {

using std::max;
using std::min;
using std::numeric_limits;
using std::ostream;
using std::runtime_error;
using std::vector;
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

bool overlaps(const Poly& f, const Poly& g) {
  const int n = max(f.length(), g.length());
  for (int i = 0; i < n; i++)
    if (!arb_overlaps(f[i], g[i]))
      return false;
  return true;
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
  if (n <= 0) {
    arb_poly_truncate(y, 0);
    return;
  }
  if (y.alias(x)) {
    Poly t;
    poly_inv_refine(t, x, n, prec);
    y = t;
    return;
  }
  // 1/y = x
  // f(y) = 1/y - x
  // f'(y) = -1/y^2
  // N(y) = y0 - f(y0) / f'(y)
  //      = y0 - (1/y0 - x) / (-1/y^2)
  //      = y0 - y0(x y0 - 1)(y/y0)^2

  arb_poly_inv_series(y, x, n, prec);
  if (skip_refine) return;
  n = min(n, y.length());
  Poly y0, t, u, dy;
  poly_mid(y0, y);
  arb_poly_mullow(u, x, y0, n, prec);      // u = xy0
  arb_poly_add_si(u, u, -1, prec);         // u = xy0 - 1
  arb_poly_mullow(t, u, y0, n, prec);      // t = y0(xy0 - 1)
  arb_poly_div_series(u, y, y0, n, prec);  // u = y/y0
  arb_poly_mullow(u, u, u, n, prec);       // u = (y/y0)^2
  arb_poly_mullow(dy, t, u, n, prec);      // dy = y0(xy0 - 1)(y/y0)^2
  poly_intersect_sub(y, y0, dy, n, prec);
}

void poly_div_refine(Poly& y, const Poly& a, const Poly& b, const slong n, const slong prec) {
  if (n <= 0) {
    arb_poly_truncate(y, 0);
    return;
  } else if (skip_refine) {
    arb_poly_div_series(y, a, b, n, prec);
    return;
  } else if (y.alias(a) || y.alias(b)) {
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

void poly_add_arb(Poly& f, const Arb& a, const slong prec) {
  if (f.length()) {
    auto f0 = f.mutable_at(0);
    arb_add(f0, f0, a, prec);
  } else
    arb_poly_set_coeff_arb(f, 0, a);
}

// h = f + z^s g
void poly_add_shift_series(Poly& h, const Poly& f, const Poly& g, const slong s, const int n, const int prec) {
  if (n <= 0)
    return;
  if (n <= s) {
    safe_poly_set_trunc(h, f, n);
    return;
  }
  arb_poly_fit_length(h, n);
  for (int i = n-1; i >= s; i--)
    arb_add(h.x->coeffs + i, f[i], g[i - s], prec);
  if (!h.alias(f))
    for (int i = s-1; i >= 0; i--)
      arb_set(h.x->coeffs + i, f[i]);
  _arb_poly_set_length(h, n);
}

// h = f - z^s g
void poly_sub_shift_series(Poly& h, const Poly& f, const Poly& g, const slong s, const int n, const int prec) {
  if (n <= 0)
    return;
  if (n <= s) {
    safe_poly_set_trunc(h, f, n);
    return;
  }
  arb_poly_fit_length(h, n);
  for (int i = n-1; i >= s; i--)
    arb_sub(h.x->coeffs + i, f[i], g[i - s], prec);
  if (!h.alias(f))
    for (int i = s-1; i >= 0; i--)
      arb_set(h.x->coeffs + i, f[i]);
  _arb_poly_set_length(h, n);
}

void poly_log_refine(Poly& y, const Poly& x, const slong n, const slong prec) {
  if (n <= 0) {
    arb_poly_truncate(y, 0);
    return;
  } else if (skip_refine) {
    arb_poly_log_series(y, x, n, prec);
    return;
  }

  // Grab constant term
  Arb c;
  arb_set(c, x[0]);
  if (!arb_is_positive(c))
    throw runtime_error(format("poly_log_refine: possibly negative constant term %.3g", c));
  const bool one = arb_equal_si(c, 1);

  // log via y' = x'/x
  Poly dx;
  arb_poly_derivative(dx, x, prec);
  poly_div_refine(y, dx, x, n-1, prec);
  arb_poly_integral(y, y, prec);

  // Handle constant term
  if (!one) {
    arb_log(c, c, prec);
    poly_add_arb(y, c, prec);
  }
}

void poly_log1p_refine(Poly& y, const Poly& x, const slong n, const slong prec) {
  if (n <= 0) {
    arb_poly_truncate(y, 0);
    return;
  } else if (skip_refine) {
    arb_poly_log1p_series(y, x, n, prec);
    return;
  } else if (!arb_contains_si(x[0], 0))
    throw runtime_error("poly_log1p_refine: expected 0 constant term");
  else if (y.alias(x)) {
    Poly t;
    poly_log1p_refine(t, x, n, prec);
    y = t;
    return;
  }
  arb_poly_add_si(y, x, 1, prec);
  poly_log_refine(y, y, n, prec);
}

void poly_exp_refine(Poly& y, const Poly& x, const slong n, const slong prec) {
  if (n <= 0) {
    arb_poly_truncate(y, 0);
    return;
  } else if (y.alias(x)) {
    Poly t;
    poly_exp_refine(t, x, n, prec);
    y = t;
    return;
  }
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

// y = z^-s log(1 + z^s x) + O(z^n)
void poly_log1p_shift_refine(Poly& y, const Poly& x, const slong s, const slong n, const slong prec) {
  if (n <= 0) {
    arb_poly_truncate(y, 0);
    return;
  } else if (s < 1)
    throw runtime_error("poly_log1p_shift_refine: need s >= 1");
  // t = 1 + z^s x
  Poly t;
  arb_poly_shift_left(t, x, s);
  arb_poly_truncate(t, n);
  arb_poly_add_si(t, t, 1, prec);
  // y = D_n x
  arb_poly_fit_length(y, n);
  for (int i = 0; i < n; i++)
    arb_mul_si(y.x->coeffs+i, x[i], s+i, prec);
  _arb_poly_set_length(y, n);
  // y = y / t
  poly_div_refine(y, y, t, n, prec);
  // y = J_n y
  const int m = y.length();
  for (int i = 0; i < m; i++) {
    auto yi = y.mutable_at(i);
    arb_div_si(yi, yi, s+i, prec);
  }
}

// Newton iteration, without refinement
template<class Step> static inline void newton_iterate(const slong n0, const slong n, Step&& step) {
  // Determine sizes
  vector<slong> ms;
  for (slong m = n; m > n0; m = (m+1)/2)
    ms.push_back(m);

  // Newton iterate
  int m0 = n0;
  for (int i = ms.size()-1; i >= 0; i--) {
    const slong m = ms[i];
    step(m0, m);
    m0 = m;
  }
}

// y = z^-s (exp(z^s x) - 1) + O(z^n)
void poly_expm1_shift_refine(Poly& y, const Poly& x, const slong a, const slong s, const slong n, slong prec) {
  if (n <= 0) {
    arb_poly_truncate(y, 0);
    return;
  } if (y.alias(x)) {
    Poly t;
    poly_expm1_shift_refine(t, x, a, s, n, prec);
    y = t;
    return;
  } else if (abs(a) != 1)
    throw runtime_error(format("poly_expm1_shift_refine: expected |a| = 1, got a = %d", a));
  else if (s < 1)
    throw runtime_error(format("poly_expm1_shift_refine: need s >= 1, got %d", s));
  // y = z^-s (exp(az^s x) - 1)
  // log1p(y, s) = z^-s log(1 + z^s (z^-s (exp(az^s x) - 1)))
  //             = z^-s log(1 + exp(az^s x) - 1)
  //             = z^-s az^s x
  //             = ax
  // f(y) = log1p(y, s) - ax
  // f'(y) = 1/(1 + z^s y)
  // N(y) = y0 - f(y0) / f'(y)
  //      = y0 - (1 + z^s y)(log1p(y0, s) - ax)

  // Base case: y = ax[0] + O(z)
  arb_poly_set_trunc(y, x, 1);
  if (a < 0)
    arb_poly_neg(y, y);

  Poly t, u;
  const auto delta = [a, s, prec, &x, &t, &u](Poly& dy, const Poly& y0, const Poly& y, const int m0, const int m) {
    // t = log1p(y0) - ax
    poly_log1p_shift_refine(t, y0, s, m, prec);
    (a > 0 ? arb_poly_sub_series : arb_poly_add_series)(t, t, x, m, prec);  // t -= ax
    t.assert_low_zero(m0);
    t >>= m0;
    // u = 1 + z^s y
    arb_poly_shift_left(u, y, s);
    arb_poly_add_si(u, u, 1, prec);
    arb_poly_truncate(u, m);
    // du = z^m0 tu
    safe_poly_mullow(dy, t, u, m-m0, prec);
    dy <<= m0;
  };

  // Approximate all coefficients
  Poly dy;
  newton_iterate(1, n, [&delta, &y, &dy, prec](const int m0, const int m) {
    delta(dy, y, y, m0, m);
    arb_poly_sub_series(y, y, dy, m, prec);
  });

  // One more Newton step
  if (!skip_refine) {
    Poly y0;
    poly_mid(y0, y);
    delta(dy, y0, y, 0, n);
    poly_intersect_sub(y, y0, dy, n, prec);
  }
}

// y = z^-s log (1 + z^s e^x) + O(z^n) = log1p(e^x, s) + O(z^n)
void poly_log1p_exp_shift_refine(Poly& y, const Poly& x, const slong s, const slong n, const slong prec) {
  if (n <= 0) {
    arb_poly_truncate(y, 0);
    return;
  } else if (y.alias(x)) {
    Poly t;
    poly_log1p_exp_shift_refine(t, x, s, n, prec);
    y = t;
    return;
  }
  // y = z^-s log(1 + z^s e^x) + O(z^n)
  // e^(z^s y) = 1 + z^s e^x + O(z^(n+s))
  // f(y) = e^(z^s y) - 1 - z^s e^x
  // f'(y) = z^s e^(z^s y)
  // N(y) = y0 - f(y0) / f'(y)
  //      = y0 - (e^(z^s y0) - 1 - z^s e^x) / (z^s e^(z^s y))
  //      = y0 - (z^-s (e^(z^s y0) - 1) - e^x) / e^(z^s y)
  //      = y0 - (z^-s (e^(z^s y0) - 1) - e^x) / e^(z^s y0) e^(z^s (y0-y))
  //      = y0 - e^(z^s (y0-y)) (z^-s (1 - e^(-z^s y0)) - e^(x-z^s y0))
  //      = y0 + e^(z^s (y0-y)) (expm1(-y0, s) + e^(x-z^s y0))

  // Base case
  if (!n) {
    arb_poly_zero(y);
    return;
  } else {
    Arb c;
    arb_exp(c, x[0], prec);
    arb_poly_set_coeff_arb(y, 0, c);
    arb_poly_truncate(y, 1);
  }

  Poly t, e;
  const auto delta = [s, prec, &x, &t, &e](Poly& dy, const Poly& y, const int m0, const int m) {
    // dy = expm1(-y, s)
    poly_expm1_shift_refine(dy, y, -1, s, m, prec);
    // e = exp(x - z^s y)
    arb_poly_shift_left(t, y, s);
    arb_poly_sub_series(t, x, t, m, prec);
    poly_exp_refine(e, t, m, prec);
    // dy = -expm1(-y, s) - exp(x - z^s y)
    arb_poly_add_series(dy, dy, e, m, prec);
    arb_poly_neg(dy, dy);
    dy.assert_low_zero(m0);
  };

  // Approximate all coefficients
  Poly dy;
  newton_iterate(1, n, [&delta, &y, &dy, prec](const int m0, const int m) {
    delta(dy, y, m0, m);
    arb_poly_sub_series(y, y, dy, m, prec);
  });

  // One more Newton step
  if (!skip_refine) {
    Poly y0;
    poly_mid(y0, y);
    delta(dy, y0, 0, n);

    // dy *= e^(z^s (y-y0))
    arb_poly_zero(t);
    arb_poly_fit_length(t, n);
    for (int i = 0; i < n-s; i++) {
      const auto& r = arb_radref(y[i]);
      arb_set_interval_neg_pos_mag(t.x->coeffs + i + s, r, r, prec);
    }
    _arb_poly_set_length(t, n);
    poly_exp_refine(e, t, n, prec);
    safe_poly_mullow(t, e, dy, n, prec);

    // Incorporate refinement step
    poly_intersect_sub(y, y0, t, n, prec);
  }
}

}  // namespace mandelbrot
