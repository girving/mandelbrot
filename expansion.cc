// Higher precision floating point numbers via expansion arithmetic

#include "expansion.h"
#include "arb_cc.h"
#include "arf_cc.h"
#include "debug.h"
#include "format.h"
#include "nearest.h"
#include "span.h"
namespace mandelbrot {

template<int n> int sign(const Expansion<n> x) {
  for (int i = 0; i < n; i++)
    if (x.x[i])
      return x.x[i] > 0 ? 1 : -1;
  return 0;
}

template<int n> Expansion<n> abs(const Expansion<n> x) {
  Expansion<n> ax;
  const double s = sign(x);
  for (int i = 0; i < n; i++)
    ax.x[i] = s * x.x[i];
  return ax;
}

template<int n> Expansion<n>::operator bool() const {
  for (int i = 0; i < n; i++)
    if (x[i])
      return true;
  return false;
}

template<int n> bool Expansion<n>::operator==(const Expansion y) const {
  return !bool(*this - y);
}

template<int n> Expansion<n> inv(const Expansion<n> x) {
  // Base case
  Expansion<n> y;
  y.x[0] = 1 / x.x[0];

  // Newton step:
  //   1/y = x
  //   f(y) = 1/y - x
  //   f'(y) = -1/y^2
  //   N(y) = y - f(y) / f'(y)
  //        = y - (1/y - x) / (-1/y^2)
  //        = y - y(xy - 1)
  //        = y(2 - xy)
  const Expansion<n> two(2);
  for (int i = 0; i < n; i++)
    y = y*(two - x*y);
  return y;
}

template<int n> Expansion<n> Expansion<n>::operator/(const Expansion b) const {
  // Compute the inverse and multiply
  const auto inv_b = inv(b);
  Expansion y = *this * inv_b;

  // One more step of Newton refinement:
  //   y = a/b
  //   f(y) = by - a
  //   f'(y) = b
  //   N(y) = y - (b*y - a)/b
  y = y - (b*y - *this) * inv_b;
  return y;
}

template<int n> Arb Expansion<n>::arb(const int prec) const {
  Arb a, t;
  for (int i = n-1; i >= 0; i--) {
    arb_set_d(t, x[i]);
    arb_add(a, a, t, prec);
  }
  return a;
}

template<int n> Arf exact_arf(const Expansion<n> x) {
  Arf a, t;
  for (int i = 0; i < n; i++) {
    arf_set_d(t, x.x[i]);
    arf_add(a, a, t, ARF_PREC_EXACT, ARF_RND_NEAR);
  }
  return a;
}

template<int n> Arb exact_arb(const Expansion<n> x) {
  Arb b;
  arb_set_arf(b, exact_arf(x));
  return b;
}

template<int n> span<const double> Expansion<n>::span() const {
  return SPAN_NAMESPACE::span<const double>(x, size_t(n));
}

template<int n> ostream& operator<<(ostream& out, const Expansion<n> e) {
  return out << exact_arf(e);
}

template<int n> string safe(const Expansion<n> x) {
  const auto a = exact_arf(x);
  for (int p = 1; p < 20*n; p++) {
    const string s = format("%.*g", p, a);
    const Expansion<n> y(s);
    if (x == y)
      return s;
  }
  die("ran out of precision trying to stringify x = %.100g (%.17g)", x, x.span());
}

template<int n> Expansion<n>::Expansion(const string& s) {
  *this = nearest<Expansion<n>>([&s](const int prec) {
    Arb a;
    arb_set_str(a, s.c_str(), prec);
    return a;
  });
}

#define N(n) \
  template int sign(const Expansion<n>); \
  template Expansion<n> abs(const Expansion<n>); \
  template Expansion<n>::operator bool() const; \
  template bool Expansion<n>::operator==(const Expansion) const; \
  template Arb Expansion<n>::arb(const int) const; \
  template Arb exact_arb(const Expansion<n>); \
  template Arf exact_arf(const Expansion<n>); \
  template span<const double> Expansion<n>::span() const; \
  template ostream& operator<<(ostream&, const Expansion<n>); \
  template Expansion<n> inv(const Expansion<n>); \
  template Expansion<n> Expansion<n>::operator/(const Expansion<n>) const; \
  template Expansion<n>::Expansion(const string&); \
  template string safe(const Expansion<n>);
N(2)
N(3)

}  // namespace mandelbrot