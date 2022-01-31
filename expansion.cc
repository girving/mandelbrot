// Higher precision floating point numbers via expansion arithmetic

#include "expansion.h"
#include "arb_cc.h"
#include "arf_cc.h"
#include "format.h"
#include <span>
namespace mandelbrot {

using std::span;

template<int n> Expansion<n>::operator bool() const {
  for (int i = 0; i < n; i++)
    if (x[i])
      return true;
  return false;
}

template<int n> bool Expansion<n>::operator==(const Expansion y) const {
  return !bool(*this - y);
}

template<int n> Arb Expansion<n>::arb(const int prec) const {
  Arb a, t;
  for (int i = n-1; i >= 0; i--) {
    arb_set_d(t, x[i]);
    arb_add(a, a, t, prec);
  }
  return a;
}

template<int n> Arf Expansion<n>::exact_arf() const {
  Arf a, t;
  for (int i = 0; i < n; i++) {
    arf_set_d(t, x[i]);
    arf_add(a, a, t, ARF_PREC_EXACT, ARF_RND_NEAR);
  }
  return a;
}

template<int n> std::span<const double> Expansion<n>::span() const {
  return std::span<const double>(x, size_t(n));
}

template<int n> ostream& operator<<(ostream& out, const Expansion<n> e) {
  return out << e.exact_arf();
}

#define N(n) \
  template Expansion<n>::operator bool() const; \
  template bool Expansion<n>::operator==(const Expansion) const; \
  template Arb Expansion<n>::arb(const int) const; \
  template Arf Expansion<n>::exact_arf() const; \
  template std::span<const double> Expansion<n>::span() const; \
  template ostream& operator<<(ostream&, const Expansion<n>);
N(2)
N(3)
N(4)

}  // namespace mandelbrot
