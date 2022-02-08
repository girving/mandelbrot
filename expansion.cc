// Higher precision floating point numbers via expansion arithmetic

#include "expansion.h"
#include "arb_cc.h"
#include "arf_cc.h"
#include "debug.h"
#include "format.h"
#include "nearest.h"
#include "print.h"
#include "span.h"
#include <random>
namespace mandelbrot {

using std::bernoulli_distribution;
using std::uniform_int_distribution;
using std::uniform_real_distribution;

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
  return std::span<const double>(x, size_t(n));
}

template<int n> ostream& operator<<(ostream& out, const Expansion<n> e) {
  return out << exact_arf(e);
}

template<int n> string maybe_nice_safe(const Expansion<n> x) {
  const auto a = exact_arf(x);
  for (int p = 15*n-2; p <= 30*n; p++) {
    auto s = format("%.*g", p, a);
    if (x == Expansion<n>(s)) return s;
  }
  return string();
}

template<int n> string safe(const Expansion<n> x) {
  auto s = maybe_nice_safe(x);
  return s.size() ? s : format("%.17g", x.span());
}

template<int n> Expansion<n>::Expansion(string_view s) : Expansion(string(s)) {}

template<int n> Expansion<n>::Expansion(const string& s) {
  slow_assert(s.size(), s);
  if (s[0] == '[') {
    const char* p = s.c_str() + 1;
    for (int i = 0; i < n; i++) {
      char* end;
      x[i] = strtod(p, &end);
      slow_assert(end > p && *end == (i+1 < n ? ',' : ']'), s);
      p = end + 1;
    }
    slow_assert(p == s.c_str() + s.size(), s);
  } else {
    *this = nearest<Expansion<n>>([&s](const int prec) {
      Arb a;
      arb_set_str(a, s.c_str(), prec);
      return a;
    });
  }
}

template<int n> Expansion<n> random_expansion_with_exponent(mt19937& mt, int e) {
  Expansion<n> a;
  for (int i = 0; i < n; i++) {
    a.x[i] = ldexp(uniform_real_distribution<double>(-1, 1)(mt), e);
    e = exponent(a.x[i]) - 52 - 1;
  }
  return a;
}

template<int n> Expansion<n> random_expansion(mt19937& mt) {
  const int e = uniform_int_distribution<int>(-20, 20)(mt);
  return random_expansion_with_exponent<n>(mt, e);
}

template<int n> Expansion<n> random_expansion_near(mt19937& mt, const Expansion<n> x) {
  const int e = exponent(x.x[0]) - uniform_int_distribution<int>(1, 52*n)(mt);
  return (bernoulli_distribution(0.5)(mt) ? x : -x) + random_expansion_with_exponent<n>(mt, e);
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
  template Expansion<n>::Expansion(const string&); \
  template Expansion<n>::Expansion(string_view); \
  template string safe(const Expansion<n>); \
  template string maybe_nice_safe(const Expansion<n>); \
  template Expansion<n> random_expansion(mt19937&); \
  template Expansion<n> random_expansion_near(mt19937&, const Expansion<n>); \
  template Expansion<n> random_expansion_with_exponent(mt19937&, int);
N(2)
N(3)

}  // namespace mandelbrot
