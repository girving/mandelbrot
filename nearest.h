// Accurate constants via arb
#pragma once

#include "acb_cc.h"
#include "arb_cc.h"
#include "arf_cc.h"
#include "complex.h"
#include "debug.h"
#include <optional>
namespace mandelbrot {

using std::is_same_v;
using std::move;
using std::optional;
using std::remove_const_t;
template<int n> struct Expansion;

template<class S> static inline const char* nice_type() {
  if constexpr (is_same_v<S,double>) return "double";
  else if constexpr (is_same_v<S,Expansion<2>>) return "Expansion<2>";
  else { static_assert(is_same_v<S,Expansion<3>>); return "Expansion<3>"; }
}

template<class S> struct RoundNear;
template<> struct RoundNear<double> {
  static double round(const arf_t c) { return arf_get_d(c, ARF_RND_NEAR); }
};
template<int n> struct RoundNear<Expansion<n>> {
  static Expansion<n> round(const arf_t c);
};
template<class S> S round_near(const arf_t c) { return RoundNear<S>::round(c); }

// Round both lower and upper bounds to the nearest value, and return nothing if they don't match
template<class S> optional<S> round_nearest(const arb_t c, const int prec);
template<class S> optional<Complex<S>> round_nearest(const acb_t z, const int prec);

// Given f : prec → S, increase prec until we get perfect rounding
template<class S,class F,class C> auto nearest(F&& f, C&& context) {
  const int max_prec = 3200;
  for (int prec = 200;; prec *= 2) {
    const auto fp = f(prec);
    const auto c = round_nearest<S>(fp, prec);
    if (c) return *c;
    if (2*prec > max_prec)
      die("%s: ran out of precision (%s, max prec = %g)\n  f(%d) = %s",
          context(), nice_type<S>(), max_prec, prec, fp.safe());
  }
}

// 𝜋
template<class S> S nearest_pi();

// sqrt(a/b)
template<class S> S nearest_sqrt(const int64_t a, const int64_t b);

// exp(2𝜋i a/b)
template<class S> Complex<S> nearest_twiddle(const int64_t a, const int64_t b);

// exp(2𝜋i a/b) for a ∈ [0,zs.size())
// We first compute factored using fast_prec, filling in holes using nonfactored evaluation.
// Returns the number of times the fallback occurs.
template<class S> int64_t nearest_twiddles(span<Complex<S>> zs, const int64_t b, const int fast_prec);

}  // namespace mandelbrot
