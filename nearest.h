// Accurate constants via arb
#pragma once

#include "acb_cc.h"
#include "arb_cc.h"
#include "arf_cc.h"
#include "complex.h"
#include "debug.h"
#include <optional>
namespace mandelbrot {

using std::move;
using std::optional;
using std::remove_const_t;
using std::span;
template<int n> struct Expansion;

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
template<class S,class F> auto nearest(F&& f) {
  const int max_prec = 1600;
  for (int prec = 200; prec <= max_prec; prec <<= 1) {
    const auto c = round_nearest<S>(f(prec), prec);
    if (c) return *c;
  }
  die("nearest ran out of precision (max prec = %g)", max_prec);
}

// 𝜋
template<class S> S nearest_pi();

// sqrt(a/b)
template<class S> S nearest_sqrt(const int64_t a, const int64_t b);

// exp(2𝜋i a/b)
template<class S> Complex<S> nearest_twiddle(const int64_t a, const int64_t b);

// exp(2𝜋i a/b) for a ∈ [0,zs.size())
template<class S> void nearest_twiddles(span<Complex<S>> zs, const int64_t b);

}  // namespace mandelbrot
