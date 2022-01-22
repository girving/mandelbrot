// Accurate constants via arb
#pragma once

#include "arb_cc.h"
#include "arf_cc.h"
#include "complex.h"
#include "debug.h"
namespace mandelbrot {

template<class S> S round_near(const arf_t c);
template<> inline double round_near(const arf_t c) { return arf_get_d(c, ARF_RND_NEAR); }

// Given f : prec → S, increase prec until we get perfect rounding
template<class S,class F> S nearest(F&& f) {
  const int max_prec = 1600;
  Arb c;
  Arf a, b;
  for (int prec = 200; prec <= max_prec; prec <<= 1) {
    f(c, prec);
    arb_get_interval_arf(a, b, c, prec);
    const auto a_ = round_near<S>(a);
    const auto b_ = round_near<S>(b);
    if (a_ == b_)
      return a_;
  } 
  die("nearest ran out of precision: c = %g, max prec = %g", c, max_prec);
}

// 𝜋
template<class S> S nearest_pi() {
  static const S pi = nearest<S>([](Arb& c, int prec) { arb_const_pi(c, prec); });
  return pi;
}

// exp(2𝜋i a/b)
template<class S> Complex<S> nearest_cis_tau(const int64_t a, const int64_t b) {
  fmpq_t t;
  fmpq_init(t);
  fmpq_set_si(t, 2*a, b);
  const auto c = nearest([t](Arb& c, int prec) { arb_cos_pi_fmpq(c, t, prec); });
  const auto s = nearest([t](Arb& c, int prec) { arb_sin_pi_fmpq(c, t, prec); });
  fmpq_clear(t);
  return Complex<S>(c, s);
}

}  // namespace mandelbrot
