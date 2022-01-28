// Accurate constants via arb
#pragma once

#include "arf_cc.h"
#include "complex.h"
#include "debug.h"
#include <optional>
namespace mandelbrot {

using std::move;
using std::optional;
using std::remove_const_t;
using std::span;

template<class S> S round_near(const arf_t c);
template<> inline double round_near(const arf_t c) { return arf_get_d(c, ARF_RND_NEAR); }

// 𝜋
template<class S> S nearest_pi();

// sqrt(a/b)
template<class S> S nearest_sqrt(const int64_t a, const int64_t b);

// exp(2𝜋i a/b)
template<class S> Complex<S> nearest_twiddle(const int64_t a, const int64_t b);

// exp(2𝜋i a/b) for a ∈ [0,zs.size())
template<class S> void nearest_twiddles(span<Complex<S>> zs, const int64_t b);

}  // namespace mandelbrot
