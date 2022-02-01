// Series arithmetic
#pragma once

#include "arith.h"
#include "bit.h"
#include "debug.h"
#include "fft.h"
#include "is_interval.h"
#include "span.h"
#include <algorithm>
#include <cstring>
#include <iostream>
#include <memory>
#include <type_traits>
namespace mandelbrot {

using std::conditional_t;
using std::initializer_list;
using std::is_const_v;
using std::is_trivially_copyable_v;
using std::is_trivially_destructible_v;
using std::ldexp;
using std::max;
using std::min;
using std::move;
using std::numeric_limits;
using std::ostream;
using std::remove_const_t;
using std::shared_ptr;
using std::swap;
template<class T> struct Series;
template<class F> struct SeriesExp;

// For approximate near zero assertions
static inline bool near_zero(const double x) { return abs(x) < 1e-5; }

template<class T> struct Series {
  static_assert(is_trivially_copyable_v<T>);
  static_assert(is_trivially_destructible_v<T>);
  typedef remove_const_t<T> Scalar;
  typedef Scalar value_type;
private:
  typedef Scalar S;
  shared_ptr<T[]> x;
  int64_t known_;  // Known terms
  int64_t nonzero_;  // Possibly nonzero terms
  int64_t limit_;  // Allocated terms

  struct Unusable {};
  template<class A> friend struct Series;
public:

  Series() : known_(0), nonzero_(0), limit_(0) {}
  explicit Series(int64_t limit)
    : known_(0), nonzero_(0), limit_(relu(limit)) {
    if (limit_)
      x.reset(static_cast<S*>(malloc(limit_ * sizeof(S))), [](S* p) { free(p); });
  }
  Series(int64_t limit, initializer_list<S>&& cs)
    : limit_(relu(limit)) {
    slow_assert(cs.size() <= size_t(limit_));
    known_ = nonzero_ = cs.size();
    if (limit_) {
      shared_ptr<S[]> y(static_cast<S*>(malloc(limit_ * sizeof(S))), [](S* p) { free(p); });
      int64_t i = 0;
      for (const auto& c : cs)
        y[i++] = c;
      x = move(y);
    }
  }
  Series(initializer_list<S>&& cs) : Series(cs.size(), move(cs)) {}
  Series(const Series& g) = default;
  Series(const Series<conditional_t<is_const_v<T>,S,Unusable>>& g)
    : x(g.x), known_(g.known()), nonzero_(g.nonzero()), limit_(g.limit()) {}
  Series(Series&& g) : x(move(g.x)), known_(g.known_), nonzero_(g.nonzero_), limit_(g.limit_) {
    g.known_ = g.nonzero_ = g.limit_ = 0;
  }
  ~Series() = default;
  void clear() { x.reset(); known_ = nonzero_ = limit_ = 0; }

  // Assignment
  template<class A> void set_scalar(int64_t known, const A& a) {
    known = relu(known);
    slow_assert(limit_ >= (known > 0));
    known_ = known;
    nonzero_ = known > 0;
    if (known_)
      x[0] = a;
  }
  void operator=(const Series& g) {
    slow_assert(limit_ >= g.nonzero_);
    known_ = g.known_;
    nonzero_ = g.nonzero_;
    memcpy(x.get(), g.x.get(), nonzero_ * sizeof(S));
  }
  void operator=(const Series<conditional_t<is_const_v<T>,Unusable,const S>>& g) {
    slow_assert(limit_ >= g.nonzero_);
    known_ = g.known_;
    nonzero_ = g.nonzero_;
    memcpy(x.get(), g.x.get(), nonzero_ * sizeof(S));
  }
  template<class F> void operator=(SeriesExp<F>&& e) { e.set(*this); }
  void swap(Series<S>& g) {
    using std::swap;
    swap(x, g.x);
    swap(known_, g.known_);
    swap(nonzero_, g.nonzero_);
    swap(limit_, g.limit_);
  }

  // Adding or removing const
  const Series<const S>& const_() const { return *reinterpret_cast<const Series<const S>*>(this); }
  const Series<S>& const_cast_() const { return *reinterpret_cast<const Series<S>*>(this); }

  // Information
  int64_t known() const { return known_; }
  int64_t nonzero() const { return nonzero_; }
  int64_t limit() const { return limit_; }
  bool valid(const int64_t i) const { return (uint64_t)i < (uint64_t)nonzero_; }
  const S& operator[](const int64_t n) const { assert(valid(n)); return x[n]; }
  S& operator[](const int64_t n) { assert(valid(n)); return x[n]; }
  T* data() const { return x.get(); }
  template<class A> bool alias(const Series<A>& f) const {
    return x && !x.owner_before(f.x) && !f.x.owner_before(x);
  }

  // Iteration
  typedef T* iterator;
  typedef T* const_iterator;
  T* begin() const { return x.get(); }
  T* end() const { return x.get() + nonzero_; }
  size_t size() const { return nonzero_; }

  // Verify that low terms vanish
  void assert_low_near_zero(const int64_t n) const {
    slow_assert(n <= known_);
    const auto nz = min(n, nonzero_);
    for (int64_t i = 0; i < nz; i++)
      slow_assert(near_zero(x[i]), "x = %.3g, x[%d] = %.3g != 0", *this, i, x[i]);
  }

  // Output
  friend ostream& operator<<(ostream& out, const Series& f) { return out << f.span(); }

  // Change the number of terms in place
  void set_unknown() { known_ = nonzero_ = 0; }
  void set_counts(const int64_t known, const int64_t nonzero) {
    slow_assert(uint64_t(nonzero) <= uint64_t(min(limit_, known)));
    known_ = known;
    nonzero_ = nonzero;
  }
  void set_known(const int64_t known) {
    slow_assert(known >= 0);
    known_ = known;
    nonzero_ = min(nonzero_, known);
  }

  // Non-aliasing copy
  Series<S> copy(const int64_t limit) const {
    Series<S> g(limit);
    g = *this;
    return g;
  }

  // Span accessors
  SPAN_NAMESPACE::span<T> span() const { return SPAN_NAMESPACE::span<T>(x.get(), nonzero_); }
  SPAN_NAMESPACE::span<T> low_span(int64_t n) const { return SPAN_NAMESPACE::span<T>(x.get(), min(relu(n), nonzero_)); }
  SPAN_NAMESPACE::span<T> high_span(int64_t n) const { return SPAN_NAMESPACE::span<T>(x.get() + n, relu(nonzero_ - n)); }

  // The low terms of a series, without copying
  Series<const S> low(const int64_t n) const {
    const auto nk = min(relu(n), known_);
    const auto nz = min(nk, nonzero_);
    Series<const S> h;
    h.known_ = nk;
    h.nonzero_ = h.limit_ = nz;
    if (nz)
      h.x = x;
    return h;
  }

  // The high terms of a series, without copying
  Series<const S> high(int64_t n) const {
    n = relu(n);
    const auto nk = relu(known_ - n);
    const auto nz = relu(nonzero_ - n);
    Series<const S> h;
    h.known_ = nk;
    h.nonzero_ = h.limit_ = nz;
    if (nz)
      h.x = shared_ptr<const S[]>(x, x.get() + n);
    return h;
  }

  // In-place arithmetic
  void operator+=(const int a) const { slow_assert(nonzero_); x[0] += a; }
  void operator-=(const int a) const { slow_assert(nonzero_); x[0] -= a; }
  void operator+=(const S& a) const { slow_assert(nonzero_); x[0] += a; }
  void operator-=(const S& a) const { slow_assert(nonzero_); x[0] -= a; }

  // self ±= z^s f
  template<int sign> void high_addsub(const int64_t s, const Series<const S>& f) {
    static_assert(!is_const_v<T>);
    static_assert(sign == 1 || sign == -1);
    const auto nk = min(known_, f.known_ + s);
    const auto nz = min(nk, max(nonzero_, f.nonzero_ ? f.nonzero_ + s : 0));
    slow_assert(!alias(f) && nz <= limit_);
    const auto xs = x.get() + s;
    // Fill in newly exposed zeros
    const auto zero = min(s, nz);
    for (int64_t i = nonzero_; i < zero; i++)
      x[i] = 0;
    // Add or subtract in the overlapping region
    const auto both = min(relu(nonzero_ - s), f.nonzero_);
    for (int64_t i = 0; i < both; i++) {
      if constexpr (sign > 0) xs[i] += f.x[i];
      else xs[i] -= f.x[i];
    }
    // Copy remaining values of z^s f
    const auto copy = min(relu(nz - s), f.nonzero_);
    for (int64_t i = both; i < copy; i++) {
      if constexpr (sign > 0) xs[i] = f.x[i];
      else xs[i] = -f.x[i];
    }
    known_ = nk;
    nonzero_ = nz;
  }

  void operator+=(const Series<const S>& f) { high_addsub<1>(0, f); }
  void operator-=(const Series<const S>& f) { high_addsub<-1>(0, f); }
  void high_add(const int64_t s, const Series<const S>& f) { high_addsub<1>(s, f); }
  void high_sub(const int64_t s, const Series<const S>& f) { high_addsub<-1>(s, f); }
};

// Unevaluated series computations
template<class F> struct SeriesExp { F set; };

template<class... Args> struct ScalarT;
template<class T, class... Rest> struct ScalarT<T,Rest...> : public ScalarT<Rest...> {};
template<class T, class... Rest> struct ScalarT<Series<T>,Rest...> { typedef remove_const_t<T> type; };
template<class T, class... Rest> struct ScalarT<const Series<T>&,Rest...> { typedef remove_const_t<T> type; };
template<class... Args> using Scalar = typename ScalarT<Args...>::type;
#define UNPAREN(...) __VA_ARGS__

#define SERIES_EXP(name, y, Ts, xs, args) \
  template<UNPAREN Ts> auto name(UNPAREN args) { \
    auto set = [=](auto& y) { set_##name(y, UNPAREN xs); }; \
    return SeriesExp<decltype(set)>{move(set)}; \
  } \
  template<UNPAREN Ts, class Dst> void set_##name(Dst& y, UNPAREN args)

// Multiplication: z = xy
SERIES_EXP(mul, z, (class A,class B), (x,y), (const Series<A>& x, const Series<B>& y)) {
  typedef remove_const_t<A> S;
  const auto nk = min(x.known(), y.known());
  const auto nz = min(nk, relu(x.nonzero() + y.nonzero() - 1));
  z.set_counts(nk, nz);
  fft_mul<S>(z.span(), x.low_span(nz), y.low_span(nz));
}

// Shifted multiplication: z = x(1 + z^s y)
SERIES_EXP(mul1p, z, (class A,class B), (x,y,s), (const Series<A>& x, const Series<B>& y, const int64_t s)) {
  typedef remove_const_t<A> S;
  slow_assert(!z.alias(x) && s > 0);
  const auto nk = min(x.known(), y.known() + s);
  const auto nz = min(nk, x.nonzero() + y.nonzero() + s - 1);
  const auto both = min(nz, x.nonzero());
  const auto copy = min(s, both);
  const auto zero = min(s, nz);
  z.set_counts(nk, nz);
  fft_mul<S>(z.high_span(s), x.low_span(nz-s), y.low_span(nz-s));
  for (int64_t i = s; i < both; i++)
    z[i] += x[i];
  for (int64_t i = 0; i < copy; i++)
    z[i] = x[i];
  for (int64_t i = copy; i < zero; i++)
    z[i] = 0;
}

// Squaring: y = x^2
SERIES_EXP(sqr, y, (class A), (x), (const Series<A>& x)) {
  typedef remove_const_t<A> S;
  const auto nk = x.known();
  const auto nz = min(nk, relu(2*x.nonzero() - 1));
  y.set_counts(nk, nz);
  fft_sqr<S>(y.span(), x.span());
}

// Number of Newton steps needed to go from n0 to n
static inline int newton_steps(const int64_t n0, const int64_t n) {
  slow_assert(0 < n0 && n0 <= n);
  return int(countr_zero(bit_ceil(uint64_t((n+n0-1)/n0))));
}

// Newton iteration with refinement
template<class Step> static inline void newton_iterate(int64_t n0, const int64_t n, Step&& step) {
  // We want to arrange for the last Newton iteration to be maximally efficient, in the sense of hitting
  // almost exactly the right size.  To do this, we compute these numbers of terms in reverse order
  //   m0 = n
  //   m1 = (n+1)/2
  //   mk = (n+2^k-1)/2^k
  // The first k s.t. mk <= n0 is
  //   mk <= n0
  //   (n+2^k-1)/2^k <= n0
  //   n+2^k-1 <= n0*2^k+2^k-1
  //   n <= n0*2^k
  //   (n+n0-1)/n0 <= 2^k
  //   k = ceil(log2((n+n0-1)/n0))
  // Examples:
  //   n0 = 1, n = 4, k = ceil(log2(4)) = 2
  if (n <= n0) return;
  slow_assert(n0 > 0);
  for (int kr = 2*newton_steps(n0, n)-1; kr >= 0; kr--) {
    // We alternate extension steps with refinement steps
    const int k = kr >> 1;
    const bool refine = !(kr & 1);
    const auto m = (n + (int64_t(1)<<k) - 1) >> k;
    step(n0, m, refine);
    if (refine) n0 = m;
  }
}

// Reciprocal: y = 1 / x
SERIES_EXP(inv, y, (class T), (x), (const Series<T>& x)) {
  typedef remove_const_t<T> S;
  const auto n = x.known();
  if (!n) return y.set_unknown();
  slow_assert(!y.alias(x) && x.nonzero());

  // Base case
  y.set_scalar(1, inv(x[0]));

  // Newton step:
  //   1/y = x
  //   f(y) = 1/y - x
  //   f'(y) = -1/y^2
  //   N(y) = y0 - f(y0) / f'(y)
  //        = y0 - (1/y0 - x) / (-1/y^2)
  //        = y0 - y0(x y0 - 1)(y/y0)^2
  Series<S> dy(n);
  newton_iterate(1, n, [&y, &x, &dy](const int64_t m0, const int64_t m, const bool refine) {
    y.set_known(m);

    // dy = y0(xy0-1)(y/y0)^2
    static_assert(!is_interval<T>);  // Ignore y/y0 for now
    dy = mul(x, y);
    dy -= 1;
    dy = mul(dy, y);

    // Update
    y.high_sub(m0, dy.high(m0));
  });
}

// Shifted reciprocal: y = z^-s (1 / (1 + z^s x) - 1)
SERIES_EXP(inv1p, y, (class T), (x, s), (const Series<T>& x, const int64_t s)) {
  typedef remove_const_t<T> S;
  const auto n = x.known();
  if (!n) return y.set_unknown();
  slow_assert(!y.alias(x) && s > 0);

  // Base case
  y.set_scalar(1, x.nonzero() ? -x[0] : 0);

  // Newton step:
  //   1/(1 + z^s y) = 1 + z^s x
  //   f(y) = 1/(1 + z^s y) - 1 - z^s x
  //   f'(y) = -z^s/(1 + z^s y)^2
  //   N(y) = y0 - f(y0) / f'(y)
  //        = y0 - (1/(1 + z^s y0) - 1 - z^s x) / (-z^s/(1 + z^s y)^2)
  //        = y0 + z^-s (1 + z^s y0)(1 - (1 + z^s y0) - z^s x(1 + z^s y0)) ((1 + z^s y)/(1 + z^s y0))^2
  //        = y0 - (1 + z^s y0)(y0 + x(1 + z^s y0)) ((1 + z^s y)/(1 + z^s y0))^2
  Series<S> dy(n), t(n);
  newton_iterate(1, n, [&y, &x, &dy, &t, s](const int64_t m0, const int64_t m, const bool refine) {
    y.set_known(m);

    // dy = y0(xy0-1)(y/y0)^2
    static_assert(!is_interval<T>);  // Ignore y/y0 for now
    t = mul1p(x, y, s);
    t += y;
    dy = mul1p(t, y, s);

    // Update
    y.high_sub(m0, dy.high(m0));
  });
}

// Division: y = a / b
SERIES_EXP(div, y, (class A,class B), (a,b), (const Series<A>& a, const Series<B>& b)) {
  typedef remove_const_t<A> S;
  slow_assert(!y.alias(a) && !y.alias(b));
  const auto n = min(a.known(), b.known());

  // Compute the inverse and multiply
  Series<S> inv_b(n);
  inv_b = inv(b.low(n));
  y = mul(a, inv_b);

  // One more step of Newton refinement:
  //   y = a/b
  //   f(y) = by - a
  //   f'(y) = b
  //   N(y) = y0 - (b*y0 - a)/b
  //        = y0 - (b*y0 - a)(1/b)
  static_assert(!is_interval<S>);  // Assume y0 = y
  Series<S> dy(n);
  dy = mul(b, y);
  dy -= a;
  dy = mul(dy, inv_b);
  y -= dy;
}

// Shifted division: y = a / (1 + z^s b)
SERIES_EXP(div1p, y, (class A,class B), (a,b,s), (const Series<A>& a, const Series<B>& b, const int64_t s)) {
  typedef remove_const_t<A> S;
  slow_assert(!y.alias(a) && !y.alias(b) && s > 0);
  const auto n = min(a.known(), b.known() + s);

  // Compute the inverse and multiply
  Series<S> inv_b(n - s);
  inv_b = inv1p(b.low(n - s), s);
  y = mul1p(a, inv_b, s);

  // One more step of Newton refinement
  //   c = 1 + z^s b
  //   y = a/c
  //   f(y) = cy - a
  //   f'(y) = c
  //   N(y) = y0 - (c*y0 - a)/c
  //        = y0 - (c*y0 - a)(1/c)
  static_assert(!is_interval<S>);  // Assume y0 = y
  Series<S> t(n), dy(n);
  t = mul1p(y, b, s);
  t -= a;
  dy = mul1p(t, inv_b, s);
  y -= dy;
}

// Shifted derivative: y = z^(1-s) (z^s x)'
// Allows y = x.
SERIES_EXP(derivative_shift, y, (class A), (x,s), (const Series<A>& x, const int64_t s)) {
  const auto nk = x.known(), nz = x.nonzero();
  y.set_counts(nk, nz);
  for (int64_t i = 0; i < nz; i++)
    y[i] = (s + i) * x[i];
}

// Shifted integral: y = z^(-s) (z^(s-1) x)
// Allows y = x.
SERIES_EXP(integral_shift, y, (class A), (x,s), (const Series<A>& x, const int64_t s)) {
  const auto nk = x.known(), nz = x.nonzero();
  y.set_counts(nk, nz);
  if (!s && nz)
    y[0] = 0;
  for (int64_t i = !s; i < nz; i++)
    y[i] = x[i] / (s + i);
}

// Logarithm: y = log x
SERIES_EXP(log, y, (class A), (x), (const Series<A>& x)) {
  typedef remove_const_t<A> S;
  const auto n = x.known();
  if (!n) return y.set_unknown();
  slow_assert(!y.alias(x) && x.nonzero() && x[0] == 1);

  // log via y' = x'/x
  Series<S> dx(n);
  dx = derivative_shift(x, 0);
  y = div(dx, x);
  y = integral_shift(y, 0);
}

// Shifted logarithm: y = z^-s log(1 + z^s x)
SERIES_EXP(log1p, y, (class A), (x,s), (const Series<A>& x, const int64_t s)) {
  typedef remove_const_t<A> S;
  const auto n = x.known();
  if (!n) return y.set_unknown();
  slow_assert(!y.alias(x) && s > 0);

  // log via y' = x'/x
  Series<S> dx(n);
  dx = derivative_shift(x, s);
  y = div1p(dx, x, s);
  y = integral_shift(y, s);
}

// Exponential: y = e^x
SERIES_EXP(exp, y, (class A), (x), (const Series<A>& x)) {
  typedef remove_const_t<A> S;
  const auto n = x.known();
  if (!n) return y.set_unknown();
  slow_assert(!y.alias(x) && (!x.nonzero() || x[0] == 0));

  // Base case
  y.set_scalar(1, 1);

  // Newton step:
  //   f(y) = log y - x
  //   f'(y) = 1/y
  //   N(y) = y0 - f(y0)/f'(y)
  //        = y0 - (log(y0) - x)/(1/y)
  //        = y0 - y*(log(y0) - x)
  Series<S> dy(n);
  newton_iterate(1, n, [&x, &y, &dy](const int64_t m0, const int64_t m, const bool refine) {
    y.set_known(m);

    // dy = y*(log(y0) - x)
    static_assert(!is_interval<S>);  // Assume y = y0 for now
    dy = log(y);
    dy -= x;
    dy = mul(y, dy);

    // Update
    y.high_sub(m0, dy.high(m0));
  });
}

// Shifted exponential: y = z^-s (e^(az^s x) - 1)
SERIES_EXP(expm1, y, (class A), (x,a,s), (const Series<A>& x, const int a, const int64_t s)) {
  typedef remove_const_t<A> S;
  slow_assert(!y.alias(x) && abs(a) == 1 && s > 0);

  // Base case
  const auto n = x.known();
  if (!n) return y.set_unknown();
  y.set_scalar(1, x.nonzero() ? a > 0 ? x[0] : -x[0] : 0);

  // Newton step:
  //   y = z^-s (exp(az^s x) - 1)
  //   log1p(y, s) = z^-s log(1 + z^s (z^-s (exp(az^s x) - 1)))
  //               = z^-s log(1 + exp(az^s x) - 1)
  //               = z^-s az^s x
  //               = ax
  //   f(y) = log1p(y, s) - ax
  //   f'(y) = 1/(1 + z^s y)
  //   N(y) = y0 - f(y0) / f'(y)
  //        = y0 - (1 + z^s y)(log1p(y0, s) - ax)
  Series<S> dy(n), t(n);
  newton_iterate(1, n, [&x, &y, &dy, &t, a, s](const int64_t m0, const int64_t m, const bool refine) {
    y.set_known(m);

    // dy = (1 + z^s y)(log1p(y0, s) - ax)
    static_assert(!is_interval<S>);  // Assume y = y0 for now
    t = log1p(y, s);
    if (a > 0) t -= x;
    else t += x;
    dy = mul1p(t, y, s);

    // Update
    y.high_sub(m0, dy.high(m0));
  });
}

// Shifted log1p_exp: y = log1p(e^x, s) = z^-s log (1 + z^s e^x)
SERIES_EXP(log1p_exp, y, (class A), (x,s), (const Series<A>& x, const int64_t s)) {
  typedef remove_const_t<A> S;
  const auto n = x.known();
  if (!n) return y.set_unknown();
  slow_assert(!y.alias(x) && s > 0 && (!x.nonzero() || !x[0]));

  // Base case
  y.set_scalar(1, 1);

  // Newton step:
  //   y = z^-s log(1 + z^s e^x)
  //   e^(z^s y) = 1 + z^s e^x
  //   f(y) = e^(z^s y) - 1 - z^s e^x
  //   f'(y) = z^s e^(z^s y)
  //   N(y) = y0 - f(y0) / f'(y)
  //        = y0 - (e^(z^s y0) - 1 - z^s e^x) / (z^s e^(z^s y))
  //        = y0 - (z^-s (e^(z^s y0) - 1) - e^x) / e^(z^s y)
  //        = y0 - (z^-s (e^(z^s y0) - 1) - e^x) / e^(z^s y0) e^(z^s (y0-y))
  //        = y0 - e^(z^s (y0-y)) (z^-s (1 - e^(-z^s y0)) - e^(x-z^s y0))
  //        = y0 + e^(z^s (y0-y)) (expm1(-y0, s) + e^(x-z^s y0))
  Series<S> ndy(n), t(n);
  newton_iterate(1, n, [&x, &y, &ndy, &t, s](const int64_t m0, const int64_t m, const bool refine) {
    y.set_known(m);

    // dy = expm1(-y0, s) + e^(x-z^s y0)
    t = x.low(m);
    t.high_sub(s, y);
    ndy = exp(t);
    t = expm1(y, -1, s);
    ndy += t;

    // Update
    y.high_add(m0, ndy.high(m0));
  });
}

// Multiplication by a power of two: y = 2^k x
SERIES_EXP(ldexp, y, (class A), (x,k), (const Series<A>& x, const int k)) {
  const auto nk = x.known(), nz = x.nonzero();
  y.set_counts(nk, nz);
  for (int64_t i = 0; i < nz; i++)
    y[i] = ldexp(x[i], k);
}

}  // namespace mandelbrot
