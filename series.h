// Series arithmetic
#pragma once

#include "debug.h"
#include "fft.h"
#include <algorithm>
#include <bit>
#include <iostream>
#include <memory>
#include <span>
#include <type_traits>
namespace mandelbrot {

using std::bit_ceil;
using std::conditional_t;
using std::countr_zero;
using std::initializer_list;
using std::is_const_v;
using std::is_trivially_copyable_v;
using std::is_trivially_destructible_v;
using std::ldexp;
using std::max;
using std::min;
using std::move;
using std::ostream;
using std::remove_const_t;
using std::shared_ptr;
using std::span;
using std::swap;
template<class T> struct Series;
template<class F> struct SeriesExp;

// For assertions that scalars aren't intervals
template<class T> struct IsIntervalT;
template<class T> struct IsIntervalT<const T> : public IsIntervalT<T> {};
template<> struct IsIntervalT<double> { static constexpr bool value = false; };
template<class T> static constexpr bool is_interval = IsIntervalT<T>::value;

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
  int64_t terms_;  // Known terms
  int64_t limit_;  // Allocated terms

  struct Unusable {};
  template<class A> friend struct Series;
public:

  Series() : terms_(0), limit_(0) {}
  explicit Series(int64_t limit)
    : terms_(0), limit_(max(limit, int64_t(0))) {
    if (limit_)
      x.reset(static_cast<S*>(malloc(limit_ * sizeof(S))), [](auto p) { free(p); });
  }
  Series(int64_t limit, initializer_list<S>&& cs)
    : limit_(limit) {
    slow_assert(cs.size() <= size_t(limit_));
    terms_ = cs.size();
    if (limit_) {
      shared_ptr<S[]> y(static_cast<S*>(malloc(limit_ * sizeof(S))), [](auto p) { free(p); });
      int64_t i = 0;
      for (const auto& c : cs)
        y[i++] = c;
      x = move(y);
    }
  }
  Series(initializer_list<S>&& cs) : Series(cs.size(), move(cs)) {}
  Series(const Series& g) = default;
  Series(const Series<conditional_t<is_const_v<T>,S,Unusable>>& g)
    : x(g.x), terms_(g.terms()), limit_(g.limit()) {}
  Series(Series&& g) : x(move(g.x)), terms_(g.terms_), limit_(g.limit_) { g.terms_ = g.limit_ = 0; }
  ~Series() = default;
  void clear() { x.reset(); terms_ = limit_ = 0; }

  // Assignment
  void operator=(const int a) {
    slow_assert(limit_);
    terms_ = max(terms_, int64_t(1));
    x[0] = a;
    for (int64_t i = 1; i < terms_; i++) x[i] = 0;
  }
  void operator=(const S& a) {
    slow_assert(limit_);
    terms_ = max(terms_, int64_t(1));
    x[0] = a;
    for (int64_t i = 1; i < terms_; i++) x[i] = 0;
  }
  void operator=(const Series& g) {
    slow_assert(limit_ >= g.limit_); terms_ = g.terms_; memcpy(x.get(), g.x.get(), terms_ * sizeof(S));
  }
  void operator=(const Series<conditional_t<is_const_v<T>,Unusable,const S>>& g) {
    slow_assert(limit_ >= g.limit_); terms_ = g.terms_; memcpy(x.get(), g.x.get(), terms_ * sizeof(S));
  }
  template<class F> void operator=(SeriesExp<F>&& e) { e.set(*this); }
  void swap(Series<S>& g) { using std::swap; swap(x, g.x); swap(terms_, g.terms_); swap(limit_, g.limit_); }

  // Adding or removing const
  const Series<const S>& const_() const { return *reinterpret_cast<const Series<const S>*>(this); }
  const Series<S>& const_cast_() const { return *reinterpret_cast<const Series<S>*>(this); }

  // Information
  int64_t terms() const { return terms_; }
  size_t size() const { return terms_; }
  int64_t limit() const { return limit_; }
  bool valid(const int64_t i) const { return (uint64_t)i < (uint64_t)terms_; }
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
  T* end() const { return x.get() + terms_; }

  // Verify that low terms vanish
  void assert_low_near_zero(const int64_t n) const {
    slow_assert(n <= terms_);
    for (int64_t i = 0; i < n; i++)
      slow_assert(near_zero(x[i]), "x = %.3g, x[%d] = %.3g != 0", *this, i, x[i]);
  }

  // Output
  friend ostream& operator<<(ostream& out, const Series& f) { return out << f.span(); }

  // Change the number of terms in place
  void set_empty() { terms_ = 0; }
  void set_terms(const int64_t n) { slow_assert(n <= limit_); terms_ = n; }
  void truncate(const int64_t n) { terms_ = min(terms_, max(n, int64_t(0))); }
  void extend(const int64_t n) {
    slow_assert(uint64_t(n) <= uint64_t(limit_));
    for (int64_t i = terms_; i < n; i++)
      x[i] = 0;
    terms_ = n;
  }
  Series<S> copy(const int64_t limit) const {
    Series<S> g(limit);
    g = *this;
    return g;
  }

  // Span accessors
  std::span<T> span() const { return std::span<T>(x.get(), terms_); }
  std::span<T> low_span(int64_t n) const { return std::span<T>(x.get(), min(max(int64_t(0), n), terms_)); }
  std::span<T> high_span(int64_t n) const { return std::span<T>(x.get() + n, max(int64_t(0), terms_ - n)); }

  // The low terms of a series, without copying
  Series<const S> low(int64_t n) const {
    n = max(n, int64_t(0));
    n = min(n, terms_);
    Series<const S> h;
    if (n > 0) {
      h.x = x;
      h.terms_ = h.limit_ = n;
    }
    return h;
  }
  Series<S> low(int64_t n) { return move(const_().low(n).const_cast_()); }

  // The high terms of a series, without copying
  Series<const S> high(int64_t n) const {
    Series<const S> h;
    if (n < terms_) {
      n = max(n, int64_t(0));
      h.x = shared_ptr<const S[]>(x, x.get() + n);
      h.terms_ = h.limit_ = terms_ - n;
    }
    return h;
  }
  Series<S> high(int64_t n) { return move(const_().high(n).const_cast_()); }

  // In-place arithmetic
  void operator+=(const int a) const { slow_assert(terms_ >= 1); x[0] += a; }
  void operator-=(const int a) const { slow_assert(terms_ >= 1); x[0] -= a; }
  void operator+=(const S& a) const { slow_assert(terms_ >= 1); x[0] += a; }
  void operator-=(const S& a) const { slow_assert(terms_ >= 1); x[0] -= a; }

  void operator+=(const Series<const S>& f) {
    terms_ = min(terms_, f.terms_);
    for (int64_t i = 0; i < terms_; i++)
      x[i] += f.x[i];
  }

  void operator-=(const Series<const S>& f) {
    terms_ = min(terms_, f.terms_);
    for (int64_t i = 0; i < terms_; i++)
      x[i] -= f.x[i];
  }
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
  const auto n = min(x.terms(), y.terms());
  z.set_terms(n);
  fft_mul<S>(z.span(), x.low_span(n), y.low_span(n));
}

// Shifted multiplication: z = x(1 + z^s y)
SERIES_EXP(mul1p, z, (class A,class B), (x,y,s), (const Series<A>& x, const Series<B>& y, const int64_t s)) {
  typedef remove_const_t<A> S;
  slow_assert(!z.alias(x) && s > 0);
  const auto n = min(x.terms(), y.terms() + s);
  z.set_terms(n);
  fft_mul<S>(z.high_span(s), x.low_span(n-s), y.low_span(n-s));
  for (int64_t i = s; i < n; i++)
    z[i] += x[i];
  const auto sn = min(s, n);
  for (int64_t i = 0; i < sn; i++)
    z[i] = x[i];
}

// Squaring: y = x^2
SERIES_EXP(sqr, y, (class A), (x), (const Series<A>& x)) {
  typedef remove_const_t<A> S;
  const auto n = x.terms();
  y.set_terms(n);
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
  for (int k = newton_steps(n0, n)-1; k >= -1; k--) {
    // We iterate down to k = -1 to do one final step of refinement with n0 = m = n
    const bool refine = k < 0;
    const auto m = refine ? n : (n + (int64_t(1)<<k) - 1) >> k;
    step(n0, m, refine);
    n0 = m;
  }
}

// Reciprocal: y = 1 / x
SERIES_EXP(inv, y, (class T), (x), (const Series<T>& x)) {
  typedef remove_const_t<T> S;
  slow_assert(!y.alias(x));

  // Base case
  const auto n = x.terms();
  if (!n) return y.set_empty();
  y = 1 / x[0];

  // Newton step:
  //   1/y = x
  //   f(y) = 1/y - x
  //   f'(y) = -1/y^2
  //   N(y) = y0 - f(y0) / f'(y)
  //        = y0 - (1/y0 - x) / (-1/y^2)
  //        = y0 - y0(x y0 - 1)(y/y0)^2
  Series<S> dy(n);
  newton_iterate(1, n, [&y, &x, &dy](const int64_t m0, const int64_t m, const bool refine) {
    y.extend(m);

    // dy = y0(xy0-1)(y/y0)^2
    static_assert(!is_interval<T>);  // Ignore y/y0 for now
    dy = mul(x, y);
    dy -= 1;
    dy = mul(dy, y);

    if (refine)
      y -= dy;
    else { // Extend
      dy.assert_low_near_zero(m0);
      y.high(m0) -= dy.high(m0);
    }
  });
}

// Shifted reciprocal: y = z^-s (1 / (1 + z^s x) - 1)
SERIES_EXP(inv1p, y, (class T), (x, s), (const Series<T>& x, const int64_t s)) {
  typedef remove_const_t<T> S;
  slow_assert(!y.alias(x) && s > 0);

  // Base case
  const auto n = x.terms();
  if (!n) return y.set_empty();
  y = -x[0];

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
    y.extend(m);

    // dy = y0(xy0-1)(y/y0)^2
    static_assert(!is_interval<T>);  // Ignore y/y0 for now
    t = mul1p(x, y, s);
    t += y;
    dy = mul1p(t, y, s);

    if (refine)
      y -= dy;
    else { // Extend
      dy.assert_low_near_zero(m0);
      y.high(m0) -= dy.high(m0);
    }
  });
}

// Division: y = a / b
SERIES_EXP(div, y, (class A,class B), (a,b), (const Series<A>& a, const Series<B>& b)) {
  typedef remove_const_t<A> S;
  slow_assert(!y.alias(a) && !y.alias(b));
  const auto n = min(a.terms(), b.terms());

  // Compute the inverse and multiply
  Series<S> inv_b(n);
  inv_b = inv(b.low(n));
  y = mul(a, inv_b);

  // One more step of Newton refinement
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
  const auto n = min(a.terms(), b.terms() + s);

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
  const auto n = x.terms();
  y.set_terms(n);
  for (int64_t i = 0; i < n; i++)
    y[i] = x[i] * (s + i);
}

// Shifted integral: y = z^(-s) (z^(s-1) x)
// Allows y = x.
SERIES_EXP(integral_shift, y, (class A), (x,s), (const Series<A>& x, const int64_t s)) {
  const auto n = x.terms();
  y.set_terms(n);
  if (!s && n > 0)
    y[0] = 0;
  for (int64_t i = !s; i < n; i++)
    y[i] = x[i] / (s + i);
}

// Logarithm: y = log x
SERIES_EXP(log, y, (class A), (x), (const Series<A>& x)) {
  typedef remove_const_t<A> S;
  const auto n = x.terms();
  if (!n) return y.set_empty();
  slow_assert(!y.alias(x));
  slow_assert(x[0] == 1);

  // log via y' = x'/x
  Series<S> dx(n);
  dx = derivative_shift(x, 0);
  y = div(dx, x);
  y = integral_shift(y, 0);
}

// Shifted logarithm: y = z^-s log(1 + z^s x)
SERIES_EXP(log1p, y, (class A), (x,s), (const Series<A>& x, const int64_t s)) {
  typedef remove_const_t<A> S;
  const auto n = x.terms();
  if (!n) return y.set_empty();
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
  const auto n = x.terms();
  if (!n) return y.set_empty();
  slow_assert(!y.alias(x) && x[0] == 0);

  // Base case
  y = 1;

  // Newton step:
  //   f(y) = log y - x
  //   f'(y) = 1/y
  //   N(y) = y0 - f(y0)/f'(y)
  //        = y0 - (log(y0) - x)/(1/y)
  //        = y0 - y*(log(y0) - x)
  Series<S> dy(n);
  newton_iterate(1, n, [&x, &y, &dy](const int64_t m0, const int64_t m, const bool refine) {
    y.extend(m);

    // dy = y*(log(y0) - x)
    static_assert(!is_interval<S>);  // Assume y = y0 for now
    dy = log(y);
    dy -= x;
    dy = mul(y, dy);

    if (refine)
      y -= dy;
    else {  // Expand
      dy.assert_low_near_zero(m0);
      y.high(m0) -= dy.high(m0);
    }
  });
}

// Shifted exponential: y = z^-s (e^(az^s x) - 1)
SERIES_EXP(expm1, y, (class A), (x,a,s), (const Series<A>& x, const int a, const int64_t s)) {
  typedef remove_const_t<A> S;
  slow_assert(!y.alias(x) && abs(a) == 1 && s > 0);

  // Base case
  const auto n = x.terms();
  if (!n) return y.set_empty();
  y = a * x[0];

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
    y.extend(m);

    // dy = (1 + z^s y)(log1p(y0, s) - ax)
    static_assert(!is_interval<S>);  // Assume y = y0 for now
    t = log1p(y, s);
    if (a > 0) t -= x;
    else t += x;
    dy = mul1p(t, y, s);

    if (refine)
      y -= dy;
    else {  // Expand
      dy.assert_low_near_zero(m0);
      y.high(m0) -= dy.high(m0);
    }
  });
}

// Shifted log1p_exp: y = log1p(e^x, s) = z^-s log (1 + z^s e^x)
SERIES_EXP(log1p_exp, y, (class A), (x,s), (const Series<A>& x, const int64_t s)) {
  typedef remove_const_t<A> S;
  const auto n = x.terms();
  slow_assert(!y.alias(x) && s > 0 && (!n || x[0] == 0));

  // Base case
  if (!n) return y.set_empty();
  y = 1;

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
    y.extend(m);

    // dy = expm1(-y0, s) + e^(x-z^s y0)
    t = x.low(m);
    t.high(s) -= y;
    ndy = exp(t);
    t = expm1(y, -1, s);
    ndy += t;

    if (refine)
      y += ndy;
    else {  // Expand
      ndy.assert_low_near_zero(m0);
      y.high(m0) += ndy.high(m0);
    }
  });
}

// Multiplication by a power of two: y = 2^k x
SERIES_EXP(ldexp, y, (class A), (x,k), (const Series<A>& x, const int k)) {
  const auto n = x.terms();
  y.set_terms(n);
  for (int64_t i = 0; i < n; i++)
    y[i] = ldexp(x[i], k);
}

}  // namespace mandelbrot
