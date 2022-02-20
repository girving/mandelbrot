// Series arithmetic
#pragma once

#include "arith.h"
#include "array.h"
#include "bit.h"
#include "codelets.h"
#include "cutil.h"
#include "debug.h"
#include "fft.h"
#include "is_interval.h"
#include "loops.h"
#include "noncopyable.h"
#include "preprocessor.h"
#include "print.h"
#include "span.h"
#include <algorithm>
#include <cstring>
#include <iostream>
#include <memory>
#include <type_traits>
namespace mandelbrot {

using std::add_const_t;
using std::conditional_t;
using std::enable_if_t;
using std::initializer_list;
using std::is_const_v;
using std::is_constructible_v;
using std::is_convertible_v;
using std::is_trivially_copyable_v;
using std::is_trivially_destructible_v;
using std::ldexp;
using std::max;
using std::min;
using std::move;
using std::numeric_limits;
using std::ostream;
using std::remove_const_t;
using std::swap;
using std::type_identity_t;
using std::unique_ptr;
struct Poly;
struct Sig;

template<class T, bool view = false> struct Series;
template<class T> using SeriesView = Series<T,true>;
template<class A> struct IsSeries { static constexpr bool value = false; };
template<class T, bool v> struct IsSeries<Series<T,v>> { static constexpr bool value = true; };
template<class F> struct SeriesExp;

// x += a
template<class T> void add_scalar(Series<T>& x, const typename Series<T>::Scalar a);

// self ±= sign z^s f for sign = ±1
template<class T> void high_addsub(Series<T>& y, const int sign, const int64_t s, SeriesView<add_const_t<T>> x);

template<class T, bool view_> struct Series : public conditional_t<view_,span<T>,Array<T>> {
  static_assert(!IsSeries<T>::value);  // Catch template bugs early
  typedef Undevice<remove_const_t<T>> Scalar;
  typedef Scalar value_type;
  typedef conditional_t<view_,span<T>,Array<T>> Base;
  using Base::data;
private:
  typedef Scalar S;
  typedef add_const_t<T> CT;

  int64_t known_;  // Known terms
  int64_t nonzero_;  // Possibly nonzero terms

  struct Unusable {};
  template<class A,bool v> friend struct Series;
public:

  Series() : known_(0), nonzero_(0) {}
  explicit Series(int64_t limit) : Base(limit), known_(0), nonzero_(0) {}
  Series(const Series&) = default;
  Series(int64_t limit, initializer_list<S>&& cs)
    : Series(limit) {
    slow_assert(cs.size() <= size_t(Base::size_));
    known_ = nonzero_ = cs.size();
    int64_t i = 0;
    for (const auto& c : cs)
      data()[i++] = c;
  }
  Series(initializer_list<S>&& cs) : Series(cs.size(), move(cs)) {}
  template<class U,bool r> Series(
      const Series<U,r>& g,
      enable_if_t<is_convertible_v<const typename Series<U,r>::Base&,Base>,Unusable> u = Unusable())
    : Base(g), known_(g.known_), nonzero_(g.nonzero_) {}
  Series(Series&& g) : Base(move(static_cast<Base&>(g))), known_(g.known_), nonzero_(g.nonzero_) {
    g.known_ = g.nonzero_ = 0;
  }
  ~Series() = default;
  void clear() { Base::clear(); known_ = nonzero_ = 0; }

  // Assignment
  void set_scalar(int64_t known, const S a) {
    known = relu(known);
    slow_assert(limit() >= (known > 0));
    known_ = known;
    nonzero_ = known > 0;
    if (known_) {
      // One slow write
      if constexpr (is_device<T>) single_host_to_device(data(), a);
      else data()[0] = a;
    }
  }
  void set_scalar(int64_t known, const int a) { set_scalar(known, S(a)); }
  void operator=(const Series& g) {
    slow_assert(limit() >= g.nonzero_);
    known_ = g.known_;
    nonzero_ = g.nonzero_;
    same_to_same(span(), g.span());
  }
  void operator=(conditional_t<is_const_v<T> && view_,Unusable,SeriesView<CT>> g) {
    slow_assert(limit() >= g.nonzero_);
    known_ = g.known_;
    nonzero_ = g.nonzero_;
    same_to_same(span(), g.span());
  }
  template<class F> void operator=(SeriesExp<F>&& e) { e.set(*this); }
  void swap(Series& g) {
    static_cast<Base&>(*this).swap(g);
    std::swap(known_, g.known_);
    std::swap(nonzero_, g.nonzero_);
  }

  // Adding or removing const
  const Series<const S>& const_() const { return *reinterpret_cast<const Series<const S>*>(this); }
  const Series<S>& const_cast_() const { return *reinterpret_cast<const Series<S>*>(this); }

  // Information
  int64_t known() const { return known_; }
  int64_t nonzero() const { return nonzero_; }
  int64_t limit() const { if constexpr (view_) return nonzero_; else return Base::size_; }
  bool valid(const int64_t i) const { return (uint64_t)i < (uint64_t)nonzero_; }
  T& operator[](const int64_t n) const { assert(valid(n)); return data()[n]; }
  bool alias(SeriesView<CT> f) const {
    static_assert(!view_);
    return data() <= f.data() && f.data() < data() + limit();
  }

  // Iteration
  typedef T* iterator;
  typedef T* const_iterator;
  T* begin() const { return data(); }
  T* end() const { return data() + nonzero_; }
  size_t size() const { return nonzero_; }

  // Verify that low terms vanish
  void assert_low_near_zero(const int64_t n) const {
    slow_assert(n <= known_);
    const auto nz = min(n, nonzero_);
    for (int64_t i = 0; i < nz; i++)
      slow_assert(bound(data()[i]) < 1e-6, "x = %.3g, x[%d] = %.3g != 0", *this, i, data()[i]);
  }

  // Output
  friend ostream& operator<<(ostream& out, const Series& f) {
    static_assert(!is_device<T>);
    return out << f.span();
  }

  // Change the number of terms in place
  void set_unknown() { known_ = nonzero_ = 0; }
  void set_counts(const int64_t known, const int64_t nonzero) {
    slow_assert(uint64_t(nonzero) <= uint64_t(min(limit(), known)));
    known_ = known;
    nonzero_ = nonzero;
  }
  void set_known(const int64_t known) {
    slow_assert(known >= 0);
    known_ = known;
    nonzero_ = min(nonzero_, known);
  }

  // Non-aliasing copy
  Series<remove_const_t<T>> copy(const int64_t limit) const {
    Series<remove_const_t<T>> g(limit);
    g = *this;
    return g;
  }

  // Span accessors
  std::span<T> span() const { return std::span<T>(data(), nonzero_); }
  std::span<T> low_span(int64_t n) const { return std::span<T>(data(), min(relu(n), nonzero_)); }
  std::span<T> high_span(int64_t n) const { return std::span<T>(data() + n, relu(nonzero_ - n)); }

  // Noncopying view
  SeriesView<CT> view() const {
    SeriesView<CT> h;
    static_cast<typename SeriesView<CT>::Base&>(h) = *this;
    h.known_ = known_;
    h.nonzero_ = nonzero_;
    return h;
  }

  // The low terms of a series, without copying
  SeriesView<CT> low(const int64_t n) const {
    const auto nk = min(relu(n), known_);
    const auto nz = min(nk, nonzero_);
    SeriesView<CT> h;
    h.known_ = nk;
    h.nonzero_ = nz;
    if (nz)
      static_cast<typename SeriesView<CT>::Base&>(h) = *this;
    return h;
  }

  // The high terms of a series, without copying
  SeriesView<CT> high(int64_t n) const {
    n = relu(n);
    const auto nk = relu(known_ - n);
    const auto nz = relu(nonzero_ - n);
    SeriesView<CT> h;
    h.known_ = nk;
    h.nonzero_ = nz;
    if (nz)
      static_cast<typename SeriesView<CT>::Base&>(h) = std::span<CT>(data() + n, nz);
    return h;
  }

  // In-place arithmetic
  void operator+=(const S a) { add_scalar(*this, a); }
  void operator-=(const S a) { add_scalar(*this, -a); }
  void operator+=(const int a) { add_scalar(*this, S(a)); }
  void operator-=(const int a) { add_scalar(*this, S(-a)); }

  void operator+=(SeriesView<CT> f) { high_addsub(*this, 1, 0, f); }
  void operator-=(SeriesView<CT> f) { high_addsub(*this, -1, 0, f); }
  void high_add(const int64_t s, SeriesView<CT> f) { high_addsub(*this, 1, s, f); }
  void high_sub(const int64_t s, SeriesView<CT> f) { high_addsub(*this, -1, s, f); }
};

// Unevaluated series computations
template<class F> struct SeriesExp { F set; };

template<class... Args> struct ScalarT;
template<class T, class... Rest> struct ScalarT<T,Rest...> : public ScalarT<Rest...> {};
template<class T, bool v, class... Rest> struct ScalarT<Series<T,v>,Rest...> {
  typedef typename Series<T,v>::Scalar type;
};
template<class T, bool v, class... Rest> struct ScalarT<const Series<T,v>&,Rest...> {
  typedef typename Series<T,v>::Scalar type;
};
template<class... Args> using Scalar = typename ScalarT<Args...>::type;

#define SERIES_EXP(name, y, Ts, xs, capture, args) \
  template<UNPAREN Ts> auto name(UNPAREN args) { \
    auto set = [UNPAREN capture](auto& y) { set_##name(y, UNPAREN xs); }; \
    return SeriesExp<decltype(set)>{move(set)}; \
  } \
  template<UNPAREN Ts, class Dst> void set_##name(Dst& y, UNPAREN args)

// Negation: y = -x
DEF_LOOP(neg_loop, n, i, (S* y, const S* x),
  y[i] = -x[i];)
SERIES_EXP(neg, y, (class A,bool v), (x), (x=x.view()), (const Series<A,v>& x)) {
  const auto nk = x.known(), nz = x.nonzero();
  y.set_counts(nk, nz);
  neg_loop(nz, y.data(), x.data());
}

// Multiplication: z = xy
SERIES_EXP(mul, z, (class A,class B,bool va,bool vb), (x,y), (x=x.view(),y=y.view()),
           (const Series<A,va>& x, const Series<B,vb>& y)) {
  const auto nk = min(x.known(), y.known());
  const auto nz = min(nk, relu(x.nonzero() + y.nonzero() - 1));
  z.set_counts(nk, nz);
  fft_mul(z.span(), x.low_span(nz), y.low_span(nz));
}

// Shifted multiplication: z = x(1 + z^s y)
template<class T> void mul1p_post(Series<T>& z, SeriesView<add_const_t<T>> x,
                                  const int64_t post, const int64_t s, const int64_t xnz);
SERIES_EXP(mul1p, z, (class A,class B,bool va,bool vb), (x,y,s), (x=x.view(),y=y.view(),s),
           (const Series<A,va>& x, const Series<B,vb>& y, const int64_t s)) {
  slow_assert(!z.alias(x) && s > 0);
  const auto nk = min(x.known(), y.known() + s);
  const auto nz = min(nk, x.nonzero() + y.nonzero() + s - 1);
  const auto post = min(nz, max(s, x.nonzero()));
  z.set_counts(nk, nz);
  fft_mul(z.high_span(s), x.low_span(nz-s), y.low_span(nz-s));
  mul1p_post(z, x, post, s, x.nonzero());
}

// Squaring: y = x^2
SERIES_EXP(sqr, y, (class A,bool v), (x), (x=x.view()), (const Series<A,v>& x)) {
  const auto nk = x.known();
  const auto nz = min(nk, relu(2*x.nonzero() - 1));
  y.set_counts(nk, nz);
  fft_sqr(y.span(), x.span());
}

// Set y[0] = exp(x[0]) on either CPU or GPU.  exp can reference args and x0.
#if CODELETS  // Make a simple base case so that codelets can run
#define SERIES_BASE(name, args, names, exp) \
  DEF_SERIAL(name##_base_serial, (S* ys, const S* xs, const int nx COMMA_UNPAREN args), \
    const S x0 = nx ? xs[0] : 0; \
    ys[0] = (exp);) \
  template<class T> void name##_base(Series<T>& y, type_identity_t<SeriesView<const T>> x COMMA_UNPAREN args) { \
    y.set_counts(1, 1); \
    name##_base_serial(y.data(), x.data(), x.nonzero() COMMA_UNPAREN names); \
  }
#else
#define SERIES_BASE(...)  // codelets.cc will generate these for us
#endif

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
SERIES_BASE(inv, (), (), inv(x0))
SERIES_EXP(inv, y, (class A,bool v), (x), (x=x.view()), (const Series<A,v>& x)) {
  typedef Scalar<Series<A>> S;
  typedef remove_const_t<A> DS;
  const auto n = x.known();
  if (!n) return y.set_unknown();
  slow_assert(!y.alias(x) && x.nonzero());
  inv_base(y, x);

  // Newton step:
  //   1/y = x
  //   f(y) = 1/y - x
  //   f'(y) = -1/y^2
  //   N(y) = y0 - f(y0) / f'(y)
  //        = y0 - (1/y0 - x) / (-1/y^2)
  //        = y0 - y0(x y0 - 1)(y/y0)^2
  Series<DS> dy(n);
  newton_iterate(y.known(), n, [&y, &x, &dy](const int64_t m0, const int64_t m, const bool refine) {
    y.set_known(m);

    // dy = y0(xy0-1)(y/y0)^2
    static_assert(!is_interval<S>);  // Ignore y/y0 for now
    dy = mul(x, y);
    dy -= 1;
    dy = mul(dy, y);

    // Update
    y.high_sub(m0, dy.high(m0));
  });
}

// Shifted reciprocal: y = z^-s (1 / (1 + z^s x) - 1)
SERIES_EXP(inv1p, y, (class A,bool v), (x,s), (x=x.view(),s), (const Series<A,v>& x, const int64_t s)) {
  typedef Scalar<Series<A>> S;
  typedef remove_const_t<A> DS;
  const auto n = x.known();
  if (!n) return y.set_unknown();
  slow_assert(!y.alias(x) && s > 0);

  // Base case:
  //   y = z^-s (1 / (1 + z^s x) - 1)
  //     = z^-s (1 - z^s x + O(z^2s) - 1)
  //     = z^-s (-z^s x + O(z^2s))
  //     = -x + O(z^s)
  y = neg(x.low(s));

  // Newton step:
  //   1/(1 + z^s y) = 1 + z^s x
  //   f(y) = 1/(1 + z^s y) - 1 - z^s x
  //   f'(y) = -z^s/(1 + z^s y)^2
  //   N(y) = y0 - f(y0) / f'(y)
  //        = y0 - (1/(1 + z^s y0) - 1 - z^s x) / (-z^s/(1 + z^s y)^2)
  //        = y0 + z^-s (1 + z^s y0)(1 - (1 + z^s y0) - z^s x(1 + z^s y0)) ((1 + z^s y)/(1 + z^s y0))^2
  //        = y0 - (1 + z^s y0)(y0 + x(1 + z^s y0)) ((1 + z^s y)/(1 + z^s y0))^2
  Series<DS> dy(n), t(n);
  newton_iterate(y.known(), n, [&y, &x, &dy, &t, s](const int64_t m0, const int64_t m, const bool refine) {
    y.set_known(m);

    // dy = y0(xy0-1)(y/y0)^2
    static_assert(!is_interval<S>);  // Ignore y/y0 for now
    t = mul1p(x, y, s);
    t += y;
    dy = mul1p(t, y, s);

    // Update
    y.high_sub(m0, dy.high(m0));
  });
}

// Division: y = a / b
SERIES_EXP(div, y, (class A,class B,bool va,bool vb), (a,b), (a=a.view(),b=b.view()),
           (const Series<A,va>& a, const Series<B,vb>& b)) {
  typedef Scalar<Series<A>> S;
  typedef remove_const_t<A> DS;
  slow_assert(!y.alias(a) && !y.alias(b));
  const auto n = min(a.known(), b.known());

  // Compute the inverse and multiply
  Series<DS> inv_b(n);
  inv_b = inv(b.low(n));
  y = mul(a, inv_b);

  // One more step of Newton refinement:
  //   y = a/b
  //   f(y) = by - a
  //   f'(y) = b
  //   N(y) = y0 - (b*y0 - a)/b
  //        = y0 - (b*y0 - a)(1/b)
  static_assert(!is_interval<S>);  // Assume y0 = y
  Series<DS> dy(n);
  dy = mul(b, y);
  dy -= a;
  dy = mul(dy, inv_b);
  y -= dy;
}

// Shifted division: y = a / (1 + z^s b)
SERIES_EXP(div1p, y, (class A,class B,bool va,bool vb), (a,b,s), (a=a.view(),b=b.view(),s),
           (const Series<A,va>& a, const Series<B,vb>& b, const int64_t s)) {
  typedef Scalar<Series<A>> S;
  typedef remove_const_t<A> DS;
  slow_assert(!y.alias(a) && !y.alias(b) && s > 0);
  const auto n = min(a.known(), b.known() + s);

  // Compute the inverse and multiply
  Series<DS> inv_b(n - s);
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
  Series<DS> t(n), dy(n);
  t = mul1p(y, b, s);
  t -= a;
  dy = mul1p(t, inv_b, s);
  y -= dy;
}

// Shifted derivative: y = z^(1-s) (z^s x)'
// Allows y = x.
DEF_LOOP(derivative_shift_loop, n, i, (S* y, const S* x, const int s),
  y[i] = (s + i) * x[i];)
SERIES_EXP(derivative_shift, y, (class A,bool v), (x,s), (x=x.view(),s), (const Series<A,v>& x, const int64_t s)) {
  const auto nk = x.known(), nz = x.nonzero();
  y.set_counts(nk, nz);
  derivative_shift_loop(nz, y.data(), x.data(), s);
}

// Shifted integral: y = z^(-s) (z^(s-1) x)
// Allows y = x.
DEF_LOOP(integral_shift_loop, n, i, (S* y, const S* x, const int s),
  y[i] = s + i ? x[i] / (s + i) : S(0);)
SERIES_EXP(integral_shift, y, (class A,bool v), (x,s), (x=x.view(),s), (const Series<A,v>& x, const int64_t s)) {
  const auto nk = x.known(), nz = x.nonzero();
  y.set_counts(nk, nz);
  integral_shift_loop(nz, y.data(), x.data(), s);
}

// Logarithm: y = log x
SERIES_EXP(log, y, (class A,bool v), (x), (x=x.view()), (const Series<A,v>& x)) {
  typedef remove_const_t<A> DS;
  const auto n = x.known();
  if (!n) return y.set_unknown();
  slow_assert(!y.alias(x) && x.nonzero());
  if constexpr (!codelets && !is_device<DS>) slow_assert(x[0] == 1);

  // log via y' = x'/x
  Series<DS> dx(n);
  dx = derivative_shift(x, 0);
  y = div(dx, x);
  y = integral_shift(y, 0);
}

// Shifted logarithm: y = z^-s log(1 + z^s x)
SERIES_EXP(log1p, y, (class A,bool v), (x,s), (x=x.view(),s), (const Series<A,v>& x, const int64_t s)) {
  typedef remove_const_t<A> DS;
  const auto n = x.known();
  if (!n) return y.set_unknown();
  slow_assert(!y.alias(x) && s > 0);

  // log via y' = x'/x
  Series<DS> dx(n);
  dx = derivative_shift(x, s);
  y = div1p(dx, x, s);
  y = integral_shift(y, s);
}

// Exponential: y = e^x
SERIES_BASE(exp, (), (), S(1))
SERIES_EXP(exp, y, (class A,bool v), (x), (x=x.view()), (const Series<A,v>& x)) {
  typedef Scalar<Series<A>> S;
  typedef remove_const_t<A> DS;
  const auto n = x.known();
  if (!n) return y.set_unknown();
  slow_assert(!y.alias(x));
  if constexpr (!codelets && !is_device<A>) slow_assert(!x.nonzero() || x[0] == 0);
  exp_base(y, x);

  // Newton step:
  //   f(y) = log y - x
  //   f'(y) = 1/y
  //   N(y) = y0 - f(y0)/f'(y)
  //        = y0 - (log(y0) - x)/(1/y)
  //        = y0 - y*(log(y0) - x)
  Series<DS> dy(n);
  newton_iterate(y.known(), n, [&x, &y, &dy](const int64_t m0, const int64_t m, const bool refine) {
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
SERIES_EXP(expm1, y, (class A,bool v), (x,a,s), (x=x.view(),a,s),
           (const Series<A,v>& x, const int a, const int64_t s)) {
  typedef Scalar<Series<A>> S;
  typedef remove_const_t<A> DS;
  slow_assert(!y.alias(x) && abs(a) == 1 && s > 0);

  // Base case:
  //   y = z^-s (e^(az^s x) - 1)
  //     = z^-s (1 + az^s x + O(z^2s) - 1)
  //     = z^-s (az^s x + O(z^2s))
  //     = ax + O(z^s)
  const auto n = x.known();
  if (a > 0) y = x.low(s);
  else y = neg(x.low(s));

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
  Series<DS> dy(n), t(n);
  newton_iterate(y.known(), n, [&x, &y, &dy, &t, a, s](const int64_t m0, const int64_t m, const bool refine) {
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
SERIES_EXP(log1p_exp, y, (class A,bool v), (x,s), (x=x.view(),s), (const Series<A,v>& x, const int64_t s)) {
  typedef remove_const_t<A> DS;
  const auto n = x.known();
  if (!n) return y.set_unknown();
  slow_assert(!y.alias(x) && s > 0);

  // Base case:
  //   y = z^-s log (1 + z^s e^x)
  //     = z^-s (z^s e^x + O(z^2s))
  //     = e^x + O(z^s)
  y = exp(x.low(s));

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
  Series<DS> ndy(n), t(n);
  newton_iterate(y.known(), n, [&x, &y, &ndy, &t, s](const int64_t m0, const int64_t m, const bool refine) {
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
DEF_LOOP(ldexp_loop, n, i, (S* y, const S* x, const int k),
  y[i] = ldexp(x[i], k);)
SERIES_EXP(ldexp, y, (class A,bool v), (x,k), (x=x.view(),k), (const Series<A,v>& x, const int k)) {
  const auto nk = x.known(), nz = x.nonzero();
  y.set_counts(nk, nz);
  ldexp_loop(nz, y.data(), x.data(), k);
}

// For test purposes
Series<double> approx(const Poly& x, const int64_t n);
double error(SeriesView<const double> x, SeriesView<const double> y, const bool relative = false);
double error(SeriesView<const double> x, initializer_list<double>&& ys, const bool relative = false);
double error(SeriesView<const double> x, const Poly& y, const bool relative = false);
#define ASSERT_TOL2(tol, x, y) { const auto e = error(x, y); ASSERT_LE(e, tol) << format("e %g, x %g, y %g", e, x, y); }
#define ASSERT_TOL(tol, x, ...) { const Series<double> _y({__VA_ARGS__}); ASSERT_TOL2(tol, x, _y); }
#define ASSERT_CLOSE(x, ...) ASSERT_TOL(3e-14, x, __VA_ARGS__)
#define ASSERT_CLOSE2(x, y) ASSERT_TOL2(3e-14, x, y)
#define ASSERT_EXACT(x, ...) ASSERT_TOL(0, x, __VA_ARGS__)
#define ASSERT_EXACT2(x, y) ASSERT_TOL2(0, x, y)

// Host to device and backwards
template<class S> void host_to_device(Series<Device<S>>& y, type_identity_t<SeriesView<const S>> x) {
  slow_assert(y.limit() >= x.nonzero());
  y.set_counts(x.known(), x.nonzero());
  host_to_device(y.span(), x.span());
}
template<class S> void device_to_host(Series<S>& y, type_identity_t<SeriesView<const Device<S>>> x) {
  slow_assert(y.limit() >= x.nonzero());
  y.set_counts(x.known(), x.nonzero());
  device_to_host(y.span(), x.span());
}
template<class T,bool v> conditional_t<is_device<T>,Series<typename Series<T,v>::Scalar>,const Series<T,v>&>
host_copy(const Series<T,v>& x) {
  if constexpr (is_device<T>) {
    typedef typename Series<T,v>::Scalar S;
    Series<S> hx(x.nonzero());
    device_to_host(hx, x);
    return hx;
  } else
    return x;
}

// Write a series to a file in a simple text format, or read it back
template<class T> void write_series(const string& path, const vector<string>& comments, SeriesView<const T> x);
template<class T> tuple<vector<string>,Series<T>> read_series(const string& path);

}  // namespace mandelbrot

// Pull in autogenerated codelets
#if !CODELETS
#include "gen-series-bases.h"
#endif
