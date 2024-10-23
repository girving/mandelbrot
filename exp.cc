// Symbolic expressions with CSE

#include "exp.h"
#include "bit.h"
#include "complex.h"
#include "join.h"
#include <variant>
namespace mandelbrot {

using std::holds_alternative;
using std::make_pair;
using std::move;
using std::nullopt;
using std::variant;

CSE* CSE::active_ = 0;

static string parens_if(const string& s, const bool p) { return p ? tfm::format("(%s)", s) : s; }

// Expression subtypes
// Precedence follows https://en.wikipedia.org/wiki/Operators_in_C_and_C++#Operator_precedence
#define EXP(prec_, args_, arbitrary_) \
  static constexpr int prec = prec_; \
  span<Exp> args() { return args_; } \
  span<const Exp> args() const { return args_; } \
  Sig arbitrary() const { using mandelbrot::arbitrary; return arbitrary_; }
#define SIMPLE(prec_, args_, arbitrary_, str) \
  EXP(prec_, args_, arbitrary_) \
  string show(const int need_prec) const { return parens_if(str, prec > need_prec); }
struct VarExp { string v; Sig sig; SIMPLE(0, span<Exp>(), sig, v) };
struct IntExp { int n; SIMPLE(0, span<Exp>(), Sig(n), tfm::format("%d", n)) };
struct FieldExp { Exp x; const char* f; SIMPLE(2, (span{&x, 1}), arbitrary(f, x.sig()), tfm::format("%s.%s", x.show(2), f)) };
struct NegExp { Exp x; SIMPLE(3, (span{&x, 1}), -x.sig(), tfm::format("-%s", x.show(3))) };
struct AddExp { Exp x, y; SIMPLE(6, (span{&x,2}), arbitrary(x.sig(), y.sig()), tfm::format("%s + %s", x.show(6), y.show(5))) };
struct SubExp { Exp x, y; SIMPLE(6, (span{&x,2}), arbitrary(x.sig(), y.sig()), tfm::format("%s - %s", x.show(6), y.show(5))) };
struct MulExp {
  Exp x, y; EXP(5, (span{&x, 2}), arbitrary(x.sig(), y.sig()))
  string show(const int need_prec) const {
    if (x.sig() == y.sig()) return tfm::format("sqr(%s)", x);
    if (x.two()) return tfm::format("twice(%s)", y);
    if (y.two()) return tfm::format("twice(%s)", x);
    return parens_if(tfm::format("%s * %s", x.show(5), y.show(4)), prec > need_prec);
  }
};
struct DivExp {
  Exp x; int a; EXP(5, (span{&x, 1}), arbitrary(x.sig(), Sig(a)))
  string show(const int need_prec) const {
    if (a == 2) return tfm::format("half(%s)", x);
    if (a > 2 && has_single_bit(unsigned(a))) return tfm::format("ldexp(%s, -%d)", x, int(countr_zero(unsigned(a))));
    return parens_if(tfm::format("%s / %s", x.show(5), a), prec > need_prec);
  }
};
template<int n> struct CallExp {
  const char* f; Exp xs[n];
  static constexpr int prec = 2;
  span<Exp> args() { return xs; }
  span<const Exp> args() const { return xs; }
  Sig arbitrary() const {
    Sig s[n];
    for (int i = 0; i < n; i++) s[i] = xs[i].sig();
    return mandelbrot::arbitrary(f, s);
  }
  string show(const int need_prec) const { return tfm::format("%s(%s)", f, join(xs)); }
};
struct OtherExp {
  string s; int prec; Sig sig;
  span<Exp> args() const { return span<Exp>(); }
  Sig arbitrary() const { return sig; }
  string show(const int need_prec) const { return parens_if(s, prec > need_prec); }
};

// Store expression details inside a shared_ptr<const ExpData> for cheap copying
struct ExpData {
  variant<VarExp,IntExp,FieldExp,NegExp,AddExp,SubExp,MulExp,DivExp,CallExp<1>,CallExp<2>,CallExp<3>,OtherExp> exp;
  Sig sig;
};

Exp::Exp(const int a) : d(new ExpData{IntExp{a}, Sig(a)}) {}
Exp::Exp(const string& x, const Sig sig) : d(new ExpData{VarExp{x, sig}, sig}) {}
int Exp::prec() const { return visit([](const auto& e) { return e.prec; }, d->exp); }
bool Exp::fast() const { return d->exp.index() < 4; }  // Don't CSE VarExp, IntExp, FieldExp, or NegExp
const Sig& Exp::sig() const { return d->sig; }

template<class E> const E* Exp::get() const {
  return holds_alternative<E>(d->exp) ? &std::get<E>(d->exp) : 0;
}

optional<int> Exp::unint() const {
  optional<int> n;
  if (const auto* e = get<IntExp>()) n = e->n;
  return n;
}

optional<Exp> Exp::unneg() const {
  optional<Exp> r;
  const auto n = unint();
  if (n && *n < 0) r = Exp(-*n);
  else if (const auto* y = get<NegExp>()) r = y->x;
  return r;
}

string Exp::show(const int need_prec) const {
  return visit([need_prec](const auto& e) { return e.show(need_prec); }, d->exp);
}

span<const Exp> Exp::args() const {
  return visit([](const auto& e) { return span<const Exp>(e.args()); }, d->exp);
}

Exp Exp::map_args(const function<Exp(const Exp&)>& f) const {
  const auto sig = this->sig();
  return visit([&f,sig](const auto& e) {
    auto me = e;
    for (Exp& x : me.args())
      x = f(x);
    return Exp(shared_ptr<const ExpData>(new ExpData{me, sig}));
  }, d->exp);
}

CSE::CSE(const bool assume_exact) : assume_exact(assume_exact) { slow_assert(!active_); active_ = this; }
CSE::~CSE() { active_ = 0; }

Exp CSE::cse(const Exp& e) {
  if (e.d->exp.index() < 4) return e;  // Don't CSE VarExp, IntExp, FieldExp, or NegExp
  if (const auto small = unsmall(e.sig())) return Exp(*small);

  // Try to pull the expression out of our CSE database
  const auto it = signatures.find(e.sig());
  if (it != signatures.end()) return it->second;

  // Not found, so add it in
  signatures.insert(make_pair(e.sig(), e));
  return e;
}

template<class E> Exp CSE::cse(const E& e, Sig sig) {
  slow_assert(active_);
  auto& cse = *active_;
  if (!cse.assume_exact) sig = e.arbitrary();
  return cse.cse(Exp(shared_ptr<const ExpData>(new ExpData{move(e), sig})));
}

#define SUB(E) template Exp CSE::cse(const E&, Sig);
SUB(VarExp) SUB(IntExp) SUB(FieldExp) SUB(NegExp) SUB(AddExp) SUB(SubExp) SUB(MulExp) SUB(DivExp)
SUB(CallExp<1>) SUB(CallExp<2>) SUB(CallExp<3>) SUB(OtherExp)
#undef SUB

Exp operator-(const Exp& e) {
  if (const auto n = e.unint()) return Exp(-*n);
  if (const auto x = e.unneg()) return *x;
  return CSE::cse(NegExp{e}, -e.sig());
}

vector<Exp> operator-(span<const Exp> xs) {
  vector<Exp> ys;
  for (const auto& x : xs)
    ys.push_back(-x);
  return ys;
}

Exp operator+(const Exp& x, const Exp& y) {
  if (x.zero()) return y;
  if (y.zero()) return x;
  if (const auto nx = x.unneg()) return y - *nx;
  if (const auto ny = y.unneg()) return x - *ny;
  return CSE::cse(AddExp{x, y}, x.sig() + y.sig());
}

Exp operator-(const Exp& x, const Exp& y) {
  if (x.zero()) return -y;
  if (y.zero()) return x;
  if (const auto nx = x.unneg()) return -(*nx + y);
  if (const auto ny = y.unneg()) return x + *ny;
  return CSE::cse(SubExp{x, y}, x.sig() - y.sig());
}

Exp operator*(const Exp& x, const Exp& y) {
  if (x.zero() || y.zero()) return 0;
  if (const auto nx = x.unneg()) return -(*nx * y);
  if (const auto ny = y.unneg()) return -(x * *ny);
  if (x.one()) return y;
  if (y.one()) return x;
  return CSE::cse(MulExp{x, y}, x.sig() * y.sig());
}

Exp operator/(const Exp& x, const int a) {
  slow_assert(a);
  if (a == 1) return x;
  if (a == -1) return -x;
  if (a < 0) return -(x / -a);
  if (const auto nx = x.unneg()) return -(*nx / a);
  return CSE::cse(DivExp{x, a}, x.sig() * inv(Sig(a)));
}

Exp call(const char* f, const Exp& x0, const Sig sig) { return CSE::cse(CallExp<1>{f, {x0}}, sig); }
Exp call(const char* f, const Exp& x0, const Exp& x1, const Sig sig) { return CSE::cse(CallExp<2>{f, {x0, x1}}, sig); }
Exp call(const char* f, const Exp& x0, const Exp& x1, const Exp& x2, const Sig sig) {
  return CSE::cse(CallExp<3>{f, {x0, x1, x2}}, sig);
}

Exp other(const string& s, const int prec, const Sig sig) { return CSE::cse(OtherExp{s, prec, sig}, sig); }

Exp fma(const Exp& x, const Exp& y, const Exp& s) {
  return call("fma", x, y, s, random_sig());
}

Exp sum(span<const Exp> xs) {
  const size_t n = xs.size();
  if (n == 0) return 0;
  if (n == 1) return xs[0];
  return sum(xs.first(n/2)) + sum(xs.last((n+1)/2));
}

Exp inv(const Exp& x) { return call("inv", x, inv(x.sig())); }

Exp Exp::field(const char* f) const {
  FieldExp e{*this, f};
  return CSE::cse(e, e.arbitrary());
}

Complex<Exp> split(const Exp& e) {
  return Complex<Exp>(e.field("r"), e.field("i"));
}

vector<Exp> Exp::grad() const {
  static const Exp unknown("unknown", random_sig());
  const auto grad = [](const auto& e) {
    #define GRAD(T, ...) \
      if constexpr (is_same_v<decltype(e),const T&>) { \
        Exp g[] = {__VA_ARGS__}; return vector<Exp>(g+0, g+sizeof(g)/sizeof(Exp)); \
      } else
    GRAD(VarExp)
    GRAD(IntExp)
    GRAD(FieldExp, unknown)
    GRAD(NegExp, -1)
    GRAD(AddExp, 1, 1)
    GRAD(SubExp, 1, -1)
    GRAD(MulExp, e.y, e.x)
    GRAD(DivExp, unknown)  // We don't need this, so be lazy
    GRAD(CallExp<1>, unknown)
    GRAD(CallExp<2>, unknown, unknown)
    GRAD(CallExp<3>, unknown, unknown, unknown)
    GRAD(OtherExp)
    static_assert(sizeof(e) == 1000);
  };
  return visit(grad, d->exp);
}

void Stats::add(const Exp& e) {
  const auto add = [this](const auto& e) {
    #define ADD(T, ...) if constexpr (is_same_v<decltype(e),const T&>) { __VA_ARGS__; } else
    ADD(VarExp)
    ADD(IntExp)
    ADD(FieldExp)
    ADD(NegExp, negs++)
    ADD(AddExp, adds++)
    ADD(SubExp, adds++)
    ADD(MulExp, muls++)
    ADD(DivExp, others++)
    ADD(CallExp<1>, calls++)
    ADD(CallExp<2>, calls++)
    ADD(CallExp<3>, calls++)
    ADD(OtherExp)
    static_assert(sizeof(e) == 1000);
  };
  visit(add, e.d->exp);
}

string Stats::show() const {
  int total = 0;
  vector<string> terms;
  const auto p = [&total, &terms](const char* name, const int n) {
    total += n;
    if (n) terms.push_back(tfm::format("%d %s%s", n, name, n == 1 ? "" : "s"));
  };
  #define P(n) p(#n, n##s);
  P(add) P(mul) P(neg) P(call) P(other)
  #undef P
  return tfm::format("%d op%s = %s", total, total == 1 ? "" : "s", join(terms, " + "));
}

}  // namespace mandelbrot
