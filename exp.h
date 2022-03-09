// Symbolic expressions with CSE
#pragma once

#include "arith.h"
#include "debug.h"
#include "format.h"
#include "is_interval.h"
#include "noncopyable.h"
#include "sig.h"
#include <functional>
#include <memory>
#include <optional>
#include <type_traits>
#include <unordered_map>
#include <vector>
namespace mandelbrot {

using std::function;
using std::optional;
using std::ostream;
using std::shared_ptr;
using std::unordered_map;
using std::vector;

struct CSE;
struct Exp;
struct ExpData;
struct Stats;
template<> struct IsIntervalT<Exp> { static constexpr bool value = false; };

struct Exp {
private:
  shared_ptr<const ExpData> d;
  Exp(const shared_ptr<const ExpData>& d) : d(d) {}
  struct Zero {};
  friend struct CSE;
  friend struct Stats;
public:

  Exp() : Exp(0) {}
  Exp(const Zero*) : Exp(0) {}
  Exp(const int a);
  Exp(const string& x, const Sig sig);
  Exp(const Exp& e) : d(e.d) {}
  Exp(Exp&& e) : d(e.d) {}  // Stay valid on move
  Exp& operator=(const Exp& f) { d = f.d; return *this; }

  optional<int> unint() const;
  bool zero() const { return unint() == 0; }
  bool one() const { return unint() == 1; }
  bool two() const { return unint() == 2; }
  int prec() const;  // Precedence follows https://en.wikipedia.org/wiki/Operators_in_C_and_C++#Operator_precedence
  bool fast() const;  // Fast enough to not CSE
  const Sig& sig() const;
  optional<Exp> unneg() const;
  Exp field(const char* f) const;
  explicit operator bool() const { return !zero(); }
  bool operator==(const Exp& e) const { return sig() == e.sig(); }
  string show(const int need_prec) const;
  span<const Exp> args() const;
  vector<Exp> grad() const;  // Gradient w.r.t. args()
  Exp map_args(const function<Exp(const Exp&)>& f) const;
  friend ostream& operator<<(ostream& out, const Exp& e) { return out << e.show(100); }
  template<class E> const E* get() const;
};

// Operation counts
struct Stats {
  int negs = 0, adds = 0, muls = 0, calls = 0, others = 0;
  void add(const Exp& e);
  string show() const;
};

// Common subexpression elimination
struct CSE : public Noncopyable {
private:
  static CSE* active_;
  const bool assume_exact;  // Whether to assume arithmetic is exact for CSE purposes
  unordered_map<Sig,Exp,SigHash> signatures;
  Exp cse(const Exp& e);
public:
  CSE(const bool assume_exact);
  ~CSE();

  // CSE a potential new expression
  template<class E> static Exp cse(const E& e, Sig sig);
};

// Arithmetic
Exp operator-(const Exp& e);
vector<Exp> operator-(span<const Exp> xs);
Exp operator+(const Exp& x, const Exp& y);
Exp operator-(const Exp& x, const Exp& y);
Exp operator*(const Exp& x, const Exp& y);
Exp operator/(const Exp& x, const int a);
Exp fma(const Exp& x, const Exp& y, const Exp& s);
Exp sum(span<const Exp> xs);
Exp inv(const Exp& x);
Complex<Exp> split(const Exp& e);

// Function calls
Exp call(const char* f, const Exp& x0, const Sig sig);
Exp call(const char* f, const Exp& x0, const Exp& x1, const Sig sig);
Exp call(const char* f, const Exp& x0, const Exp& x1, const Exp& x2, const Sig sig);
Exp other(const string& s, const int prec, const Sig sig);

}  // namespace mandelbrot
