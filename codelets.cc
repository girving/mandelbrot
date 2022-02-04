// Code generation

#include "arith.h"
#include "debug.h"
#include "format.h"
#include "noncopyable.h"
#include "print.h"
#include <array>
#include <deque>
#include <fstream>
#include <memory>
#include <optional>
#include <tuple>
#include <type_traits>
#include <unordered_set>
#include <vector>
namespace mandelbrot {

using std::any_of;
using std::array;
using std::common_type_t;
using std::cout;
using std::deque;
using std::endl;
using std::move;
using std::nullopt;
using std::ofstream;
using std::optional;
using std::ostream;
using std::string_view;
using std::tie;
using std::tuple;
using std::unique_ptr;
using std::unordered_set;
using std::vector;

// Container utilities
template<class P> auto& deref(const P& p) { slow_assert(p); return *p; }
template<class C> auto pop_front(C& xs) { auto x = xs.front(); xs.pop_front(); return x; }
template<class C> auto pop_back(C& xs) { auto x = xs.back(); xs.pop_back(); return x; }
template<class C,class A> bool contains(const C& xs, const A& x) { return xs.find(x) != xs.end(); }

// Concatenate a bunch of containers into a vector
template<class... Args> auto concat(const Args&... args) {
  vector<common_type_t<typename Args::value_type...>> xs;
  (xs.insert(xs.end(), args.begin(), args.end()), ...);
  return xs;
}

// String utilities
template<class C> string join(const C& ss, const string& sep = ", ") {
  string j;
  for (const auto& s : ss) {
    if (j.size()) j += sep;
    j += format("%s", s);
  }
  return j;
}
bool startswith(const string_view s, const string_view start) { return s.substr(0, start.size()) == start; }
bool endswith(const string_view s, const string_view end) { return s.substr(relu(s.size() - end.size())) == end; }

// Function object for nonzeroness
struct Nonzero { template<class T> bool operator()(T&& x) const { return bool(x); } };
constexpr Nonzero nonzero;
template<class C> bool any(const C& xs) { return std::any_of(xs.begin(), xs.end(), nonzero); }

// Indented printing
unique_ptr<ostream> out;
int indent = 0;
const bool debug = false;

struct Indent : public Noncopyable {
  const int outer;
  Indent() : outer(indent) { indent += 2; }
  ~Indent() { indent = outer; }
};

void line() {
  if (debug) cout << endl;
  deref(out) << endl;
}
template<class T> void line(T&& x) {
  const string dent(indent, ' ');
  if (debug) cout << dent << x << endl;
  deref(out) << dent << x << endl;
}
template<class... Args> void line(const Args&... args) { line(format(args...)); }

struct Scope : public Noncopyable {
  unique_ptr<Indent> i;
  template<class... Args> Scope(const Args&... args) {
    line(format(args...) + " {");
    i.reset(new Indent);
  }
  ~Scope() { i.reset(); line("}"); }
};

struct Blank : public Noncopyable { ~Blank() { line(); }};

struct Header : public Noncopyable {
  Header(const string& path, const string& comment, const vector<string>& includes) {
    slow_assert(!out);
    out.reset(new ofstream(path));
    line("// " + comment);
    line("// Autogenerated by codelets.cc");
    line("#pragma once");
    line();
    for (const auto& i : includes)
      line("#include \"%s\"", i);
    line("namespace mandelbrot {");
    line();
  }
  ~Header() {
    line("}  // namespace mandelbrot");
    out.reset();
  }
};

struct Var {
  string name;

  explicit Var(const string& name) : name(name) {}
  Var(const Var& x) = default;
  bool operator==(const Var& x) const { return name == x.name; }
  friend ostream& operator<<(ostream& out, const Var& x) { return out << x.name; }
};

struct VarHash {
  auto operator()(const Var& x) const { return std::hash<string>()(x.name); }
};

// List of variables
vector<Var> vars(const string& prefix, const int n) {
  vector<Var> xs;
  for (int i = 0; i < n; i++)
    xs.push_back(Var(format("%s%d", prefix, i)));
  return xs;
}

struct Exp {
  string exp;
  int prec;  // Precedence, following https://en.wikipedia.org/wiki/Operators_in_C_and_C++#Operator_precedence
  vector<Exp> args;
  unordered_set<Var,VarHash> deps;

  Exp(const Var& x) : exp(x.name), prec(0), deps{x} {}

  Exp(const string& exp, const int prec, const vector<Exp>& args)
    : exp(exp), prec(prec), args(args) {
    for (const auto& e : args)
      for (const auto& v : e.deps)
        this->deps.insert(v);
  }

  // Wrap in parens
  Exp parens() const {
    Exp e = *this;
    if (!(startswith(exp, "(") && endswith(exp, ")"))) {
      e.exp = format("(%s)", e.exp);
      e.prec = 2;
    }
    return e;
  }

  // Ensure we have precedence <= prec
  Exp ensure(const int prec) const { return this->prec <= prec ? *this : parens(); }

  Exp operator-() const {
    const auto e = ensure(3);
    return Exp(format("-%s", e.exp), 3, {e});
  }

  friend ostream& operator<<(ostream& out, const Exp& e) { return out << e.exp; }
};

vector<Exp> exps(const vector<Var>& x) { return vector<Exp>(x.begin(), x.end()); }

#define LEFT_TO_RIGHT_OP(op, prec) \
  Exp operator op(const Exp& x, const Exp& y) { \
    return Exp{format("%s " #op " %s", x.ensure(prec), y.ensure(prec-1)), prec, {x, y}}; \
  }
LEFT_TO_RIGHT_OP(+, 6)
LEFT_TO_RIGHT_OP(-, 6)
LEFT_TO_RIGHT_OP(*, 5)
#undef LEFT_TO_RIGHT_OP

Exp fma(const Exp& x, const Exp& y, const Exp& s) {
  return Exp{format("__builtin_fma(%s, %s, %s)", x, y, s), 2, {x, y, s}};
}

// const type name = exp;
struct Stmt {
  Var var;
  string type;
  Exp exp;

  friend ostream& operator<<(ostream& out, const Stmt& s) {
    return out << "const " << s.type << ' ' << s.var << " = " << s.exp << ';';
  }
};

// A basic block
struct Block {
  vector<Stmt> stmts;
  unordered_set<string> names;

  Block(const vector<Var>& inputs) {
    for (const auto& x : inputs)
      names.insert(x.name);
  }

  Var fresh(const string& prefix) {
    for (int n = 0;; n++) {
      const auto s = n ? format("%s%d", prefix, n) : prefix;
      if (names.insert(s).second)
        return Var(s);
    }
  }

  Exp add(const string& prefix, const string& type, Exp&& exp) {
    const auto name = fresh(prefix);
    stmts.push_back(Stmt{name, type, move(exp)});
    return name;
  }
  Exp add(const string& prefix, Exp&& exp) { return add(prefix, "auto", move(exp)); }

  void strip(const vector<Exp>& live_exps) {
    unordered_set<Var, VarHash> live;
    for (const auto& e : live_exps)
      for (const auto& x : e.deps)
        live.insert(x);
    for (int i = int(stmts.size()) - 1; i >= 0; i--) {
      const auto s = stmts.begin() + i;
      if (contains(live, s->var))
        for (const auto& d : s->exp.deps)
          live.insert(d);
      else
        stmts.erase(s);
    }
  }

  void lines() const {
    for (const auto& s : stmts)
      line(s);
  }
};

struct ExpansionBlock : public Block {
  using Block::Block;

  // For each i, ulp(x[i]) >= |x[i+1]|
  typedef tuple<Exp,Exp> E2;
  typedef tuple<optional<Exp>,optional<Exp>> Ez2;
  typedef vector<Exp> Es;

  Es neg(const Es& x) {
    Es nx;
    for (const auto& a : x)
      nx.push_back(-a);
    return nx;
  }

  // Turn a + b with no overlap properties into an expansion.
  // This is Theorem 7 in Shewchuck.
  E2 two_sum(const Exp& a, const Exp& b) {
    const auto x = add("x", a + b);
    const auto v = add("v", x - a);
    const auto y = add("y", (a - (x - v)).parens() + (b - v));
    return {x, y};
  }

  Ez2 two_sum(const optional<Exp>& a, const optional<Exp>& b) {
    if (!a) return {b, nullopt};
    if (!b) return {a, nullopt};
    return Ez2(two_sum(*a, *b));
  }

  // a * b as an expansion, using fused multiply-add
  E2 two_prod(const Exp& a, const Exp& b) {
    const auto x = add("x", a * b);
    const auto y = add("y", fma(a, b, -x));
    return {x, y};
  }

  // Sum two expansions, producing an expansion of the same size.
  // This is Figure 2 of Collange et al.
  Es collange_add(const Es& x, const Es& y) {
    const int n = int(x.size());
    slow_assert(n == int(y.size()));
    deque<Exp> rest;
    for (int i = 0; i < n; i++) {
      rest.push_back(x[i]);
      rest.push_back(y[i]);
    }
    vector<Exp> state;
    while (rest.size()) {
      auto t = pop_front(rest);
      for (int i = 0; i < int(state.size()); i++)
        tie(state[i], t) = two_sum(t, state[i]);
      state.push_back(t);
    }
    slow_assert(int(state.size()) == 2*n);
    return {state.begin(), state.begin() + n};
  }

  // Sum via negation + symbolic add
  Es collange_sub(const Es& x, const Es& y) {
    return collange_add(x, neg(y));
  }

  // Multiply two expansions, producing an expansion of the same size
  // This is Algorithm 2 of Collange et al.
  Es collange_mul(const Es& x, const Es& y) {
    const int n = int(x.size());
    slow_assert(n == int(y.size()));
    vector<Exp> pi;
    deque<optional<Exp>> s(n);
    for (int i = 0; i < n; i++) {
      deque<optional<Exp>> e(n), ep(n);
      for (int j = 0; j < n; j++) {
        const auto [p, ej] = two_prod(x[j], y[i]);
        e[j] = ej;
        tie(s[j], ep[j]) = two_sum(s[j], p);
      }
      pi.push_back(deref(pop_front(s)));
      s.emplace_back();
      while (any(e)) {
        for (int j = 0; j < n; j++)
          tie(s[j], e[j]) = two_sum(s[j], e[j]);
        e.pop_back();
        e.emplace_front();
      }
      while (any(ep)) {
        for (int j = 0; j < n; j++)
          tie(s[j], ep[j]) = two_sum(s[j], ep[j]);
        ep.pop_back();
        ep.emplace_front();
      }
    }
    slow_assert(n == int(pi.size()));
    return pi;
  }
};

void expansion_arithmetic(const string& path) {
  Header h(path, "Expansion arithmetic codelets", {});
  for (const int n : {2, 3, 4}) {
    // Negation
    {
      Blank b;
      Scope fun("__host__ __device__ static inline Expansion<%d> operator-(const Expansion<%d> x)", n, n);
      vector<string> nx;
      for (int i = 0; i < n; i++)
        nx.push_back(format("-x.x[%d]", i));
      line("return Expansion<%d>(%s, nonoverlap);", n, join(nx));
    }

    // Add, subtract, multiply
    const auto binary = [n](const string& op, const auto& body) {
      Blank b;
      Scope fun("__host__ __device__ static inline Expansion<%d>\n"
                "operator%s(const Expansion<%d> x, const Expansion<%d> y)",
                n, op, n, n);
      line("#ifdef __clang__");
      line("#pragma clang fp reassociate(off)");
      line("#endif  // __clang__");
      const auto x = vars("x", n), y = vars("y", n);
      line("const auto [%s] = x.x;", join(x));
      line("const auto [%s] = y.x;", join(y));
      ExpansionBlock B(concat(vector<Var>{Var("x"), Var("y")}, x, y));
      const auto s = body(B, exps(x), exps(y));
      B.strip(s);
      B.lines();
      line("return Expansion<%d>(%s, nonoverlap);", n, join(s));
    };
    binary("+", [](auto& B, const auto& x, const auto& y) { return B.collange_add(x, y); });
    binary("-", [](auto& B, const auto& x, const auto& y) { return B.collange_sub(x, y); });
    binary("*", [](auto& B, const auto& x, const auto& y) { return B.collange_mul(x, y); });
  }
}

}  // namespace mandelbrot
using namespace mandelbrot;

int main(int argc, char** argv) {
  try {
    const vector<string> paths(argv + 1, argv + argc);
    for (const auto& path : paths) {
      if (endswith(path, "gen-expansion.h"))
        expansion_arithmetic(path);
      else
        die("Unmatched path '%s'", path);
    }
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }
}
