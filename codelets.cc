// Code generation

#define CODELETS 1

#include "arith.h"
#include "debug.h"
#include "format.h"
#include "noncopyable.h"
#include "print.h"
#include "series.h"
#include <array>
#include <deque>
#include <fstream>
#include <memory>
#include <optional>
#include <random>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>
namespace mandelbrot {

using std::any_of;
using std::array;
using std::common_type_t;
using std::cout;
using std::declval;
using std::deque;
using std::endl;
using std::independent_bits_engine;
using std::make_pair;
using std::make_tuple;
using std::mt19937;
using std::nullopt;
using std::ofstream;
using std::optional;
using std::ostream;
using std::remove_cvref_t;
using std::string_view;
using std::tie;
using std::tuple;
using std::unique_ptr;
using std::unordered_map;
using std::unordered_set;
using std::vector;

struct Unit {};

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
bool endswith(const string_view s, const string_view end) {
  if (s.size() < end.size()) return false;
  return s.substr(s.size() - end.size()) == end;
}

// Cache a function
template<class A,class F> struct Cache {
  typedef decltype(declval<F>()(declval<A>())) B;
  const F f;
  unordered_map<A,B> cache;
  B operator()(const A& x) {
    const auto it = cache.find(x);
    if (it != cache.end()) return it->second;
    const auto y = f(x);
    cache.insert(make_pair(x, y));
    return y;
  }
};
template<class A,class F> auto cache(F&& f) { return Cache<A,remove_cvref_t<F>>{f, {}}; }

// Signatures for polynomials (values modulo a variety of small primes) for Schwartz–Zippel testing
static const uint64_t primes[] = {
  9223372036854775837ul,
  9223372036854775907ul,
  9223372036854775931ul,
  9223372036854775939ul,
};
struct Sig {
  static constexpr int n = sizeof(primes) / sizeof(uint64_t);
  uint64_t x[n];

  Sig() : x{0} {}

  explicit Sig(const int a) {
    for (int i = 0; i < n; i++)
      x[i] = a >= 0 ? a : a + primes[i];
  }

  bool operator==(const Sig s) const {
    for (int i = 0; i < n; i++)
      if (x[i] != s.x[i]) return false;
    return true;
  }

  friend ostream& operator<<(ostream& out, const Sig s) {
    const bool first = true;
    if (first) out << s.x[0];
    else out << span<const uint64_t>(s.x);
    return out;
  }

  Sig operator-() const {
    Sig r;
    for (int i = 0; i < n; i++) r.x[i] = x[i] ? primes[i] - x[i] : 0;
    return r;
  }

  Sig operator+(const Sig s) const {
    Sig r;
    for (int i = 0; i < n; i++) {
      const auto t = __uint128_t(x[i]) + s.x[i];
      r.x[i] = t < primes[i] ? t : t - primes[i];
    }
    return r;
  }

  Sig operator-(const Sig s) const { return *this + (-s); }

  Sig operator*(const Sig s) const {
    Sig r;
    for (int i = 0; i < n; i++)
      r.x[i] = (__uint128_t(x[i]) * s.x[i]) % primes[i];
    return r;
  }
};

// Find s,t s.t. sx + tp = 1, so that a^{-1} = s (mod p)
// From https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Pseudocode.
uint64_t inv_mod(const uint64_t x, const uint64_t p) {
  // Work in signed 128 bits for super laziness.
  __int128_t rp = x, r = p;
  __int128_t sp = 1, s = 0;
  __int128_t tp = 0, t = 1;
  while (r) {
    const auto q = rp / r;
    tie(rp, r) = make_tuple(r, rp - q*r);
    tie(sp, s) = make_tuple(s, sp - q*s);
    tie(tp, t) = make_tuple(t, tp - q*t);
  }
  slow_assert(sp*x + tp*p == 1);
  // 1 = sx + tp = (s+p)x + (t-x)p
  if (sp < 0) sp += p;
  slow_assert(sp >= 0);
  return sp;
}

Sig inv(const Sig s) {
  Sig r;
  for (int i = 0; i < Sig::n; i++)
    r.x[i] = inv_mod(s.x[i], primes[i]);
  slow_assert(r * s == Sig(1));
  return r;
}

Sig random_sig() {
  static independent_bits_engine<mt19937,64,uint64_t> mt(7);
  Sig s;
  for (int i = 0; i < Sig::n; i++) {
    const auto bits = mt();
    static_assert(is_same_v<decltype(bits),const uint64_t>);
    s.x[i] = bits;
  }
  return s;
}

struct SigHash {
  // The first value is random enough for a hash
  auto operator()(const Sig& s) const { return std::hash<uint64_t>()(s.x[0]); }
};

// Detect small signatures
static unordered_map<Sig,int,SigHash> make_smalls() {
  unordered_map<Sig,int,SigHash> smalls;
  for (int i = 0; i <= 32; i++) {
    smalls[Sig(i)] = i;
    smalls[-Sig(i)] = -i;
  }
  return smalls;
}
optional<int> unsmall(const Sig s) {
  static const auto smalls = make_smalls();
  const auto it = smalls.find(s);
  return it != smalls.end() ? optional<int>(it->second) : nullopt;
}

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
  const string close;
  unique_ptr<Indent> i;
  template<class... Args> Scope(const string& close, const Args&... args)
    : close(close) {
    line(format(args...));
    i.reset(new Indent);
  }
  ~Scope() { i.reset(); line(close); }
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
  Sig sig;

  Var(const string& name, const Sig sig) : name(name), sig(sig) {}
  Var(const Var& x) = default;
  bool operator==(const Var& x) const { return name == x.name; }
  friend ostream& operator<<(ostream& out, const Var& x) { return out << x.name; }
};

struct VarHash {
  auto operator()(const Var& x) const { return std::hash<string>()(x.name); }
};

// Variables with random signatures
Var input_var(const string& name) { return Var(name, random_sig()); }
vector<Var> input_vars(const string& prefix, const int n) {
  vector<Var> xs;
  for (int i = 0; i < n; i++)
    xs.push_back(input_var(format("%s%d", prefix, i)));
  return xs;
}

struct Exp {
private:
  struct Data {
    string exp;
    bool atom;
    int prec;  // Precedence, following https://en.wikipedia.org/wiki/Operators_in_C_and_C++#Operator_precedence
    Sig sig;
    vector<Exp> args;
  };
  shared_ptr<const Data> e;
public:

  Exp() : Exp(Var("nonsense", random_sig())) {}
  Exp(const Var& x) : e(new Data{x.name, true, 0, x.sig, {}}) {}
  Exp(const int a) : e(new Data{format("%d", a), true, 0, Sig(a), {}}) {}
  Exp(const string& exp, const int prec, const Sig s, const vector<Exp>& args)
    : e(new Data{exp, false, prec, s, args}) {}
  Exp(const Exp& e) : e(e.e) {}
  Exp(Exp&& e) : e(e.e) {}  // Stay valid on move
  Exp& operator=(const Exp& f) { e = f.e; return *this; }

  bool zero() const { return e->exp == "0"; }
  bool one() const { return e->exp == "1"; }
  bool atom() const { return e->atom; }
  int prec() const { return e->prec; }
  const Sig& sig() const { return e->sig; }
  const string& exp() const { return e->exp; }
  bool negation() const { return prec() == 3 && startswith(exp(), "-"); }
  bool parens() const { return prec() == 2 && startswith(exp(), "(") && endswith(exp(), ")"); }
  const Exp& arg() const { slow_assert(e->args.size() == 1, exp()); return e->args[0]; }

  unordered_set<string> deps() const {
    unordered_set<string> deps;
    for (const auto& a : e->args)
      for (const auto& v : a.deps())
        deps.insert(v);
    if (atom())
      deps.insert(exp());
    return deps;
  }

  // Wrap in parens
  Exp add_parens() const {
    if (atom() || parens()) return *this;
    return Exp(format("(%s)", exp()), 2, sig(), {*this});
  }

  // Ensure we have precedence <= prec
  Exp ensure(const int prec) const { return this->prec() <= prec ? *this : add_parens(); }

  Exp operator-() const {
    if (zero()) return *this;
    if (negation()) {
      if (atom()) {
        const int i = atoi(exp().c_str());
        slow_assert(i < 0);
        return Exp(-i);
      }
      return arg();
    }
    const auto e = ensure(3);
    return Exp(format("-%s", e.exp()), 3, -sig(), {e});
  }

  friend ostream& operator<<(ostream& out, const Exp& e) { return out << e.exp(); }
};

template<> struct IsIntervalT<Exp> { static constexpr bool value = false; };

vector<Exp> exps(const vector<Var>& x) { return vector<Exp>(x.begin(), x.end()); }

Exp operator+(const Exp& x, const Exp& y);
Exp operator-(const Exp& x, const Exp& y);

Exp operator+(const Exp& x, const Exp& y) {
  if (x.zero()) return y;
  if (y.zero()) return x;
  if (x.negation()) return y - (-x);
  if (y.negation()) return x - (-y);
  return Exp{format("%s + %s", x.ensure(6), y.ensure(5)), 6, x.sig() + y.sig(), {x, y}};
}

Exp operator-(const Exp& x, const Exp& y) {
  if (x.zero()) return -y;
  if (y.zero()) return x;
  if (x.negation()) return -((-x) + y);
  if (y.negation()) return x + (-y);
  return Exp{format("%s - %s", x.ensure(6), y.ensure(5)), 6, x.sig() - y.sig(), {x, y}};
}

Exp operator*(const Exp& x, const Exp& y) {
  if (x.zero() || y.zero()) return 0;
  if (x.negation()) return -((-x) * y);
  if (y.negation()) return -(x * (-y));
  if (x.one()) return y;
  if (y.one()) return x;
  const auto sig = x.sig() * y.sig();
  if (x.exp() == y.exp()) return Exp{format("sqr(%s)", x), 2, sig, {x}};
  return Exp{format("%s * %s", x.ensure(5), y.ensure(4)), 5, sig, {x, y}};
}

Exp operator/(const Exp& x, const int a) {
  slow_assert(a);
  if (a == 1) return x;
  if (a == -1) return -x;
  return Exp{format("%s / %d", x.ensure(5), a), 5, x.sig() * inv(Sig(a)), {x}};
}

Exp fma(const Exp& x, const Exp& y, const Exp& s) {
  return Exp{format("__builtin_fma(%s, %s, %s)", x, y, s), 2, random_sig(), {x, y, s}};
}

template<class Op> Exp reduce(const Op& op, const vector<Exp>& xs) {
  const size_t n = xs.size();
  slow_assert(n);
  if (n == 1) return xs[0];
  const auto p = xs.begin();
  return op(reduce(op, vector<Exp>(p, p + n/2)),
            reduce(op, vector<Exp>(p + n/2, p + n)));
}
Exp sum(const vector<Exp>& xs) { return reduce(std::plus<void>(), xs); }

// const type name = exp;
struct Stmt {
  optional<Var> var;
  string type;
  Exp exp;

  Stmt(const Var& var, const string& type, const Exp& exp) : var(var), type(type), exp(exp) {}
  explicit Stmt(const string& s) : type("void"), exp(s, 100, Sig(), {}) {}
  Stmt() : Stmt("") {}

  friend ostream& operator<<(ostream& out, const Stmt& s) {
    const bool sigs = true;
    if (s.var) {
      out << "const " << s.type << ' ' << *s.var << " = " << s.exp << ';';
      if (sigs) out << "  // sig = " << s.exp.sig();
    } else
      if (!s.exp.exp().empty())
        out << s.exp << ';';
    return out;
  }
};

// A basic block
struct Block : public Noncopyable {
private:
  static Block* active_;
public:
  const bool sz_cse;  // Whether to use Schwartz-Zippel for CSE
  vector<Stmt> stmts;
  unordered_map<string,Var> expressions;  // For CSE on string expressions
  unordered_map<Sig,Var,SigHash> signatures;  // For CSE via Schwartz–Zippel
  unordered_set<string> names;

  Block(const vector<Var>& inputs, const bool sz_cse)
    : sz_cse(sz_cse) {
    slow_assert(!active_);
    active_ = this;
    for (const auto& x : inputs)
      names.insert(x.name);
  }
  ~Block() { active_ = 0; }

  static Block& active() { slow_assert(active_); return *active_; }

  string fresh(const string& prefix) {
    for (int n = -1;; n++) {
      const auto s = n < 0 ? prefix : format("%s%d", prefix, n);
      if (names.insert(s).second)
        return s;
    }
  }

  Exp add(const string& prefix, const string& type, const Exp& exp) {
    if (exp.atom()) return exp;
    if (exp.negation()) return -add(prefix, type, exp.arg());
    if (exp.parens()) return add(prefix, type, exp.arg());
    // CSE
    if (sz_cse) {
      const auto small = unsmall(exp.sig());
      if (small) return Exp(*small);
      const auto it = signatures.find(exp.sig());
      if (it != signatures.end()) return it->second;
    } else {
      const auto it = expressions.find(exp.exp());
      if (it != expressions.end()) return it->second;
    }
    // Fine, compute it
    const Var var(fresh(prefix), exp.sig());
    stmts.push_back(Stmt(var, type, exp));
    // Remember
    signatures.insert(make_pair(exp.sig(), var));
    expressions.insert(make_pair(exp.exp(), var));
    return var;
  }
  Exp add(const string& prefix, const Exp& exp) { return add(prefix, "auto", exp); }

  Exp sum(const string& prefix, const vector<Exp>& exps) {
    return reduce([this,&prefix](const auto& x, const auto& y) { return add(prefix, x + y); }, exps);
  }

  void strip(span<const Exp> roots) {
    unordered_set<string> live;
    for (const auto& e : roots)
      for (const auto& x : e.deps())
        live.insert(x);
    for (int i = int(stmts.size()) - 1; i >= 0; i--) {
      const auto s = stmts.begin() + i;
      slow_assert(s->var);
      if (contains(live, s->var->name))
        for (const auto& d : s->exp.deps())
          live.insert(d);
      else
        stmts.erase(s);
    }
  }

  // Reorder statements so that each exp can be returned in order
  void multi_strip(span<const Exp> roots, span<const vector<Stmt>> returns) {
    // Sort statements into classed based on the earliest root that depends on them
    const int rn = roots.size();
    slow_assert(returns.empty() || rn == int(returns.size()));
    unordered_map<string,int> earliest;
    const auto depend = [&earliest](const string& s, int r) {
      const auto it = earliest.find(s);
      r = it != earliest.end() ? min(r, it->second) : r;
      earliest[s] = r;
    };
    for (int r = 0; r < rn; r++)
      for (const auto& x : roots[r].deps())
        depend(x, r);
    vector<vector<Stmt>> reorder(rn);
    for (int i = int(stmts.size()) - 1; i >= 0; i--) {
      const auto& s = stmts[i];
      slow_assert(s.var);
      const auto it = earliest.find(s.var->name);
      if (it != earliest.end()) {
        const int r = it->second;
        for (const auto& d : s.exp.deps())
          depend(d, r);
        reorder[r].push_back(s);
      }
    }

    // Output reordered statements
    stmts.clear();
    for (int r = 0; r < rn; r++) {
      for (int i = int(reorder[r].size()) - 1; i >= 0; i--)
        stmts.push_back(reorder[r][i]);
      if (returns.size())
        for (const auto& s : returns[r])
          stmts.push_back(s);
    }
  }

  void lines() const {
    for (const auto& s : stmts)
      line(s);
  }
};

Block* Block::active_ = 0;

struct ExpansionBlock : public Block {
  // No Schwartz-Zippel for expansion arithmetic!
  static constexpr bool sz_cse = false;
  ExpansionBlock(const vector<Var>& inputs) : Block(inputs, sz_cse) {}

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
    const auto y = add("y", (a - (x - v)) + (b - v));
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
      Scope fun("}", "__host__ __device__ static inline Expansion<%d> operator-(const Expansion<%d> x) {", n, n);
      vector<string> nx;
      for (int i = 0; i < n; i++)
        nx.push_back(format("-x.x[%d]", i));
      line("return Expansion<%d>(%s, nonoverlap);", n, join(nx));
    }

    // Add, subtract, multiply
    const auto binary = [n](const string& op, const auto& body) {
      Blank b;
      Scope fun("}", "__host__ __device__ static inline Expansion<%d>\n"
                "operator%s(const Expansion<%d> x, const Expansion<%d> y) {",
                n, op, n, n);
      line("#ifdef __clang__");
      line("#pragma clang fp reassociate(off)");
      line("#endif  // __clang__");
      const auto x = input_vars("x", n), y = input_vars("y", n);
      line("const auto [%s] = x.x;", join(x));
      line("const auto [%s] = y.x;", join(y));
      ExpansionBlock B(concat(vector<Var>{input_var("x"), input_var("y")}, x, y));
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

// Zero extension
auto extend(const string& x) {
  return cache<int>([x](const int i) {
    const auto v = format("%s%d", x, i);
    line("const auto %s = %d < n%s ? %s[%d] : S(0);", v, i, x, x, i);
    return Exp(input_var(v));
  });
};

// Early exit
const auto early = [](const string& n) {
  return cache<int>([n](const int i) { line("if (%s <= %d) return;", n, i); return Unit(); });
};

// Base cases for series multiplication and squaring
void mul_bases(const string& path) {
  Header h(path, "Multiplication and squaring base cases", {"loops.h"});
  const int n = 8;
  line("static constexpr int mul_base_n = %d;", n);
  line("static constexpr int sqr_base_n = %d;\n", n);

  // Multiplication
  {
    Blank b;
    Scope f(")", "DEF_SERIAL(mul_base, (S* z, const int nz, const S* x, const int nx, const S* y, const int ny),");
    auto exit = early("nz");
    exit(0);
    auto xi = extend("x");
    auto yi = extend("y");
    // Cache inputs so that aliasing works
    for (int i = 0; i < n; i++) xi(i);
    for (int i = 0; i < n; i++) yi(i);
    line();
    for (int i = 0; i < n; i++) {
      exit(i);
      vector<Exp> ts;
      for (int j = 0; j <= i; j++)
        ts.push_back(xi(j) * yi(i - j));
      line("z[%d] = %s;", i, sum(ts));
    }
  };

  // Squaring
  {
    Blank b;
    Scope f(")", "DEF_SERIAL(sqr_base, (S* y, const int ny, const S* x, const int nx),");
    auto exit = early("ny");
    exit(0);
    auto xi = extend("x");
    for (int i = 0; i < n; i++) xi(i); // Cache inputs so that aliasing works
    line();
    for (int i = 0; i < n; i++) {
      exit(i);
      vector<Exp> ts;
      for (int j = 0; j <= i; j++)
        ts.push_back(xi(j) * xi(i - j));
      line("y[%d] = %s;", i, sum(ts));
    }
  };
}

// Override series arithmetic
void add_scalar(Series<Exp>& x, const Exp a) {
  slow_assert(x.nonzero());
  auto& B = Block::active();
  x[0] = B.add("s", x[0] + a);
}
void high_addsub(Series<Exp>& y, const int sign, const int64_t s, SeriesView<const Exp> x) {
  auto& B = Block::active();
  const auto ynz = y.nonzero(), xnz = x.nonzero();
  const auto nk = min(y.known(), x.known() + s);
  const auto nz = min(nk, max(ynz, xnz ? xnz + s : 0));
  slow_assert(abs(sign) == 1 && nz <= y.limit());
  const auto x_ = x.copy(xnz);  // Watch for aliasing
  y.set_counts(nk, nz);
  for (int i = 0; i < nz; i++) {
    const auto yi = i < ynz ? y[i] : Exp(0);
    auto xi = uint32_t(i-s) < uint32_t(xnz) ? x_[i-s] : Exp(0);
    if (sign < 0) xi = -xi;
    y[i] = B.add("s", yi + xi);
  }
}
void fft_mul(span<Exp> z, span<const Exp> x, span<const Exp> y) {
  // Actually, we use naive multiplication
  auto& B = Block::active();
  const int nz = z.size(), nx = x.size(), ny = y.size();
  for (int i = 0; i < nz; i++) {
    vector<Exp> ts;
    for (int j = 0; j <= i; j++)
      if (j < nx && i - j < ny)
        ts.push_back(B.add("t", x[j] * y[i - j]));
    z[i] = B.sum("z", ts);
  }
}

Exp inv(const Exp& x) { return Block::active().add("r", Exp{format("inv(%s)", x), 2, inv(x.sig()), {x}}); }

// Base cases for series functions
void series_bases(const string& path) {
  Header h(path, "Series function base cases", {"loops.h"});
  const int n = 4;
  const bool sz_cse = true;

  // Unary
  const auto unary = [n](const string& name, const int leading, auto&& set) {
    {
      Blank b;
      Scope f(")", "DEF_SERIAL(%s_base_serial, (S* y, const int ny, const S* x, const int nx),", name);
      vector<Var> inputs{input_var("y"), input_var("ny"), input_var("x"), input_var("nx")};
      for (int i = 0; i < leading; i++) inputs.push_back(input_var(format("x%d", i)));  // Make names nicer
      Block B(inputs, sz_cse);
      Series<Exp> x(n), y(n);
      x.set_counts(n, n);
      for (int i = 0; i < n; i++)
        x[i] = i < leading ? Exp(0) : B.add("x", Exp(format("%d < nx ? x[%d] : S(0)", i, i), 16, random_sig(), {}));
      set(y, x.view());
      vector<vector<Stmt>> returns(n);
      for (int i = 0; i < n; i++) {
        returns[i].push_back(Stmt(format("y[%d] = %s", i, y[i])));
        if (i+1 < n) {
          returns[i].push_back(Stmt());
          returns[i].push_back(Stmt(format("if (ny == %d) return", i+1)));
        }
      }
      B.multi_strip(y, returns);
      B.lines();
    } {
      Blank b;
      Scope f("}", "template<class T> void %s_base(Series<T>& y, type_identity_t<SeriesView<const T>> x) {", name);
      line("const int ny = min(%d, int(x.known()));", n);
      line("y.set_counts(ny, ny);");
      line("if (ny) %s_base_serial(y.data(), ny, x.data(), x.nonzero());", name);
    }
  };
  #define UNARY(name, leading) unary(#name, leading, [=](auto& y, const auto x) { y = name(x); });
  UNARY(inv, 0)
  UNARY(exp, 1)
}

}  // namespace mandelbrot
using namespace mandelbrot;

int main(int argc, char** argv) {
  try {
    const vector<string> paths(argv + 1, argv + argc);
    for (const auto& path : paths) {
      if (endswith(path, "gen-expansion.h")) expansion_arithmetic(path);
      else if (endswith(path, "gen-mul-bases.h")) mul_bases(path);
      else if (endswith(path, "gen-series-bases.h")) series_bases(path);
      else die("Unmatched path '%s'", path);
    }
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }
}
