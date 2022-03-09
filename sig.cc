// Modular arithmetic signatures for Schwartz-Zippel testing

#define OPENSSL_SUPPRESS_DEPRECATED  // Use SHA256_* without warnings

#include "sig.h"
#include "bit.h"
#include "print.h"
#include <openssl/sha.h>
#include <random>
#include <tuple>
#include <unordered_map>
namespace mandelbrot {

using std::independent_bits_engine;
using std::is_same_v;
using std::make_tuple;
using std::min;
using std::mt19937;
using std::nullopt;
using std::tie;
using std::unordered_map;

// a^n in F
template<class Field> static typename Field::T pow(const Field F, typename Field::T a, uint64_t n) {
  typename Field::T r = 1;
  while (n) {
    if (n & 1) r = F.mul(r, a);
    a = F.mul(a, a);
    n >>= 1;
  }
  return r;
}

// Fast primality test for 64-bit integers.  Modified from
// 1. Forisek and Jancina, Fast primality testing for integers that fit into a machine word.
//    http://ceur-ws.org/Vol-1326/020-Forisek.pdf
// 2. https://miller-rabin.appspot.com

// Check whether n is a strong pseudoprime to base a
static bool is_sprp(const uint64_t n, const uint64_t a) {
  if (n % a == 0) return n == a;
  const int s = countr_zero(n - 1);
  const auto d = (n - 1) >> s;
  auto x = pow(Fp{n}, a, d);
  if (x == 1) return true;
  for (int r = 0; r < s; r++) {
    if (x == n-1) return true;
    x = Fp{n}.mul(x, x);
  }
  return false;
}
static bool is_prime(const uint64_t n) {
  for (const uint32_t p : {2, 3, 5, 7})
    if (n % p == 0)
      return n == p;
  for (const uint64_t a : {2, 325, 9375, 28178, 450775, 9780504, 1795265022})
    if (!is_sprp(n, a))
      return false;
  return true;
}

// Find s,t s.t. sx + tp = 1, so that a^{-1} = s (mod p)
// From https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Pseudocode.
uint64_t Fp::inv(const uint64_t x) const {
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
  slow_assert(mul(sp, x) == 1);
  return sp;
}

static uint64_t random(const Fp F) {
  static independent_bits_engine<mt19937,64,uint64_t> mt(7);
  const auto bits = mt();
  static_assert(is_same_v<decltype(bits),const uint64_t>);
  return bits % F.p;
}

bool Fp::has_sqrt(const uint64_t a) const {
  if (!a) return true;
  const auto s = pow(*this, a, (p-1)>>1);
  slow_assert(s == 1 || s == neg(1));
  return s == 1;
}

// Fp + Fp sqrt(b)
struct Fp2 {
  Fp F;
  uint64_t b;

  // u + v sqrt(b)
  struct T {
    uint64_t u, v;
    T(unsigned u) : u(u), v(0) {}
    T(uint64_t u, uint64_t v) : u(u), v(v) {}
  };

  T mul(const T x, const T y) const {
    // x y = (xu + xv sqrt(b)) (yu + yv sqrt(b))
    //     = (xu yu + xv yv b) + (xu yv + xv yu) sqrt(b)
    return T(F.add(F.mul(x.u, y.u), F.mul(F.mul(x.v, y.v), b)),
             F.add(F.mul(x.u, y.v), F.mul(x.v, y.u)));
  };
};

// Modular square root via https://en.wikipedia.org/wiki/Cipolla's_algorithm
uint64_t Fp::sqrt(const uint64_t n) const {
  if (!n) return n;
  const auto& F = *this;

  // Find a s.t. a^2 - n is not a square
  uint64_t a, b;
  for (;;) {
    a = random(F);
    b = sub(mul(a, a), n);
    if (!has_sqrt(b)) break;
  }

  // Compute x = (a + sqrt(b))^((p+1)/2)
  const Fp2 F2{F, b};
  const auto x = pow(F2, Fp2::T{a, 1}, (F.p+1)/2);
  slow_assert(x.v == 0);
  const auto r = min(x.u, F.neg(x.u));
  slow_assert(F.mul(r, r) == n);
  return r;
}

const array<Fp,4> fields = []() {
  const int two_roots = 4;
  array<Fp,4> fields;
  int count = 0;
  for (uint64_t p = uint64_t(0) - 1;; p--) {
    if (p % 8 != 7 || !is_prime(p)) continue;
    const Fp F{p};
    if (F.has_sqrt(F.neg(1))) continue;  // Make sure i ∉ F so that Complex<F> works
    uint64_t r = 2;
    for (int i = 0; i < two_roots; i++) {
      if (!F.has_sqrt(r)) goto skip;
      r = F.sqrt(r);
    }
    fields[count++] = F;
    if (count == fields.size()) break;
    skip:;
  }
  return fields;
}();

ostream& operator<<(ostream& out, const Sig s) {
  const bool first = true;
  if (first) out << s.x[0];
  else out << span<const uint64_t>(s.x);
  return out;
}

#define OP(exp) \
  Sig r; \
  for (int i = 0; i < Sig::n; i++) { const Fp F = fields[i]; r.x[i] = (exp); } \
  return r;

Sig Sig::operator-() const { OP(F.neg(x[i])) }
Sig Sig::operator+(const Sig s) const { OP(F.add(x[i], s.x[i])) }
Sig Sig::operator-(const Sig s) const { OP(F.sub(x[i], s.x[i])) }
Sig Sig::operator*(const Sig s) const { OP(F.mul(x[i], s.x[i])) }
Sig inv(const Sig s) { OP(F.inv(s.x[i])); }
Sig sqrt(const Sig s) { OP(F.sqrt(s.x[i])); }

Sig random_sig() {
  Sig s;
  for (int i = 0; i < Sig::n; i++)
    s.x[i] = random(fields[i]);
  return s;
}

optional<int> unsmall(const Sig s) {
  static const auto smalls = []() {
    unordered_map<Sig,int,SigHash> smalls;
    for (int i = 0; i <= 32; i++) {
      smalls[Sig(i)] = i;
      smalls[-Sig(i)] = -i;
    }
    return smalls;
  }();
  const auto it = smalls.find(s);
  return it != smalls.end() ? optional<int>(it->second) : nullopt;
}

Sig arbitrary(const char* f, span<const Sig> ss) {
  static_assert(sizeof(Sig) == SHA256_DIGEST_LENGTH);
  SHA256_CTX ctx;
  SHA256_Init(&ctx);
  if (f) SHA256_Update(&ctx, f, strlen(f));
  SHA256_Update(&ctx, reinterpret_cast<const unsigned char*>(ss.data()), sizeof(Sig)*ss.size());
  Sig r;
  SHA256_Final(reinterpret_cast<unsigned char*>(&r), &ctx);
  return r;
}

}  // namespace mandelbrot
