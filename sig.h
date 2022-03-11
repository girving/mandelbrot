// Modular arithmetic signatures for Schwartz-Zippel testing
#pragma once

#include "complex.h"
#include "debug.h"
#include "span.h"
#include <array>
#include <cstdint>
#include <iostream>
#include <optional>
namespace mandelbrot {

using std::array;
using std::optional;
using std::ostream;

// Arithmetic modulo a prime
struct Fp {
  typedef uint64_t T;
  T p;

  T neg(const T a) const { return a ? p - a : 0; }
  T add(const T a, const T b) const { return (__uint128_t(a) + b) % p; }
  T sub(const T a, const T b) const { return (__uint128_t(a) + p - b) % p; }
  T mul(const T a, const T b) const { return (__uint128_t(a) * b) % p; }

  T inv(const T a) const;
  T sqrt(const T a) const;
  bool has_sqrt(const T a) const;
};

// We compute in the Cartesian product of these fields
extern const array<Fp,4> fields;

// Signatures for expressions (values modulo a variety of small primes) for Schwartz–Zippel testing
struct Sig {
  static constexpr int n = fields.size();
  uint64_t x[n];

  Sig() : x{0} {}

  Sig(const int a) {
    for (int i = 0; i < n; i++)
      x[i] = a >= 0 ? a : a + fields[i].p;
  }

  bool operator==(const Sig s) const {
    for (int i = 0; i < n; i++)
      if (x[i] != s.x[i]) return false;
    return true;
  }

  Sig operator-() const;
  Sig operator+(const Sig s) const;
  Sig operator-(const Sig s) const;
  Sig operator*(const Sig s) const;

  Sig operator<<(const int k) const { return *this * Sig(1 << k); }
};

ostream& operator<<(ostream& out, const Sig s);
Sig inv(const Sig s);
Sig sqrt(const Sig s);
Sig random_sig();

Sig pow(const Sig s, const int n);
Complex<Sig> pow(const Complex<Sig> s, const int n);

// Detect signatures of small integers
optional<int> unsmall(const Sig s);

struct SigHash {
  // The first value is random enough for a hash
  auto operator()(const Sig& s) const { return std::hash<uint64_t>()(s.x[0]); }
};

// Deterministic pseudorandom signatures for taking rounding error into account
Sig arbitrary(const char* f, span<const Sig> ss);
static inline Sig arbitrary(const char* f, const Sig s) { return arbitrary(f, span<const Sig>(&s, 1)); }
static inline Sig arbitrary(span<const Sig> ss) { return arbitrary(0, ss); }
static inline Sig arbitrary(const Sig s0, const Sig s1) { const Sig ss[] = {s0, s1}; return arbitrary(0, ss); }

}  // namespace mandelbrot
