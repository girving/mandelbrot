// C++ interface to flint_rand_t
#pragma once

#include <flint/flint.h>
namespace mandelbrot {

struct Rand {
  flint_rand_t x;

  Rand() { flint_randinit(x); }
  Rand(ulong seed0, ulong seed1) { flint_randinit(x); seed(seed0, seed1); }
  ~Rand() { flint_randclear(x); }

  // Noncopyable
  Rand(const Rand& a) = delete;
  Rand(Rand&& a) = delete;
  Rand& operator=(const Rand& a) = delete;
  Rand& operator=(Rand&& a) = delete;

  // Seeding
  void seed(ulong seed0, ulong seed1) { flint_randseed(x, seed0, seed1); }

  // Implicit converson
  operator flint_rand_s*() { return x; }
  operator const flint_rand_s*() const { return x; }
};

}  // namespace mandelbrot
