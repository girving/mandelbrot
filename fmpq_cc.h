// C++ interface to fmpq_t
#pragma once

#include <flint/fmpq.h>
#include <string>
namespace mandelbrot {

using std::string;

struct Fmpq {
  fmpq_t x;

  Fmpq() { fmpq_init(x); }
  Fmpq(const Fmpq& a) = delete;
  explicit Fmpq(const fmpq_t a) { fmpq_init(x); fmpq_set(x, a); }
  Fmpq(Fmpq&& a) { *x = *a.x; fmpq_init(a.x); }
  ~Fmpq() { fmpq_clear(x); }

  // Assignment
  Fmpq& operator=(const slong a) { fmpq_set_si(x, a, 1); return *this; }
  Fmpq& operator=(const Fmpq& a) { fmpq_set(x, a.x); return *this; }
  Fmpq& operator=(Fmpq&& a) { fmpq_swap(x, a.x); fmpq_zero(a.x); return *this; }

  // Implicit converson
  operator fmpq*() { return x; }
  operator const fmpq*() const { return x; }

  // Printing
  friend std::ostream& operator<<(std::ostream& out, const Fmpq& a);
};

}  // namespace mandelbrot
