// C++ interface to arb_t
#pragma once

#include <flint/arb.h>
#include <string>
namespace mandelbrot {

using std::string;

struct Arb {
  arb_t x;

  Arb() { arb_init(x); }
  Arb(const Arb& a) = delete;
  explicit Arb(const arb_t a) { arb_init(x); arb_set(x, a); }
  Arb(Arb&& a) { *x = *a.x; arb_init(a.x); }
  ~Arb() { arb_clear(x); }

  // Assignment
  void operator=(const slong a) { arb_set_si(x, a); }
  void operator=(const Arb& a) { arb_set(x, a.x); }
  void operator=(Arb&& a) { arb_swap(x, a.x); arb_zero(a.x); }

  // Implicit converson
  operator arb_ptr() { return x; }
  operator arb_srcptr() const { return x; }

  // Upper bound for |x|
  friend double bound(const Arb& x);

  // Printing
  friend std::ostream& operator<<(std::ostream& out, const Arb& a);
  string safe() const;
};

static inline Arb exact_arb(const double a) {
  Arb x;
  arb_set_d(x, a);
  return x;
}

}  // namespace mandelbrot
