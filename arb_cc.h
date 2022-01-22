// C++ interface to arb_t
#pragma once

#include <arb.h>
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
  Arb& operator=(const slong a) { arb_set_si(x, a); return *this; }
  Arb& operator=(const Arb& a) { arb_set(x, a.x); return *this; }
  Arb& operator=(Arb&& a) { arb_swap(x, a.x); arb_zero(a.x); return *this; }

  // Implicit converson
  operator arb_ptr() { return x; }
  operator arb_srcptr() const { return x; }

  // Printing
  friend std::ostream& operator<<(std::ostream& out, const Arb& a);
  string safe(const slong n) const { return arb_get_str(x, n, 0); }
};

}  // namespace mandelbrot
