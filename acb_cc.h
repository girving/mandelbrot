// C++ interface to acb_t
#pragma once

#include <flint/acb.h>
#include <string>
namespace mandelbrot {

using std::string;

struct Acb {
  acb_t x;

  Acb() { acb_init(x); }
  Acb(const Acb& a) = delete;
  explicit Acb(const acb_t a) { acb_init(x); acb_set(x, a); }
  Acb(Acb&& a) { *x = *a.x; acb_init(a.x); }
  ~Acb() { acb_clear(x); }

  // Assignment
  Acb& operator=(const slong a) { acb_set_si(x, a); return *this; }
  Acb& operator=(const Acb& a) { acb_set(x, a.x); return *this; }
  Acb& operator=(Acb&& a) { acb_swap(x, a.x); acb_zero(a.x); return *this; }

  // Implicit converson
  operator acb_ptr() { return x; }
  operator acb_srcptr() const { return x; }

  // Printing
  friend std::ostream& operator<<(std::ostream& out, const Acb& a);
  string safe() const;
};

}  // namespace mandelbrot
