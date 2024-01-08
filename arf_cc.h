// C++ interface to arf_t
#pragma once

#include <flint/arf.h>
#include <string>
namespace mandelbrot {

using std::string;

struct Arf {
  arf_t x;

  Arf() { arf_init(x); }
  Arf(const Arf& a) = delete;
  Arf(Arf&& a) { *x = *a.x; arf_init(a.x); }
  ~Arf() { arf_clear(x); }

  // Assignment
  Arf& operator=(const Arf& a) { arf_set(x, a.x); return *this; }
  Arf& operator=(Arf&& a) { arf_swap(x, a.x); arf_zero(a.x); return *this; }

  // Implicit converson
  operator arf_ptr() { return x; }
  operator arf_srcptr() const { return x; }

  // Printing
  string safe() const;  // Integer mantissa and exponent
  friend std::ostream& operator<<(std::ostream& out, const Arf& a);
};

}  // namespace mandelbrot
