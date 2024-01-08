// C++ interface to mag_t
#pragma once

#include <flint/mag.h>
namespace mandelbrot {

struct Mag {
  mag_t x;

  Mag() { mag_init(x); }
  Mag(const Mag& a) = delete;
  Mag(Mag&& a) { *x = *a.x; mag_init(a.x); }
  ~Mag() { mag_clear(x); }

  // Assignment
  Mag& operator=(const Mag& a) { mag_set(x, a.x); return *this; }
  Mag& operator=(Mag&& a) { mag_swap(x, a.x); mag_zero(a.x); return *this; }

  // Implicit converson
  operator mag_ptr() { return x; }
  operator mag_srcptr() const { return x; }
};

}  // namespace mandelbrot
