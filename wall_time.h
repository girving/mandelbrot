// Wall clock time
#pragma once

#include <sys/time.h>
#include <cstdint>
namespace mandelbrot {

struct wall_time_t {
  int64_t us; // Microseconds since January 1, 1970, or relative microseconds

  wall_time_t() : us(0) {}
  explicit wall_time_t(int64_t us) : us(us) {}

  double seconds() const { return 1e-6*us; }
  double milliseconds() const { return 1e-3*us; }

  wall_time_t& operator+=(wall_time_t t) { us += t.us; return *this; }
  wall_time_t& operator-=(wall_time_t t) { us -= t.us; return *this; }
  wall_time_t operator-(wall_time_t t) const { return wall_time_t(us-t.us); }
};

static inline wall_time_t wall_time() {
  timeval tv;
  gettimeofday(&tv,0);
  return wall_time_t(tv.tv_sec*1000000+tv.tv_usec);
}

}  // namespace mandelbrot
