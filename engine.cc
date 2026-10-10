// Parallel execution of resumable orbits and simple loops, on CPU threads or the GPU

#include "engine.h"
#include <cstdlib>
#include <numeric>
namespace mandelbrot {

int64_t scramble_stride(const int64_t n) {
  slow_assert(0 <= n && n < (int64_t(1) << 31), "scramble_stride: n = %d too large", n);
  if (n <= 2) return 1;
  int64_t s = std::max<int64_t>(1, int64_t(0.6180339887498949 * double(n)));
  while (std::gcd(s, n) != 1) s++;
  return s;
}

int env_int(const char* name, const int fallback) {
  const char* s = getenv(name);
  return s ? atoi(s) : fallback;
}

int cpu_threads() {
  static const int n = std::max(1, env_int("MANDELBROT_THREADS", int(std::thread::hardware_concurrency())));
  return n;
}

}  // namespace mandelbrot
