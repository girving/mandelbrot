// Per-iteration cost of escape-time loops on the GPU (or CPU): a bare z^2 + c loop as the practical peak,
// then the same loop plus each check Orbit::run does, then Orbit::run itself.  All orbits are interior
// points near 0 whose orbits stay bounded, iterated for a fixed number of steps, so every variant does the
// same arithmetic and there is no divergence.  Parameters are c = μ/2 - μ^2/4 with μ = r e^{2πiθ}, θ the golden
// mean and 1 - r ≈ 1e-9: just inside the main cardioid next to a Siegel disk, where orbits neither escape nor
// converge (so neither the escape test nor Brent fires) for about 1e9 steps.

#include "engine.h"
#include "orbit.h"
#include "print.h"
#include <chrono>
#include <cmath>
namespace mandelbrot {
namespace {

const int64_t kSteps = 1 << 14;

// Variant 0: bare loop.  1: + escape test.  2: + atom-domain minimum.  3: + Brent check.  4: Orbit::run.
// 5: OrbitDE::run without Newton.  6: OrbitDE::run with Newton attempts from step 64.
template<int variant> __host__ __device__ double orbit(const double x, const double y) {
  if constexpr (variant == 5 || variant == 6) {
    // OrbitDE::start stops early inside the cardioid, so set up the iteration state by hand
    OrbitDE o;
    o.x = x; o.y = y; o.zx = x; o.zy = y; o.dx = 1; o.dy = 0; o.dexp = 0;
    o.min_r2 = x * x + y * y; o.candidate = 1; o.next_newton = variant == 6 ? 64 : 1 << 30;
    o.cx = x; o.cy = y; o.check_n = 1; o.next_check = 16; o.n = 1;
    o.run(kSteps, kSteps);
    return double(o.n) + o.zx + o.r.dist;
  } else if constexpr (variant == 4) {
    Orbit<double> o;
    o.start(x, y, int64_t(1) << 40);  // No Newton.  Reports the cardioid, but initializes the state first.
    o.status = 0; o.n = 1;
    o.run(kSteps, kSteps, 256);
    return double(o.n) + o.zx;
  } else {
    double zx = x, zy = y, min_r2 = 1e300, cx = x, cy = y;
    int64_t candidate = 0, next_check = 16, n = 1;
    for (; n <= kSteps; n++) {
      if constexpr (variant >= 1) {
        if (zx * zx + zy * zy > 18446744073709551616.0) break;
      }
      const double t = zx * zx - zy * zy + x;
      zy = 2 * zx * zy + y;
      zx = t;
      if constexpr (variant >= 2) {
        const double r2 = zx * zx + zy * zy;
        if (r2 < min_r2) { min_r2 = r2; candidate = n; }
      }
      if constexpr (variant >= 3) {
        const double dx = zx - cx, dy = zy - cy;
        if (dx * dx + dy * dy < 1e-300) break;  // Never true for these orbits
        if (n == next_check) { cx = zx; cy = zy; next_check *= 2; }
      }
    }
    return zx + zy + double(candidate) + double(n);
  }
}

template<int variant> struct Bench {
  double* out;
  __host__ __device__ void operator()(const int64_t i) const {
    const double theta = 2 * M_PI * 0.6180339887498949, r = 1 - 1e-9 * (1 + double(i % 1024) / 1024);
    const double mx = r * std::cos(theta), my = r * std::sin(theta);
    const double x = mx / 2 - (mx * mx - my * my) / 4, y = my / 2 - mx * my / 2;
    out[i] = orbit<variant>(x, y);
  }
};

// The same fixed-length orbits through run_orbits (persistent threads, claims, bursts), for engine overhead
struct SiegelTask {
  typedef Orbit<double> State;
  int64_t burst;
  int min_blocks;
  double* out;
  __host__ __device__ bool start(State& o, const int64_t i) const {
    const double theta = 2 * M_PI * 0.6180339887498949, r = 1 - 1e-9 * (1 + double(i % 1024) / 1024);
    const double mx = r * std::cos(theta), my = r * std::sin(theta);
    o.start(mx / 2 - (mx * mx - my * my) / 4, my / 2 - mx * my / 2, int64_t(1) << 40);
    o.status = 0; o.n = 1;
    return false;  // Iterate even though start() reports the cardioid
  }
  __host__ __device__ bool run(State& o) const { return o.run(kSteps, burst, 256); }
  __host__ __device__ int64_t iters(const State& o) const { return o.iters(); }
  __host__ __device__ int64_t progress(const State& o) const { return o.n; }
  __host__ __device__ void finish(const State& o, const int64_t i) const { out[i] = o.zx; }
};

void run_engine(const int64_t n, const bool cuda, const int64_t burst) {
  Mem<double> out(n, cuda);
  run_orbits(SiegelTask{burst, 3, out.p}, n, cuda);  // Warm up
  const auto t0 = std::chrono::steady_clock::now();
  const auto stats = run_orbits(SiegelTask{burst, 3, out.p}, n, cuda);
  const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  print("  run_orbits, burst %-4d        %.3g iterations/s  (%.3g s, %d overflowed)", burst,
        double(n) * kSteps / secs, secs, stats.overflow);
}

template<int variant> void run(const int64_t n, const bool cuda, const char* name) {
  Mem<double> out(n, cuda);
  for_each(n, Bench<variant>{out.p}, cuda);  // Warm up
  double h;
  out.to_host(&h, 1);
  const auto t0 = std::chrono::steady_clock::now();
  for_each(n, Bench<variant>{out.p}, cuda);
  out.to_host(&h, 1);
  const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  print("  %-28s %.3g iterations/s  (%.3g s)", name, double(n) * kSteps / secs, secs);
}

}  // namespace
}  // namespace mandelbrot

int main(const int argc, const char** argv) {
  using namespace mandelbrot;
  const bool cuda = argc > 1 && std::string(argv[1]) == "--cuda";
  const int64_t n = cuda ? 8 << 20 : int64_t(cpu_threads()) << 10;
  print("orbit_bench: %d orbits of %d steps on %s", n, kSteps, cuda ? "cuda" : "cpu");
  run<0>(n, cuda, "bare z^2 + c");
  run<1>(n, cuda, "+ escape test");
  run<2>(n, cuda, "+ atom-domain minimum");
  run<3>(n, cuda, "+ Brent check");
  run<4>(n, cuda, "Orbit::run (no Newton)");
  run<5>(n, cuda, "OrbitDE::run (no Newton)");
  run<6>(n, cuda, "OrbitDE::run (Newton from 64)");
  for (const int64_t burst : {64, 1024}) run_engine(n, cuda, burst);
  return 0;
}
