// Monte Carlo area of the fattened Mandelbrot sets {c : g_M(c) < 2^-k}
//
// Samples a jittered grid over [-2, 0.5] × [0, 1.2] (M is symmetric about the real axis), classifies each
// point with escape(), and counts points below each threshold 2^-k.  Each grid cell gets 8 independent
// jittered points, one per replica, so each replica is a full jittered-grid estimate and their spread
// estimates the statistical error.

#include "debug.h"
#include "escape.h"
#include "print.h"
#include "wall_time.h"
#include <atomic>
#include <cmath>
#include <random>
#include <thread>
#include <vector>
namespace mandelbrot {
namespace {

using std::vector;

void run(const int64_t G, const int64_t max_iter, const vector<int>& ks, const uint64_t seed) {
  const double x0 = -2, x1 = 0.5, y0 = 0, y1 = 1.2;
  const double hx = (x1 - x0) / G, hy = (y1 - y0) / G, area = 2 * (x1 - x0) * (y1 - y0);
  const int R = 8;  // Replicas
  const int K = ks.size();
  const auto t0 = wall_time();
  // counts[r][k] for replica r and threshold index k; last column counts never-escaping points
  vector<std::atomic<int64_t>> counts(R * (K + 1));
  vector<std::atomic<int64_t>> iters(1);
  std::atomic<int64_t> next(0);
  vector<std::thread> pool;
  const int threads = std::max(1u, std::thread::hardware_concurrency());
  for (int t = 0; t < threads; t++)
    pool.emplace_back([&]() {
      vector<int64_t> local(R * (K + 1));
      for (int64_t i; (i = next++) < G;) {
        std::mt19937_64 rng(seed * 1000003 + i);
        std::uniform_real_distribution<double> u(0, 1);
        std::fill(local.begin(), local.end(), 0);
        int64_t it = 0;
        for (int64_t j = 0; j < G; j++)
          for (int r = 0; r < R; r++) {
            const double x = x0 + (j + u(rng)) * hx, y = y0 + (i + u(rng)) * hy;
            const auto e = escape(x, y, max_iter);
            it += e.steps < 0 ? 0 : e.steps;
            for (int k = 0; k < K; k++) local[r * (K + 1) + k] += below(e, ks[k]);
            local[r * (K + 1) + K] += e.steps < 0;
          }
        for (int k = 0; k < R * (K + 1); k++) counts[k] += local[k];
        iters[0] += it;
      }
    });
  for (auto& t : pool) t.join();
  const double secs = (wall_time() - t0).seconds();
  print("G = %d (%.3g samples), max_iter = %d, seed %d: %.1f s, %.3g escape iterations total", G,
        double(G) * G * R, max_iter, seed, secs, double(iters[0]));
  print("   k      area{g < 2^-k}      std err    (replica spread)");
  const double per = area / (double(G) * G);
  for (int k = 0; k <= K; k++) {
    double s = 0, s2 = 0;
    for (int r = 0; r < R; r++) {
      const double a = counts[r * (K + 1) + k] * per;
      s += a; s2 += a * a;
    }
    const double mean = s / R, sd = std::sqrt(std::max(0.0, s2 / R - mean * mean) / (R - 1));
    if (k < K) print("%6d   %.10f   %.2e", ks[k], mean, sd);
    else print("   inf   %.10f   %.2e   (never escaped within max_iter)", mean, sd);
  }
}

}  // namespace
}  // namespace mandelbrot

int main(const int argc, const char** argv) {
  using namespace mandelbrot;
  try {
    slow_assert(argc >= 4, "usage: %s <grid> <max_iter> <seed> [k...]", argv[0]);
    const int64_t G = atoll(argv[1]), max_iter = atoll(argv[2]);
    const uint64_t seed = atoll(argv[3]);
    vector<int> ks;
    for (int i = 4; i < argc; i++) ks.push_back(atoi(argv[i]));
    if (ks.empty()) for (int k = 8; k <= max_iter - 8; k++) ks.push_back(k);
    for (const int k : ks) slow_assert(k + 8 <= max_iter, "need max_iter ≥ k + 8 for k = %d", k);
    run(G, max_iter, ks, seed);
    return 0;
  } catch (const std::exception& e) {
    die(e.what());
  }
}
