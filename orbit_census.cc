// Census of long leaf orbits: which samples run long, how they end, and what they cost
//
// Samples points near the boundary, like escape_tree's leaves (16 points in a 4e-5 box around random points
// that escape after at least 1024 steps), and runs each to max_iter with the leaf pipeline's Newton
// schedule (first attempt at first_newton, deferred and settled as in the GPU rounds).  Reports total work,
// the orbits longer than max_iter / 16 by outcome (escaped, certified interior, hit max_iter) with the step
// they ended at, and how many samples are below the deepest threshold.  Classification counts must not depend
// on the Newton settings; work may.

#include "argparse.hpp"
#include "debug.h"
#include "engine.h"
#include "escape.h"
#include "print.h"
#include <atomic>
#include <cmath>
#include <mutex>
#include <random>
#include <thread>
#include <vector>
using namespace mandelbrot;

int main(const int argc, const char** argv) {
  try {
    argparse::ArgumentParser program("orbit_census");
    program.add_argument("--samples").scan<'i', int64_t>().default_value(int64_t(1000000));
    program.add_argument("--max-iter").scan<'i', int64_t>().default_value(int64_t(1) << 28);
    program.add_argument("--first-newton").scan<'i', int64_t>().default_value(int64_t(16384));
    program.add_argument("--max-period").scan<'i', int>().default_value(256);
    program.add_argument("--newton-iters").scan<'i', int>().default_value(30);
    program.add_argument("--newton-tol").scan<'g', double>().default_value(-1.0);
    program.add_argument("--newton-margin").scan<'g', double>().default_value(-1.0);
    program.add_argument("--seed").scan<'i', int64_t>().default_value(int64_t(100));
    program.parse_args(argc, argv);
    const int64_t samples = program.get<int64_t>("--samples"), max_iter = program.get<int64_t>("--max-iter"),
                  first_newton = program.get<int64_t>("--first-newton");
    const int max_period = program.get<int>("--max-period");
    const double tol = program.get<double>("--newton-tol");
    const NewtonOptions nw{program.get<int>("--newton-iters"), INFINITY, tol < 0 ? -1 : tol * tol, program.get<double>("--newton-margin")};
    slow_assert(max_iter < (int64_t(1) << 30), "max_iter must be below 2^30");
    const auto t0 = std::chrono::steady_clock::now();

    // Sample points, the same for any Newton settings
    std::vector<std::pair<double, double>> pts;
    std::mt19937_64 rng(program.get<int64_t>("--seed"));
    std::uniform_real_distribution<double> ux(-2, 0.5), uy(0, 1.2), u(-0.5, 0.5);
    while (int64_t(pts.size()) < samples) {
      const double x0 = ux(rng), y0 = uy(rng);
      // Iterations without Newton, so that the sample set does not depend on Newton's behavior
      Orbit<double> f;
      if (!f.start(x0, y0, int64_t(1) << 40)) f.finish(1 << 12);
      if (f.status != 1 || f.iters() < 1024) continue;  // Escaping slowly: near the boundary
      for (int s = 0; s < 16 && int64_t(pts.size()) < samples; s++) pts.push_back({x0 + 4e-5 * u(rng), y0 + 4e-5 * u(rng)});
    }

    const int64_t long_steps = max_iter / 16;
    const int octaves = int(std::log2(double(max_iter))) + 2;
    std::atomic<int64_t> next(0);
    std::mutex mu;
    int64_t below = 0, n_long = 0;
    std::vector<int64_t> esc(octaves), cert(octaves), capped(1);
    double work = 0, long_work = 0, newton_secs = 0, total_secs = 0;
    int64_t newtons = 0, certified = 0;
    std::vector<std::thread> pool;
    for (int t = 0; t < cpu_threads(); t++)
      pool.emplace_back([&]() {
        int64_t b = 0, l = 0, cp = 0;
        std::vector<int64_t> e(octaves), c(octaves);
        double w = 0, lw = 0, ns = 0;
        int64_t nn = 0, nc = 0;
        const auto s0 = std::chrono::steady_clock::now();
        for (int64_t i; (i = next.fetch_add(1)) < samples;) {
          Orbit<double> o;
          if (o.start(pts[i].first, pts[i].second, first_newton)) { b++; continue; }
          bool done = false;
          while (!done) {
            done = o.run(max_iter, 1 << 20);
            if (o.pending()) {
              const bool due = o.status == 4;
              const auto n0 = std::chrono::steady_clock::now();
              done = o.settle(max_iter, max_period, nw);
              if (due) {
                ns += std::chrono::duration<double>(std::chrono::steady_clock::now() - n0).count();
                nn++;
                nc += o.status == 2;
              }
            }
          }
          w += double(o.n);
          b += o.status != 1 || escaped_below(o.n, double(o.cx), int(max_iter - 8));  // Below 2^-(max_iter - 8)
          if (o.n > long_steps) {
            l++; lw += double(o.n);
            const int oct = int(std::log2(double(o.n)));
            if (o.status == 1) e[oct]++;
            else if (o.status == 2) c[oct]++;
            else cp++;
          }
        }
        std::lock_guard<std::mutex> g(mu);
        total_secs += std::chrono::duration<double>(std::chrono::steady_clock::now() - s0).count();
        newton_secs += ns; newtons += nn; certified += nc;
        below += b; n_long += l; capped[0] += cp; work += w; long_work += lw;
        for (int k = 0; k < octaves; k++) { esc[k] += e[k]; cert[k] += c[k]; }
      });
    for (auto& t : pool) t.join();
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    print("orbit_census: %d samples, max_iter %d, first Newton %d, max period %d, Newton tol %g, margin %g, "
          "%d threads: %.1f s", samples, max_iter, first_newton, max_period, tol, nw.margin, cpu_threads(), secs);
    print("  work %.4g steps, %.1f%% in %d orbits longer than %d steps (%d hit max_iter); %d samples below "
          "2^-(max_iter - 8)", work, 100 * long_work / work, n_long, long_steps, capped[0], below);
    print("  Newton: %d attempts, %d certified, %.1f%% of thread time", newtons, certified,
          100 * newton_secs / total_secs);
    for (int k = int(std::log2(double(long_steps))); k < octaves; k++)
      if (esc[k] || cert[k]) print("    ended at 2^%d: %d escaped, %d certified", k, esc[k], cert[k]);
    return 0;
  } catch (const std::exception& e) {
    die(e.what());
  }
}
