// Where do long orbits spend their steps?  An upper bound on what local approximations of many-iteration
// dynamics could skip.
//
// Samples leaf-like points near the boundary (as in orbit_census), keeps orbits longer than --long steps, and
// replays each in windows of --window steps.  At the start of each window it finds the smallest period
// q ≤ --max-period at which the current point returns within --ret (or the best return if --ret is 0), Newton-solves for the q-cycle nearby, and classifies the window:
//   repelling:  the orbit is within --near of a repelling cycle (candidates for Koenigs linearization jumps)
//   parabolic:  ... of a nearly neutral cycle with multiplier close to a root of unity (Fatou-coordinate gates)
//   attracting: ... of an attracting cycle (slow interior convergence)
//   critical:   the orbit passes within --near of the critical point 0 during the window (baby-copy regimes)
//   other:      none of these (chaotic wandering)
// It also reports runs of consecutive windows near the same cycle: steps inside long runs are what a jump
// through that cycle's dynamics could replace.

#include "argparse.hpp"
#include "debug.h"
#include "engine.h"
#include "escape.h"
#include "print.h"
#include <atomic>
#include <cmath>
#include <chrono>
#include <complex>
#include <mutex>
#include <random>
#include <thread>
#include <vector>
using namespace mandelbrot;
using std::string;
using std::vector;
typedef std::complex<double> C;

namespace {

enum Regime { kRepelling, kParabolic, kAttracting, kCritical, kOther, kRegimes };
const char* names[kRegimes] = {"repelling", "parabolic", "attracting", "critical", "other"};

// Newton for a q-cycle point near z; returns false if it does not converge
bool cycle(const C c, const C z, const int q, C& w, C& lambda) {
  w = z;
  for (int it = 0; it < 40; it++) {
    C f = w, d = 1;
    for (int k = 0; k < q; k++) { d = 2.0 * f * d; f = f * f + c; }
    const C step = (f - w) / (d - 1.0);
    w -= step;
    if (!(std::abs(step) < 1e2)) return false;
    if (std::abs(step) < 1e-13 * (1 + std::abs(w))) {
      C g = w; lambda = 1;
      for (int k = 0; k < q; k++) { lambda = 2.0 * g * lambda; g = g * g + c; }
      return true;
    }
  }
  return false;
}

// Distance from x to the nearest rational r/s with s ≤ 12
double rational_gap(const double x) {
  double best = 1;
  for (int s = 1; s <= 12; s++) best = std::min(best, std::abs(x * s - std::round(x * s)) / s);
  return best;
}

}  // namespace

int main(const int argc, const char** argv) {
  try {
    argparse::ArgumentParser program("orbit_regimes");
    program.add_argument("--samples").scan<'i', int64_t>().default_value(int64_t(2000000));
    program.add_argument("--max-iter").scan<'i', int64_t>().default_value(int64_t(1) << 24);
    program.add_argument("--long").help("analyze orbits longer than this").scan<'i', int64_t>()
        .default_value(int64_t(1) << 16);
    program.add_argument("--window").scan<'i', int>().default_value(4096);
    program.add_argument("--max-period").scan<'i', int>().default_value(1024);
    program.add_argument("--near").help("distance to a cycle point (or 0) that counts as near").scan<'g', double>()
        .default_value(0.05);
    program.add_argument("--ret").help("use the smallest period returning within this distance (0: best return)")
        .scan<'g', double>().default_value(1e-3);
    program.add_argument("--seed").scan<'i', int64_t>().default_value(int64_t(7));
    program.parse_args(argc, argv);
    const int64_t samples = program.get<int64_t>("--samples"), max_iter = program.get<int64_t>("--max-iter"),
                  long_steps = program.get<int64_t>("--long");
    const int window = program.get<int>("--window"), max_period = program.get<int>("--max-period");
    const double near = program.get<double>("--near"), ret = program.get<double>("--ret");
    const auto t0 = std::chrono::steady_clock::now();

    // Leaf-like sample points
    vector<std::pair<double, double>> pts;
    std::mt19937_64 rng(program.get<int64_t>("--seed"));
    std::uniform_real_distribution<double> ux(-2, 0.5), uy(0, 1.2), u(-0.5, 0.5);
    while (int64_t(pts.size()) < samples) {
      const double x0 = ux(rng), y0 = uy(rng);
      if (escape(x0, y0, 1 << 12).iters < 1024) continue;
      for (int s = 0; s < 16 && int64_t(pts.size()) < samples; s++) pts.push_back({x0 + 4e-5 * u(rng), y0 + 4e-5 * u(rng)});
    }

    std::atomic<int64_t> next(0);
    std::mutex mu;
    double steps[kRegimes] = {}, total = 0, all_work = 0, run_steps[5] = {};  // run_steps: in runs of ≥ 1, 4, 16, 64, 256 windows
    int64_t n_long = 0;
    vector<double> lam_hist(40);  // log10(| |λ| - 1 |) for repelling/parabolic windows, weighted by steps
    const int qb = 12, eb = 31;  // Period buckets 1, 2, 3-4, ..., > 1024; length octaves
    vector<double> by_q(eb * qb);  // Escaped long orbits longer than 2^e by bucket of the last hugged period
    vector<int64_t> nsteps(samples);  // Per-sample orbit length, for box clustering
    vector<std::thread> pool;
    for (int t = 0; t < cpu_threads(); t++)
      pool.emplace_back([&]() {
        double st[kRegimes] = {}, tot = 0, all = 0, rs[5] = {};
        int64_t nl = 0;
        vector<double> lh(40), bq(eb * qb);
        for (int64_t i; (i = next.fetch_add(1)) < samples;) {
          const double x = pts[i].first, y = pts[i].second;
          Orbit<double> o;
          if (o.start(x, y, 16384)) continue;
          o.finish(max_iter, 1024, NewtonOptions{30, INFINITY, 1e-20, 1e-6});
          all += double(o.n);
          nsteps[i] = o.n;
          if (o.n <= long_steps) continue;
          nl++;
          // Replay in windows
          const C c(x, y);
          C z = c;
          int64_t n = 1;
          int run = 0;
          C run_w = 0;
          int run_q = 0, last_q = 0;
          double run_len = 0;
          auto close_run = [&]() {
            const int thresholds[5] = {1, 4, 16, 64, 256};
            for (int k = 0; k < 5; k++) if (run >= thresholds[k]) rs[k] += run_len;
            run = 0; run_len = 0;
          };
          while (n < o.n) {
            const int64_t len = std::min<int64_t>(window, o.n - n);
            // Best return of the window's first point
            C w = z; double bd = INFINITY; int q = 0;
            for (int k = 1; k <= max_period; k++) {
              w = w * w + c;
              const double d = std::abs(w - z);
              if (ret > 0 ? d < ret : d < bd) { bd = d; q = k; if (ret > 0) break; }
            }
            C cw, lam;
            Regime r = kOther;
            const bool ok = q && cycle(c, z, q, cw, lam) && std::abs(cw - z) < near;
            if (ok) {
              last_q = q;
              const double a = std::abs(lam);
              const double rot = std::arg(lam) / (2 * M_PI);
              if (std::abs(a - 1) < 1e-3 && rational_gap(rot) < 1e-3) r = kParabolic;
              else if (a > 1) r = kRepelling;
              else r = kAttracting;
              if (a != 1) lh[std::min(39, std::max(0, int(20 + std::log10(std::abs(a - 1)))))] += double(len);
            }
            // Iterate the window, noting close approaches to 0
            double min_abs = INFINITY;
            for (int64_t k = 0; k < len; k++) { z = z * z + c; min_abs = std::min(min_abs, std::abs(z)); }
            n += len;
            if (r == kOther && min_abs < near) r = kCritical;
            st[r] += double(len);
            tot += double(len);
            // Runs near the same cycle
            if (ok && run && q == run_q && std::abs(cw - run_w) < 1e-6 * (1 + std::abs(cw))) {
              run++; run_len += double(len);
            } else {
              close_run();
              if (ok) { run = 1; run_len = double(len); run_q = q; run_w = cw; }
            }
          }
          close_run();
          if (o.status == 1 && last_q) {
            const int b = std::min(qb - 1, int(std::ceil(std::log2(double(last_q)))));
            for (int e = 0; e < eb && (int64_t(1) << e) < o.n; e++) bq[e * qb + b]++;
          }
        }
        std::lock_guard<std::mutex> g(mu);
        for (int k = 0; k < kRegimes; k++) steps[k] += st[k];
        for (int k = 0; k < 5; k++) run_steps[k] += rs[k];
        for (int k = 0; k < 40; k++) lam_hist[k] += lh[k];
        for (int k = 0; k < eb * qb; k++) by_q[k] += bq[k];
        total += tot; all_work += all; n_long += nl;
      });
    for (auto& t : pool) t.join();
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    print("orbit_regimes: %d samples, max_iter %d, %d orbits longer than %d (%.1f%% of all steps), window %d, "
          "near %g, ret %g, %d threads: %.1f s", samples, max_iter, n_long, long_steps, 100 * total / all_work, window, near, ret,
          cpu_threads(), secs);
    for (int k = 0; k < kRegimes; k++) print("  %-10s %5.1f%% of long-orbit steps", names[k], 100 * steps[k] / total);
    const int thresholds[5] = {1, 4, 16, 64, 256};
    for (int k = 0; k < 5; k++)
      print("  in runs of ≥ %3d windows near one cycle: %5.1f%%", thresholds[k], 100 * run_steps[k] / total);
    double lt = 0; for (double v : lam_hist) lt += v;
    print("  | |λ| - 1 | near cycles, weighted by steps:");
    for (int k = 0; k < 40; k++)
      if (lam_hist[k] > 0.005 * lt) print("    1e%-3d..1e%-3d %5.1f%%", k - 20, k - 19, 100 * lam_hist[k] / lt);
    // Two-phase sampling: do survivors cluster in the 16-sample boxes?
    print("  box clustering (16 samples per box, %d boxes):", samples / 16);
    for (int e = 12; (int64_t(1) << e) < max_iter; e += 2) {
      const int64_t T = int64_t(1) << e;
      int64_t boxes = 0, surv = 0;
      double work_in = 0, work_all = 0;
      for (int64_t b = 0; b + 16 <= samples; b += 16) {
        int s = 0; double w = 0;
        for (int k = 0; k < 16; k++) { s += nsteps[b + k] > T; w += double(nsteps[b + k]); }
        boxes += s > 0; surv += s; work_all += w; if (s) work_in += w;
      }
      print("    > 2^%-2d: %8.4f%% of boxes hold survivors, %5.2f per such box (random: %5.2f), %5.1f%% of work",
            e, 100.0 * boxes / (samples / 16), boxes ? double(surv) / boxes : 0.0,
            surv ? double(surv) / (samples / 16) * 16 / (1 - std::pow(1 - double(surv) / samples, 16)) / 16 : 0.0,
            100 * work_in / work_all);
    }
    // Per-component tails: which periods do deep escaping survivors hug?
    print("  escaped survivors by last hugged period (%% per row):\n    %-8s %8s  1     2     3-4   5-8   9-16  -32   -64   -128  -256  -512  -1024 >1024", "length", "count");
    for (int e = 12; e < eb; e += 2) {
      double t = 0; for (int b = 0; b < qb; b++) t += by_q[e * qb + b];
      if (t < 20) break;
      string line = tfm::format("    > 2^%-4d %8d ", e, int64_t(t));
      for (int b = 0; b < qb; b++) line += tfm::format(" %5.1f", 100 * by_q[e * qb + b] / t);
      print(line);
    }
    return 0;
  } catch (const std::exception& e) {
    die(e.what());
  }
}
