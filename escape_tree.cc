// Adaptive quadtree Monte Carlo for the areas of {c : g_M(c) < 2^-k}
//
// Only cells straddling the boundary of these sets contribute variance, and at fine scales they are a
// tiny fraction of all cells.  Starting from a base grid over [-2, 0.5] × [0, 1.2], each node classifies
// its center with distance estimates (escape_de).  If the cell fits well inside the certified disk, it is
// decided exactly: interior cells count fully; exterior cells count per threshold using Harnack bounds on
// g.  Otherwise the node splits into 4 children, down to `depth` levels, where leaves are estimated from
// m independent uniform points.  Certified leaves have no variance.  Each sampled leaf gives an unbiased
// estimate with an unbiased variance estimate s^2 / m, and leaves are independent, so the total variance
// is the sum over leaves.
//
// Besides each area A(k), we report the differences D = A(k_i) - A(k_{i+1}) between consecutive thresholds
// with their own variances, for multilevel estimates A(K) = A(k_0) - ∑ D computed in separate runs.

#include "debug.h"
#include "escape.h"
#include "print.h"
#include "wall_time.h"
#include <atomic>
#include <cmath>
#include <mutex>
#include <random>
#include <thread>
#include <vector>
namespace mandelbrot {
namespace {

using std::vector;

struct Params {
  int64_t base;      // Base grid size per axis
  int depth;         // Maximum refinement levels below the base grid
  double safety;     // Certify a cell if its half-diagonal is at most dist / safety
  int m;             // Fresh samples per sampled leaf
  int pilots;        // Pilot samples per uncertified leaf (0 disables roulette)
  double q;          // Roulette probability for leaves whose pilots agree with the center
  int64_t max_iter;
  vector<int> ks;    // Thresholds 2^-k, increasing
};

struct Stats {
  vector<double> est, var, dest, dvar;  // Areas and consecutive differences, with variances
  vector<int64_t> leaves, exact;        // Sampled and certified leaves per depth
  int64_t samples = 0, iters = 0, center_iters = 0, leaf_iters = 0;
  void init(const int K, const int D) {
    est.assign(K, 0); var.assign(K, 0); dest.assign(K, 0); dvar.assign(K, 0);
    leaves.assign(D + 1, 0); exact.assign(D + 1, 0);
  }
  void add(const Stats& o) {
    for (size_t i = 0; i < est.size(); i++) {
      est[i] += o.est[i]; var[i] += o.var[i]; dest[i] += o.dest[i]; dvar[i] += o.dvar[i];
    }
    for (size_t i = 0; i < leaves.size(); i++) { leaves[i] += o.leaves[i]; exact[i] += o.exact[i]; }
    samples += o.samples; iters += o.iters; center_iters += o.center_iters; leaf_iters += o.leaf_iters;
  }
};

struct Worker {
  const Params& p;
  std::mt19937_64 rng;
  std::uniform_real_distribution<double> u{0, 1};
  Stats st;
  vector<double> xs;  // Scratch: per-sample indicators

  Worker(const Params& p, const uint64_t seed) : p(p), rng(seed) {
    st.init(p.ks.size(), p.depth);
    xs.resize(p.m * p.ks.size());
  }

  // Exact contribution of a certified cell: below_k[k] says whether the whole cell is below threshold k
  void certified(const double a, const vector<bool>& below_k) {
    const int K = p.ks.size();
    for (int k = 0; k < K; k++) {
      st.est[k] += a * below_k[k];
      if (k + 1 < K) st.dest[k] += a * (double(below_k[k]) - double(below_k[k + 1]));
    }
  }

  void node(const double x, const double y, const double w, const double h, const int d) {
    const int K = p.ks.size();
    const double a = w * h;
    // Certify the whole cell from its center
    const auto c = escape_de(x + w / 2, y + h / 2, p.max_iter);
    st.samples++;
    st.iters += c.e.iters;
    st.center_iters += c.e.iters;
    const double r = 0.5 * std::hypot(w, h);
    if (c.dist > 0 && r * p.safety <= c.dist) {
      vector<bool> below_k(K, true);
      bool certain = true;
      if (c.e.steps > 0) {
        // Exterior: Harnack on the disk of radius dist gives g ∈ [lo, hi] · g0 on the cell
        const double t = r / c.dist, lo = std::log2((1 - t) / (1 + t)), hi = -lo;
        for (int k = 0; k < K; k++) {
          below_k[k] = c.e.log2g + hi < -p.ks[k];
          certain &= below_k[k] || c.e.log2g + lo >= -p.ks[k];
        }
      }
      if (certain) {
        st.exact[d]++;
        certified(a, below_k);
        return;
      }
    }
    if (d < p.depth) {
      const double w2 = w / 2, h2 = h / 2;
      node(x, y, w2, h2, d + 1);
      node(x + w2, y, w2, h2, d + 1);
      node(x, y + h2, w2, h2, d + 1);
      node(x + w2, y + h2, w2, h2, d + 1);
      return;
    }
    // Uncertified leaf.  Pilot points decide between full sampling (if they disagree with the center at
    // any threshold) and roulette with the center as control.  Estimates use only fresh points, which are
    // independent of the decision, so they are unbiased.
    st.leaves[d]++;
    bool agree = p.pilots > 0;
    for (int i = 0; i < p.pilots && agree; i++) {
      const auto e = escape(x + u(rng) * w, y + u(rng) * h, p.max_iter);
      st.samples++;
      st.iters += e.iters;
      st.leaf_iters += e.iters;
      for (int k = 0; k < K; k++) agree &= below(e, p.ks[k]) == below(c.e, p.ks[k]);
    }
    const double qq = agree ? p.q : 1;
    if (agree && !(u(rng) < qq)) {
      // Skipped by roulette: the estimate is the control
      for (int k = 0; k < K; k++) {
        st.est[k] += a * below(c.e, p.ks[k]);
        if (k + 1 < K) st.dest[k] += a * (double(below(c.e, p.ks[k])) - double(below(c.e, p.ks[k + 1])));
      }
      return;
    }
    for (int i = 0; i < p.m; i++) {
      const auto e = escape(x + u(rng) * w, y + u(rng) * h, p.max_iter);
      st.samples++;
      st.iters += e.iters;
      st.leaf_iters += e.iters;
      for (int k = 0; k < K; k++) xs[i * K + k] = below(e, p.ks[k]);
    }
    for (int k = 0; k < K; k++)
      for (int diff = 0; diff < (k + 1 < K ? 2 : 1); diff++) {
        const double ctrl = diff ? double(below(c.e, p.ks[k])) - double(below(c.e, p.ks[k + 1]))
                                 : double(below(c.e, p.ks[k]));
        double s = 0, s2 = 0;
        for (int i = 0; i < p.m; i++) {
          const double v = diff ? xs[i * K + k] - xs[i * K + k + 1] : xs[i * K + k];
          s += v; s2 += v * v;
        }
        const double mean = s / p.m;
        double est, var;
        if (qq < 1) {
          // Roulette: ctrl + (mean - ctrl) / q, with the conservative variance (mean - ctrl)^2 / q^2
          est = ctrl + (mean - ctrl) / qq;
          var = (mean - ctrl) * (mean - ctrl) / (qq * qq);
        } else {
          est = mean;
          var = (s2 - s * mean) / (p.m - 1) / p.m;
        }
        (diff ? st.dest : st.est)[k] += a * est;
        (diff ? st.dvar : st.var)[k] += a * a * var;
      }
  }
};

void run(const Params& p, const uint64_t seed) {
  const double x0 = -2, x1 = 0.5, y0 = 0, y1 = 1.2;
  const double w = (x1 - x0) / p.base, h = (y1 - y0) / p.base;
  const int K = p.ks.size();
  const auto t0 = wall_time();
  std::atomic<int64_t> next(0);
  std::mutex mu;
  Stats total;
  total.init(K, p.depth);
  vector<std::thread> pool;
  for (int t = 0; t < int(std::thread::hardware_concurrency()); t++)
    pool.emplace_back([&]() {
      for (int64_t i; (i = next++) < p.base;) {
        Worker wk(p, seed * 1000003 + i);
        for (int64_t j = 0; j < p.base; j++) wk.node(x0 + j * w, y0 + i * h, w, h, 0);
        std::lock_guard<std::mutex> g(mu);
        total.add(wk.st);
      }
    });
  for (auto& t : pool) t.join();
  const double secs = (wall_time() - t0).seconds();
  print("base %d, depth %d (effective grid %d), safety %g, %d samples/leaf, %d pilots, q %g, max_iter %d, "
        "seed %d: %.1f s", p.base, p.depth, p.base << p.depth, p.safety, p.m, p.pilots, p.q, p.max_iter, seed, secs);
  print("  %.3g samples, %.3g iterations (centers %.3g, leaves %.3g)", double(total.samples),
        double(total.iters), double(total.center_iters), double(total.leaf_iters));
  string ls = "  sampled leaves by depth:", es = "  certified leaves by depth:";
  for (const auto n : total.leaves) ls += tfm::format(" %.3g", double(n));
  for (const auto n : total.exact) es += tfm::format(" %.3g", double(n));
  print(es);
  print(ls);
  // Areas use the symmetry about the real axis: total = 2 × upper half
  print("   k      area{g < 2^-k}      std err");
  for (int k = 0; k < K; k++) print("%8d   %.10f   %.2e", p.ks[k], 2 * total.est[k], 2 * std::sqrt(total.var[k]));
  print("   k → k'           A(k) - A(k')     std err");
  for (int k = 0; k + 1 < K; k++)
    print("%8d → %-8d  %.6e   %.2e", p.ks[k], p.ks[k + 1], 2 * total.dest[k], 2 * std::sqrt(total.dvar[k]));
}

}  // namespace
}  // namespace mandelbrot

int main(const int argc, const char** argv) {
  using namespace mandelbrot;
  try {
    slow_assert(argc >= 10, "usage: %s <base> <depth> <safety> <samples/leaf> <pilots> <q> <max_iter> <seed> <k...>",
                argv[0]);
    Params p{atoll(argv[1]), atoi(argv[2]), atof(argv[3]), atoi(argv[4]), atoi(argv[5]), atof(argv[6]),
             atoll(argv[7]), {}};
    const uint64_t seed = atoll(argv[8]);
    for (int i = 9; i < argc; i++) p.ks.push_back(atoi(argv[i]));
    for (const int k : p.ks) slow_assert(k + 8 <= p.max_iter, "need max_iter ≥ k + 8 for k = %d", k);
    for (size_t i = 0; i + 1 < p.ks.size(); i++) slow_assert(p.ks[i] < p.ks[i + 1], "thresholds must increase");
    slow_assert(p.m >= 2, "need at least 2 samples per leaf for variance estimates");
    run(p, seed);
    return 0;
  } catch (const std::exception& e) {
    die(e.what());
  }
}
