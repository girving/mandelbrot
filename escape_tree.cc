// Adaptive quadtree Monte Carlo for the areas of {c : g_M(c) < 2^-k}
//
// Only cells straddling the boundary of these sets contribute variance, and at fine scales they are a
// tiny fraction of all cells.  Starting from a base grid over [-2, 0.5] × [0, 1.2], each node classifies
// its center with distance estimates (escape_de).  If the cell fits well inside the certified disk, it is
// decided exactly: interior cells count fully; exterior cells count per threshold using Harnack bounds on
// g.  Otherwise the node splits into 4 children, down to `depth` levels, where leaves are estimated from R
// independent points (one per replica).  Certified leaves have no variance; sampled leaves give an
// unbiased estimate, and the replica spread estimates the statistical error.

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

constexpr int R = 8;

struct Params {
  int64_t base;      // Base grid size per axis
  int depth;         // Maximum refinement levels below the base grid
  double safety;     // Certify a cell if its half-diagonal is at most dist / safety
  int64_t max_iter;
  vector<int> ks;    // Thresholds 2^-k
};

struct Stats {
  vector<double> sum;         // sum[r * K + k]: replica r's area estimate for threshold k
  vector<int64_t> leaves;     // Sampled leaves per depth
  vector<int64_t> exact;      // Certified leaves per depth
  int64_t samples = 0, iters = 0, center_iters = 0, leaf_iters = 0, leaf_exterior_iters = 0;
  int64_t why[5] = {};  // At max depth: interior dist small, exterior dist small, exterior band, no certificate, other
  void add(const Stats& o) {
    for (size_t i = 0; i < sum.size(); i++) sum[i] += o.sum[i];
    for (size_t i = 0; i < leaves.size(); i++) { leaves[i] += o.leaves[i]; exact[i] += o.exact[i]; }
    samples += o.samples; iters += o.iters;
    center_iters += o.center_iters; leaf_iters += o.leaf_iters; leaf_exterior_iters += o.leaf_exterior_iters;
    for (int i = 0; i < 5; i++) why[i] += o.why[i];
  }
};

struct Worker {
  const Params& p;
  std::mt19937_64 rng;
  std::uniform_real_distribution<double> u{0, 1};
  Stats st;

  Worker(const Params& p, const uint64_t seed) : p(p), rng(seed) {
    st.sum.assign(R * p.ks.size(), 0);
    st.leaves.assign(p.depth + 1, 0);
    st.exact.assign(p.depth + 1, 0);
  }

  Escape sample(const double x, const double y, const double w, const double h) {
    const auto e = escape(x + u(rng) * w, y + u(rng) * h, p.max_iter);
    st.samples++;
    st.iters += e.iters;
    return e;
  }

  void node(const double x, const double y, const double w, const double h, const int d) {
    const int K = p.ks.size();
    // Certify the whole cell from its center
    const auto c = escape_de(x + w / 2, y + h / 2, p.max_iter);
    st.samples++;
    st.iters += c.e.iters;
    st.center_iters += c.e.iters;
    const double r = 0.5 * std::hypot(w, h), a = w * h;
    if (c.dist > 0 && r * p.safety <= c.dist) {
      if (c.e.steps < 0) {  // Whole cell in one hyperbolic component: below every threshold
        st.exact[d]++;
        for (int rr = 0; rr < R; rr++) for (int k = 0; k < K; k++) st.sum[rr * K + k] += a;
        return;
      }
      // Exterior: Harnack on the disk of radius dist gives g ∈ [lo, hi] · g0 on the cell
      const double t = r / c.dist, lo = std::log2((1 - t) / (1 + t)), hi = -lo;
      bool certain = true;
      for (int k = 0; k < K; k++) certain &= c.e.log2g + hi < -p.ks[k] || c.e.log2g + lo >= -p.ks[k];
      if (certain) {
        st.exact[d]++;
        for (int k = 0; k < K; k++)
          if (c.e.log2g + hi < -p.ks[k])
            for (int rr = 0; rr < R; rr++) st.sum[rr * K + k] += a;
        return;
      }
    }
    if (d == p.depth) {
      const int reason = !(c.dist > 0) ? 3 : r * p.safety > c.dist ? (c.e.steps < 0 ? 0 : 1) : c.e.steps > 0 ? 2 : 4;
      st.why[reason]++;
    }
    if (d < p.depth) {
      const double w2 = w / 2, h2 = h / 2;
      node(x, y, w2, h2, d + 1);
      node(x + w2, y, w2, h2, d + 1);
      node(x, y + h2, w2, h2, d + 1);
      node(x + w2, y + h2, w2, h2, d + 1);
      return;
    }
    // Leaf: one independent point per replica
    st.leaves[d]++;
    for (int rr = 0; rr < R; rr++) {
      const auto e = sample(x, y, w, h);
      st.leaf_iters += e.iters;
      if (e.steps > 0) st.leaf_exterior_iters += e.iters;
      for (int k = 0; k < K; k++) st.sum[rr * K + k] += a * below(e, p.ks[k]);
    }
  }
};

void run(const Params& p, const uint64_t seed) {
  const double x0 = -2, x1 = 0.5, y0 = 0, y1 = 1.2;
  const double w = (x1 - x0) / p.base, h = (y1 - y0) / p.base;
  const int K = p.ks.size();
  slow_assert(K <= 32);
  const auto t0 = wall_time();
  std::atomic<int64_t> next(0);
  std::mutex mu;
  Stats total;
  total.sum.assign(R * K, 0);
  total.leaves.assign(p.depth + 1, 0);
  total.exact.assign(p.depth + 1, 0);
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
  print("base %d, depth %d (effective grid %d), safety %g, max_iter %d, seed %d: %.1f s", p.base, p.depth,
        p.base << p.depth, p.safety, p.max_iter, seed, secs);
  print("  %.3g samples, %.3g iterations: centers %.3g, leaf samples %.3g (exterior %.3g)", double(total.samples),
        double(total.iters), double(total.center_iters), double(total.leaf_iters), double(total.leaf_exterior_iters));
  string ls = "  sampled leaves by depth:", es = "  certified leaves by depth:";
  for (const auto n : total.leaves) ls += tfm::format(" %.3g", double(n));
  for (const auto n : total.exact) es += tfm::format(" %.3g", double(n));
  print(es);
  print("  uncertified at max depth: interior dist small %.3g, exterior dist small %.3g, exterior band %.3g, "
        "no certificate %.3g, other %.3g", double(total.why[0]), double(total.why[1]), double(total.why[2]),
        double(total.why[3]), double(total.why[4]));
  print(ls);
  print("   k      area{g < 2^-k}      std err");
  for (int k = 0; k < K; k++) {
    double s = 0, s2 = 0;
    for (int r = 0; r < R; r++) {
      const double a = 2 * total.sum[r * K + k];  // Symmetry about the real axis
      s += a; s2 += a * a;
    }
    const double mean = s / R, sd = std::sqrt(std::max(0.0, s2 / R - mean * mean) / (R - 1));
    print("%8d   %.10f   %.2e", p.ks[k], mean, sd);
  }
}

}  // namespace
}  // namespace mandelbrot

int main(const int argc, const char** argv) {
  using namespace mandelbrot;
  try {
    slow_assert(argc >= 7, "usage: %s <base> <depth> <safety> <max_iter> <seed> <k...>", argv[0]);
    Params p{atoll(argv[1]), atoi(argv[2]), atof(argv[3]), atoll(argv[4]), {}};
    const uint64_t seed = atoll(argv[5]);
    for (int i = 6; i < argc; i++) p.ks.push_back(atoi(argv[i]));
    for (const int k : p.ks) slow_assert(k + 8 <= p.max_iter, "need max_iter ≥ k + 8 for k = %d", k);
    run(p, seed);
    return 0;
  } catch (const std::exception& e) {
    die(e.what());
  }
}
