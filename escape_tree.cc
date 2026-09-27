// Adaptive quadtree Monte Carlo for the areas of {c : g_M(c) < 2^-k}
//
// Only cells straddling the boundary of these sets contribute variance, and at fine scales they are a
// tiny fraction of all cells.  Starting from a base grid over [-2, 0.5] × [0, 1.2], each node classifies
// its center with distance estimates (escape_de).  If the cell fits well inside the certified disk, it is
// decided exactly: interior cells count fully; exterior cells count per threshold using Harnack bounds on
// g.  Otherwise the node splits into 4 children, down to `depth` levels.  Uncertified leaves at max depth are
// collected and sampled in batches (escape_batch.h, on CPU threads or the GPU): m points per leaf, in groups
// of strata^2 points jittered on a strata × strata grid.  Group means are iid, so each leaf gives an
// unbiased estimate with an unbiased variance estimate, and leaves are independent, so the total variance
// is the sum over leaves.
//
// Besides each area A(k), we report the differences D = A(k_i) - A(k_{i+1}) between consecutive thresholds
// with their own variances.  With --prec compare, every leaf sample is classified with both float and double
// orbits, and we also report the paired difference A_float(k) - A_double(k), whose variance comes only from
// samples that flip, to measure the bias of low precision.

#include "argparse.hpp"
#include "debug.h"
#include "escape.h"
#include "escape_batch.h"
#include "print.h"
#include "wall_time.h"
#include <atomic>
#include <cmath>
#include <mutex>
#include <thread>
#include <vector>
namespace mandelbrot {
namespace {

using std::vector;

struct Params {
  int64_t base;      // Base grid size per axis
  int depth;         // Maximum refinement levels below the base grid
  double safety;     // Certify a cell if its half-diagonal is at most dist / safety
  int m;             // Samples per uncertified leaf
  int strata;        // Each group of strata^2 leaf samples is jittered on a strata × strata grid
  int64_t max_iter;
  uint64_t seed;
  string prec;       // double, float, or compare
  bool cuda;
  int64_t batch;     // Leaves per sampling batch
  vector<int> ks;    // Thresholds 2^-k, increasing
};

// Sums of estimates and variances for areas and consecutive differences
struct Sums {
  vector<double> est, var, dest, dvar;
  void init(const int K) { est.assign(K, 0); var.assign(K, 0); dest.assign(K, 0); dvar.assign(K, 0); }
};

// Tree traversal: certified cells go straight into the sums, uncertified leaves into a list
struct Walker {
  const Params& p;
  vector<double> est;            // Certified area below each threshold
  vector<int64_t> exact;         // Certified cells per depth
  vector<Leaf> leaves;
  int64_t samples = 0, iters = 0;

  Walker(const Params& p) : p(p), est(p.ks.size()), exact(p.depth + 1) {}

  void node(const double x, const double y, const double w, const double h, const int d) {
    const int K = p.ks.size();
    const auto c = escape_de(x + w / 2, y + h / 2, p.max_iter);
    samples++;
    iters += c.e.iters;
    const double r = 0.5 * std::hypot(w, h);
    if (c.dist > 0 && r * p.safety <= c.dist) {
      bool certain = true;
      vector<bool> below_k(K, true);
      if (c.e.steps > 0) {
        // Exterior: Harnack on the disk of radius dist gives g ∈ [lo, hi] · g0 on the cell
        const double t = r / c.dist, lo = std::log2((1 - t) / (1 + t)), hi = -lo;
        for (int k = 0; k < K; k++) {
          below_k[k] = c.e.log2g + hi < -p.ks[k];
          certain &= below_k[k] || c.e.log2g + lo >= -p.ks[k];
        }
      }
      if (certain) {
        exact[d]++;
        for (int k = 0; k < K; k++) est[k] += w * h * below_k[k];
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
    leaves.push_back({x, y, w, h});
  }
};

// Add one leaf's samples to the sums.  f(i, k) is sample i's value for threshold k.
template<class F> void reduce_leaf(Sums& s, const Params& p, const double a, F&& f) {
  const int K = p.ks.size(), ss = p.strata * p.strata, groups = p.m / ss;
  for (int k = 0; k < K; k++)
    for (int diff = 0; diff < (k + 1 < K ? 2 : 1); diff++) {
      double s1 = 0, s2 = 0;
      for (int gi = 0; gi < groups; gi++) {
        double v = 0;
        for (int i = gi * ss; i < (gi + 1) * ss; i++) v += diff ? f(i, k) - f(i, k + 1) : f(i, k);
        v /= ss;
        s1 += v; s2 += v * v;
      }
      const double mean = s1 / groups, var = (s2 - s1 * mean) / (groups - 1) / groups;
      (diff ? s.dest : s.est)[k] += a * mean;
      (diff ? s.dvar : s.var)[k] += a * a * var;
    }
}

void print_sums(const string& title, const Sums& s, const Params& p, const bool diffs) {
  // Areas use the symmetry about the real axis: total = 2 × upper half
  const int K = p.ks.size();
  print("  %s:", title);
  print("       k      area{g < 2^-k}      std err");
  for (int k = 0; k < K; k++) print("    %8d   %.10f   %.2e", p.ks[k], 2 * s.est[k], 2 * std::sqrt(s.var[k]));
  if (!diffs) return;
  print("       k → k'           A(k) - A(k')     std err");
  for (int k = 0; k + 1 < K; k++)
    print("    %8d → %-8d  %.6e   %.2e", p.ks[k], p.ks[k + 1], 2 * s.dest[k], 2 * std::sqrt(s.dvar[k]));
}

void run(const Params& p) {
  const double x0 = -2, x1 = 0.5, y0 = 0, y1 = 1.2;
  const double w = (x1 - x0) / p.base, h = (y1 - y0) / p.base;
  const int K = p.ks.size();
  const bool compare = p.prec == "compare", single = p.prec == "float";
  SampleParams sp{p.m, p.strata, p.seed, p.max_iter, K, {}};
  for (int k = 0; k < K; k++) sp.ks[k] = p.ks[k];

  // Double (or float) results, and with compare, float results and float - double differences
  Sums main, fl, delta;
  main.init(K); fl.init(K); delta.init(K);
  vector<double> certified(K);
  vector<int64_t> exact(p.depth + 1);
  int64_t n_leaves = 0, center_samples = 0, center_iters = 0, leaf_iters = 0, leaf_iters_f = 0, flips = 0;
  double tree_secs = 0, sample_secs = 0;
  const auto t0 = wall_time();

  int64_t rows_per_batch = 16;
  for (int64_t row0 = 0; row0 < p.base;) {
    const int64_t row1 = std::min(p.base, row0 + rows_per_batch);
    // Walk rows [row0, row1) in parallel
    const auto t1 = wall_time();
    std::atomic<int64_t> next(row0);
    std::mutex mu;
    vector<Leaf> leaves;
    vector<std::thread> pool;
    for (int t = 0; t < int(std::thread::hardware_concurrency()); t++)
      pool.emplace_back([&]() {
        for (int64_t i; (i = next++) < row1;) {
          Walker wk(p);
          for (int64_t j = 0; j < p.base; j++) wk.node(x0 + j * w, y0 + i * h, w, h, 0);
          std::lock_guard<std::mutex> g(mu);
          for (int k = 0; k < K; k++) certified[k] += wk.est[k];
          for (int d = 0; d <= p.depth; d++) exact[d] += wk.exact[d];
          center_samples += wk.samples;
          center_iters += wk.iters;
          leaves.insert(leaves.end(), wk.leaves.begin(), wk.leaves.end());
        }
      });
    for (auto& t : pool) t.join();
    tree_secs += (wall_time() - t1).seconds();

    // Sample the leaves
    const auto t2 = wall_time();
    const auto sample = [&](const bool f, vector<uint32_t>& bits) {
      bits.resize(leaves.size() * p.m);
      return p.cuda ? (f ? sample_leaves_cuda<float>(leaves, sp, bits) : sample_leaves_cuda<double>(leaves, sp, bits))
                    : (f ? sample_leaves_cpu<float>(leaves, sp, bits) : sample_leaves_cpu<double>(leaves, sp, bits));
    };
    vector<uint32_t> bits, fbits;
    leaf_iters += sample(single, bits);
    if (compare) leaf_iters_f += sample(true, fbits);
    sample_secs += (wall_time() - t2).seconds();
    for (size_t l = 0; l < leaves.size(); l++) {
      const double a = leaves[l].w * leaves[l].h;
      const uint32_t* b = bits.data() + l * p.m;
      reduce_leaf(main, p, a, [b](int i, int k) { return double((b[i] >> k) & 1); });
      if (compare) {
        const uint32_t* bf = fbits.data() + l * p.m;
        reduce_leaf(fl, p, a, [bf](int i, int k) { return double((bf[i] >> k) & 1); });
        reduce_leaf(delta, p, a, [b, bf](int i, int k) { return double((bf[i] >> k) & 1) - double((b[i] >> k) & 1); });
        for (int i = 0; i < p.m; i++) flips += b[i] != bf[i];
      }
    }
    n_leaves += leaves.size();

    // Aim for about p.batch leaves per batch
    const double per_row = double(leaves.size()) / double(row1 - row0);
    rows_per_batch = std::max<int64_t>(1, int64_t(p.batch / std::max(1.0, per_row)));
    row0 = row1;
  }
  for (int k = 0; k < K; k++) { main.est[k] += certified[k]; fl.est[k] += certified[k]; }

  const double secs = (wall_time() - t0).seconds();
  print("base %d, depth %d (effective grid %d), safety %g, %d samples/leaf, strata %d, max_iter %d, seed %d, "
        "prec %s, %s: %.1f s (tree %.1f s, sampling %.1f s)", p.base, p.depth, p.base << p.depth, p.safety, p.m,
        p.strata, p.max_iter, p.seed, p.prec, p.cuda ? "cuda" : "cpu", secs, tree_secs, sample_secs);
  print("  centers: %.3g samples, %.3g iterations; leaves: %.3g leaves, %.3g samples, %.3g iterations",
        double(center_samples), double(center_iters), double(n_leaves), double(n_leaves * p.m),
        double(leaf_iters + leaf_iters_f));
  string es = "  certified cells by depth:";
  for (const auto n : exact) es += tfm::format(" %.3g", double(n));
  print(es);
  print_sums(compare ? "double" : p.prec, main, p, true);
  if (compare) {
    print_sums("float", fl, p, false);
    print("  float - double (%d flipped samples of %d):", flips, n_leaves * p.m);
    for (int k = 0; k < K; k++)
      print("    %8d   %+.3e   %.2e", p.ks[k], 2 * delta.est[k], 2 * std::sqrt(delta.var[k]));
  }
}

}  // namespace
}  // namespace mandelbrot

int main(const int argc, const char** argv) {
  using namespace mandelbrot;
  try {
    argparse::ArgumentParser program("escape_tree");
    program.add_argument("ks").help("thresholds g < 2^-k, increasing").remaining().scan<'i', int>();
    program.add_argument("--base").help("base grid size per axis").scan<'i', int64_t>().default_value(int64_t(1000));
    program.add_argument("--depth").help("refinement levels").scan<'i', int>().default_value(5);
    program.add_argument("--safety").help("certify if half-diagonal ≤ dist / safety").scan<'g', double>()
        .default_value(4.0);
    program.add_argument("--m").help("samples per leaf").scan<'i', int>().default_value(16);
    program.add_argument("--strata").help("jitter groups of strata^2 samples").scan<'i', int>().default_value(2);
    program.add_argument("--max-iter").scan<'i', int64_t>().default_value(int64_t(1) << 20);
    program.add_argument("--seed").scan<'i', int64_t>().default_value(int64_t(1));
    program.add_argument("--prec").help("leaf orbit precision: double, float, or compare")
        .default_value(string("double"));
    program.add_argument("--cuda").help("sample leaves on the GPU").default_value(false).implicit_value(true);
    program.add_argument("--batch").help("leaves per sampling batch").scan<'i', int64_t>()
        .default_value(int64_t(1) << 22);
    program.parse_args(argc, argv);

    Params p{program.get<int64_t>("--base"), program.get<int>("--depth"), program.get<double>("--safety"),
             program.get<int>("--m"), program.get<int>("--strata"), program.get<int64_t>("--max-iter"),
             uint64_t(program.get<int64_t>("--seed")), program.get<string>("--prec"), program.get<bool>("--cuda"),
             program.get<int64_t>("--batch"), program.get<vector<int>>("ks")};
    slow_assert(p.prec == "double" || p.prec == "float" || p.prec == "compare", "bad --prec %s", p.prec);
    slow_assert(p.ks.size() <= 32, "at most 32 thresholds");
    for (const int k : p.ks) slow_assert(k + 8 <= p.max_iter, "need max_iter ≥ k + 8 for k = %d", k);
    for (size_t i = 0; i + 1 < p.ks.size(); i++) slow_assert(p.ks[i] < p.ks[i + 1], "thresholds must increase");
    slow_assert(p.strata >= 1 && p.m % (p.strata * p.strata) == 0 && p.m / (p.strata * p.strata) >= 2,
                "need m a multiple of strata^2 with at least 2 groups for variance estimates");
    run(p);
    return 0;
  } catch (const std::exception& e) {
    die(e.what());
  }
}
