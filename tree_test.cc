// Tree pipeline tests

#include "tree.h"
#include "engine.h"
#include "escape.h"
#include "tests.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <tuple>
namespace mandelbrot {
namespace {

using std::map;
using std::tuple;

const double X0 = -2, X1 = 0.5, Y0 = 0, Y1 = 1.2;

TreeParams small_params() {
  TreeParams p;
  p.base = 40;
  p.depth = 3;
  p.m = 8;
  p.strata = 2;
  p.max_iter = 1 << 14;
  p.seed = 7;
  p.ks = {16, 256, 4096};
  return p;
}

// Straightforward recursive reference for the tree and the leaf estimates
struct Reference {
  const TreeParams& p;
  vector<int64_t> certified, exact;
  vector<tuple<int, int>> leaves;
  Reference(const TreeParams& p) : p(p), certified((p.depth + 1) * p.ks.size()), exact(p.depth + 1) {
    for (int iy = 0; iy < p.base; iy++)
      for (int ix = 0; ix < p.base; ix++)
        node(ix, iy, 0);
  }
  void node(const int ix, const int iy, const int d) {
    const int K = p.ks.size();
    const double w = (X1 - X0) / double(p.base << d), h = (Y1 - Y0) / double(p.base << d), r = 0.5 * std::hypot(w, h);
    const auto e = escape_de(X0 + (ix + 0.5) * w, Y0 + (iy + 0.5) * h, p.max_iter);
    if (e.dist > 0 && r * p.safety <= e.dist) {
      bool certain = true;
      vector<bool> below(K, true);
      if (e.e.steps > 0) {
        const double t = r / e.dist, lo = std::log2((1 - t) / (1 + t)), hi = -lo;
        for (int k = 0; k < K; k++) {
          below[k] = e.e.log2g + hi < -p.ks[k];
          certain &= below[k] || e.e.log2g + lo >= -p.ks[k];
        }
      }
      if (certain) {
        exact[d]++;
        for (int k = 0; k < K; k++) certified[d * K + k] += below[k];
        return;
      }
    }
    if (d < p.depth) {
      for (int j = 0; j < 4; j++) node(2 * ix + j % 2, 2 * iy + j / 2, d + 1);
      return;
    }
    leaves.push_back({ix, iy});
  }

  // Sample s of leaf (ix, iy), as in tree.cc
  tuple<double, double> sample(const int ix, const int iy, const int s) const {
    const double w = (X1 - X0) / double(p.base << p.depth), h = (Y1 - Y0) / double(p.base << p.depth);
    const int j = s % (p.strata * p.strata), jx = j % p.strata, jy = j / p.strata;
    const uint64_t key = mix64(uint64_t(uint32_t(ix)) | uint64_t(uint32_t(iy)) << 32) + uint64_t(s);
    return {X0 + (ix + (jx + uniform(p.seed, key, 0)) / p.strata) * w,
            Y0 + (iy + (jy + uniform(p.seed, key, 1)) / p.strata) * h};
  }

  // Doubled area estimate and variance for threshold k, directly from per-leaf group means
  tuple<double, double> area(const int k) const {
    const int K = p.ks.size(), ss = p.strata * p.strata, G = p.m / ss;
    double est = 0, var = 0;
    const double a = (X1 - X0) / double(p.base << p.depth) * ((Y1 - Y0) / double(p.base << p.depth));
    for (int d = 0; d <= p.depth; d++)
      est += (X1 - X0) / double(p.base << d) * ((Y1 - Y0) / double(p.base << d)) * certified[d * K + k];
    for (const auto& [ix, iy] : leaves) {
      double s1 = 0, s2 = 0;
      for (int g = 0; g < G; g++) {
        double v = 0;
        for (int t = 0; t < ss; t++) {
          const auto [x, y] = sample(ix, iy, g * ss + t);
          v += below(escape(x, y, p.max_iter), p.ks[k]);
        }
        v /= ss;
        s1 += v; s2 += v * v;
      }
      const double mean = s1 / G;
      est += a * mean;
      var += a * a * (s2 - s1 * mean) / (G - 1) / G;
    }
    return {2 * est, 4 * var};
  }
};

TEST(scramble) {
  for (const int64_t n : {1, 2, 3, 16, 1000, 1024, 65537, 1 << 20}) {
    const int64_t s = scramble_stride(n);
    vector<bool> seen(n);
    for (int64_t j = 0; j < n; j++) {
      const int64_t i = scramble(j, s, n);
      ASSERT_TRUE(0 <= i && i < n && !seen[i]) << tfm::format("n %d, stride %d, j %d -> %d", n, s, j, i);
      seen[i] = true;
    }
    if (n >= 1000) ASSERT_LE(n / 4, std::min(s, n - s)) << tfm::format("n %d, stride %d", n, s);
  }
}

// Visits each item once: engine claims must cover [0, n) exactly, including n > 2^30 on the GPU
struct VisitTask {
  struct State { int64_t unused; };
  int64_t burst = 1;
  int min_blocks = 3;
  uint8_t* visits;
  __host__ __device__ bool start(State&, const int64_t) const { return true; }
  __host__ __device__ bool run(State&) const { return true; }
  __host__ __device__ int64_t iters(const State&) const { return 1; }
  __host__ __device__ int64_t progress(const State&) const { return 0; }
  __host__ __device__ void finish(const State&, const int64_t i) const { visits[i]++; }
};

struct CountBad {
  const uint8_t* visits;
  int64_t n, chunk;
  int64_t* out;
  __host__ __device__ void operator()(const int64_t c) const {
    int64_t bad = 0;
    for (int64_t i = c * chunk; i < n && i < (c + 1) * chunk; i++) bad += visits[i] != 1;
    out[c] = bad;
  }
};

void check_claims(const int64_t n, const bool cuda) {
  Mem<uint8_t> visits(n, cuda);
  visits.zero();
  const auto stats = run_orbits(VisitTask{1, 3, visits.p}, n, cuda);
  ASSERT_EQ(stats.iters, n);
  const int64_t chunk = 1 << 20, chunks = (n + chunk - 1) / chunk;
  Mem<int64_t> out(chunks, cuda);
  for_each(chunks, CountBad{visits.p, n, chunk, out.p}, cuda);
  vector<int64_t> h(chunks);
  out.to_host(h.data(), chunks);
  int64_t bad = 0;
  for (const auto b : h) bad += b;
  ASSERT_EQ(bad, 0) << tfm::format("n %d, %s", n, cuda ? "cuda" : "cpu");
}

TEST(engine_claims) {
  for (const int64_t n : {1, 17, 1000, 1 << 20}) check_claims(n, false);
  IF_CUDA(for (const int64_t n : {int64_t(1), int64_t(17), int64_t(1000), int64_t(1) << 20, int64_t(1200) << 20})
            check_claims(n, true);)
}

TEST(tree_matches_reference) {
  const auto p = small_params();
  const auto R = run_tree(p);
  const Reference ref(p);
  ASSERT_EQ(R.exact, ref.exact);
  ASSERT_EQ(R.certified, ref.certified);
  ASSERT_EQ(R.leaves, int64_t(ref.leaves.size()));
  ASSERT_TRUE(R.leaves > 100) << R.leaves;
  for (int k = 0; k < int(p.ks.size()); k++) {
    const auto [est, var] = ref.area(k);
    const double e = R.area_estimate(k), v = R.variance(R.area, k);
    ASSERT_TRUE(std::abs(e - est) <= 1e-11 && std::abs(v - var) <= 1e-11 * var)  // Summation order differs
        << tfm::format("k %d: estimate %.17g vs %.17g, variance %.17g vs %.17g", p.ks[k], e, est, v, var);
    ASSERT_TRUE(v > 0) << v;
  }
  // Consecutive differences agree with differences of areas
  for (int k = 0; k + 1 < int(p.ks.size()); k++) {
    const double d = R.diff_estimate(k), dd = R.area_estimate(k) - R.area_estimate(k + 1);
    ASSERT_TRUE(std::abs(d - dd) <= 1e-13) << tfm::format("k %d: %.17g vs %.17g", k, d, dd);
  }
}

TEST(batching_invariant) {
  auto p = small_params();
  p.batch = 1 << 30;
  const auto a = run_tree(p);
  p.batch = 50;  // Many small batches
  const auto b = run_tree(p);
  ASSERT_TRUE(b.batches > 5) << b.batches;
  ASSERT_EQ(a.certified, b.certified);
  ASSERT_EQ(a.leaves, b.leaves);
  ASSERT_EQ(a.leaf_iters, b.leaf_iters);
  for (int k = 0; k < int(p.ks.size()); k++) {
    ASSERT_EQ(a.area[k].s, b.area[k].s);
    ASSERT_EQ(a.area[k].q, b.area[k].q);
    ASSERT_EQ(a.area[k].p, b.area[k].p);
  }
}

TEST(compare_float) {
  auto p = small_params();
  p.prec = "compare";
  const auto R = run_tree(p);
  ASSERT_TRUE(R.flips > 0 && R.flips < R.leaves * p.m / 10) << tfm::format("%d flips of %d", R.flips, R.leaves * p.m);
  for (int k = 0; k < int(p.ks.size()); k++) {
    const double d = R.estimate(R.delta, k, false), f = R.estimate(R.float_area, k, true);
    ASSERT_TRUE(std::abs(f - R.area_estimate(k) - d) <= 1e-13) << tfm::format("k %d", k);
    ASSERT_TRUE(std::abs(d) <= 6 * std::sqrt(R.variance(R.delta, k)) + 1e-12) << tfm::format("k %d: delta %g", k, d);
  }
}

TEST(cuda_matches_cpu) {
  IF_CUDA({
    auto p = small_params();
    p.max_iter = 1 << 16;
    p.ks = {16, 256, 4096, 65000};
    for (const string prec : {"double", "float"}) {
      p.prec = prec;
      p.cuda = false;
      const auto c = run_tree(p);
      p.cuda = true;
      const auto g = run_tree(p);
      int64_t cert_diff = 0;
      for (size_t i = 0; i < c.certified.size(); i++) cert_diff += std::abs(c.certified[i] - g.certified[i]);
      print("cuda vs cpu (%s): leaves %d vs %d, certified differences %d, leaf iterations %d vs %d", prec,
            g.leaves, c.leaves, cert_diff, g.leaf_iters, c.leaf_iters);
      ASSERT_LE(std::abs(c.leaves - g.leaves), c.leaves / 10000 + 1);
      ASSERT_LE(cert_diff, c.leaves / 10000 + 1);
      for (int k = 0; k < int(p.ks.size()); k++) {
        const double ec = c.area_estimate(k), eg = g.area_estimate(k), sd = std::sqrt(c.variance(c.area, k));
        print("  k %d: cpu %.12f, cuda %.12f, std err %.2e", p.ks[k], ec, eg, sd);
        ASSERT_LE(std::abs(ec - eg), 0.1 * sd + 1e-12);
      }
    }
  })
}

}  // namespace
}  // namespace mandelbrot
