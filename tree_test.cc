// Tree pipeline tests

#include "tree.h"
#include "rounded.h"
#include "orbit_expansion.h"
#include "engine.h"
#include "escape.h"
#include "tests.h"
#include <algorithm>
#include <cstdio>
#include <unistd.h>
#include <cmath>
#include <map>
#include <random>
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
  __host__ __device__ bool pending(const State&) const { return false; }
  __host__ __device__ bool settle(State&) const { return true; }
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

TEST(tree_newton_matches_reference) {
  // Leaf orbits long enough for the Newton schedule: certificates change cost, never classification
  auto p = small_params();
  p.max_iter = 1 << 16;
  p.first_newton = 1024;
  p.ks = {16, 256, 4096, 65000};
  const auto R = run_tree(p);
  const Reference ref(p);
  for (int k = 0; k < int(p.ks.size()); k++) {
    const auto [est, var] = ref.area(k);
    const double e = R.area_estimate(k);
    ASSERT_TRUE(std::abs(e - est) <= 1e-11) << tfm::format("k %d: estimate %.17g vs %.17g", p.ks[k], e, est);
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

TEST(shards_merge_exactly) {
  // Shards run separately, saved, loaded, and merged give exactly the unsharded result
  auto p = small_params();
  p.prec = "compare";
  const auto a = run_tree(p);
  auto m = empty_result(p);
  p.shards = 3;
  for (p.shard = 0; p.shard < p.shards; p.shard++) {
    const auto path = tfm::format("/tmp/tree_test_shard_%d_%d.txt", getpid(), p.shard);
    save_result(run_tree(p), path);
    auto q = p;
    q.shard = 0; q.shards = 1;
    const auto S = load_result(path, q);
    ASSERT_EQ(S.p.shard, p.shard);
    merge(m, S);
    std::remove(path.c_str());
  }
  ASSERT_EQ(a.certified, m.certified);
  ASSERT_EQ(a.exact, m.exact);
  ASSERT_EQ(a.leaves, m.leaves);
  ASSERT_EQ(a.leaf_iters, m.leaf_iters);
  ASSERT_EQ(a.flips, m.flips);
  for (int k = 0; k < int(p.ks.size()); k++)
    for (const auto v : {&TreeResult::area, &TreeResult::diff, &TreeResult::float_area, &TreeResult::delta}) {
      ASSERT_EQ((a.*v)[k].s, (m.*v)[k].s);
      ASSERT_EQ((a.*v)[k].q, (m.*v)[k].q);
      ASSERT_EQ((a.*v)[k].p, (m.*v)[k].p);
    }
}

TEST(split_batches) {
  // Batches whose tree levels reach max_level_cells split in half and retry, without changing any result
  // (2000 splits deeper levels' cell lists; 30 also splits base cell ranges)
  auto p = small_params();
  const auto a = run_tree(p);
  for (const int64_t limit : {2000, 30}) {
    p.max_level_cells = limit;
    const auto b = run_tree(p);
    ASSERT_LT(a.batches, b.batches);
    ASSERT_EQ(a.leaves, b.leaves);
    ASSERT_EQ(a.leaf_iters, b.leaf_iters);
    ASSERT_EQ(a.centers, b.centers);
    for (size_t i = 0; i < a.certified.size(); i++) ASSERT_EQ(a.certified[i], b.certified[i]);
    for (int k = 0; k < int(p.ks.size()); k++) {
      ASSERT_EQ(a.area[k].s, b.area[k].s);
      ASSERT_EQ(a.area[k].q, b.area[k].q);
      ASSERT_EQ(a.area[k].p, b.area[k].p);
    }
  }
}

TEST(split_depth) {
  // Building shallow levels first for runs of base cells, then batching the cells at split_depth, changes no
  // result (including with level splits inside the second phase)
  auto p = small_params();
  const auto a = run_tree(p);
  for (const int sd : {1, 2, 3}) {
    for (const int64_t limit : {int64_t(1) << 31, int64_t(500)}) {
      p.split_depth = sd;
      p.max_level_cells = limit;
      const auto b = run_tree(p);
      ASSERT_EQ(a.leaves, b.leaves);
      ASSERT_EQ(a.leaf_iters, b.leaf_iters);
      ASSERT_EQ(a.centers, b.centers);
      ASSERT_EQ(a.center_iters, b.center_iters);
      for (size_t i = 0; i < a.certified.size(); i++) ASSERT_EQ(a.certified[i], b.certified[i]);
      for (size_t i = 0; i < a.exact.size(); i++) ASSERT_EQ(a.exact[i], b.exact[i]);
      for (int k = 0; k < int(p.ks.size()); k++) {
        ASSERT_EQ(a.area[k].s, b.area[k].s);
        ASSERT_EQ(a.area[k].q, b.area[k].q);
        ASSERT_EQ(a.area[k].p, b.area[k].p);
      }
    }
  }
}

TEST(deep_queue) {
  // Suspending long samples into deep passes changes no result: every sum is bit-identical, plain and with roulette
  for (const bool roulette : {false, true}) {
    auto p = small_params();
    if (roulette) {
      p.first_newton = 512;
      p.max_iter = (1 << 16) + 8;
      p.ks = {1024, 2048, 4096, 8192, 16384, 32768, 65536};
      p.roulette_from = 2048;
    }
    const auto a = run_tree(p);
    p.deep_from = roulette ? 2064 : 1032;
    p.deep_batch = roulette ? 20 : 200;
    p.deep_cap = 5;  // Too little room, to exercise reruns
    const auto b = run_tree(p);
    print("  roulette %d: %d samples suspended, %d deep passes", int(roulette), b.deep_samples, b.deep_passes);
    ASSERT_LT(1, b.deep_passes);
    ASSERT_EQ(a.leaf_iters, b.leaf_iters);
    for (int k = 0; k < int(p.ks.size()); k++) {
      ASSERT_EQ(a.area[k].s, b.area[k].s);
      ASSERT_EQ(a.area[k].q, b.area[k].q);
      ASSERT_EQ(a.area[k].p, b.area[k].p);
      if (k + 1 < int(p.ks.size())) {
        ASSERT_EQ(a.diff[k].s, b.diff[k].s);
        ASSERT_EQ(a.diff[k].q, b.diff[k].q);
        ASSERT_EQ(a.diff[k].p, b.diff[k].p);
      }
    }
  }
}

TEST(roulette_unbiased) {
  // Russian roulette reweights escapes past the reference threshold without bias: estimates match the full run
  // within their extra noise (identically below the first decision)
  auto p = small_params();
  p.first_newton = 512;
  p.max_iter = (1 << 18) + 8;
  p.ks = {1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144};
  const auto a = run_tree(p);
  p.roulette_from = 2048;
  const auto b = run_tree(p);
  ASSERT_LT(b.leaf_iters, a.leaf_iters);
  ASSERT_EQ(a.area[0].s, b.area[0].s);
  ASSERT_EQ(a.area[1].s, b.area[1].s);
  for (int k = 0; k + 1 < int(p.ks.size()); k++) {
    // The runs share samples, so the difference's std error is at most the sum of theirs
    const double d = b.diff_estimate(k) - a.diff_estimate(k),
                 s = std::sqrt(a.variance(a.diff, k)) + std::sqrt(b.variance(b.diff, k));
    if (k < 1) ASSERT_EQ(d, 0);
    ASSERT_LT(std::abs(d), 4 * s + 1e-15) << tfm::format("k %d: %g vs %g", p.ks[k], d, s);
  }
}

TEST(dd_matches_compare) {
  // prec dd classifies exactly as the double-double half of comparedd
  auto p = small_params();
  p.prec = "comparedd";
  const auto a = run_tree(p);
  p.prec = "dd";
  const auto b = run_tree(p);
  for (int k = 0; k < int(p.ks.size()); k++) {
    ASSERT_EQ(a.float_area[k].s, b.area[k].s);
    ASSERT_EQ(a.float_area[k].q, b.area[k].q);
  }
}

TEST(box) {
  // A box inside the period-2 disk |c + 1| < 1/4 is certified interior everywhere: area exactly twice the box
  auto p = small_params();
  p.x0 = -1.1; p.x1 = -0.9; p.y0 = 0; p.y1 = 0.1;
  const auto R = run_tree(p);
  ASSERT_EQ(R.leaves, 0);
  for (int k = 0; k < int(p.ks.size()); k++)
    ASSERT_TRUE(std::abs(R.area_estimate(k) - 0.04) < 1e-12) << tfm::format("k %d: %.17g", k, R.area_estimate(k));
}

TEST(tiles) {
  // Tile-resolved differences partition the global ones exactly
  auto p = small_params();
  p.tiles = 5;
  const auto R = run_tree(p);
  const int K = p.ks.size(), T2 = p.tiles * p.tiles;
  for (int k = 0; k + 1 < K; k++) {
    GroupSums g;
    double D = 0, V = 0;
    for (int t = 0; t < T2; t++) { g += R.tile_diff[t * K + k]; D += R.tile_diff_estimate(t, k); V += R.tile_diff_variance(t, k); }
    ASSERT_EQ(g.s, R.diff[k].s);
    ASSERT_EQ(g.q, R.diff[k].q);
    ASSERT_EQ(g.p, R.diff[k].p);
    ASSERT_TRUE(std::abs(D - R.diff_estimate(k)) <= 1e-12 * std::abs(R.diff_estimate(k))) << tfm::format("k %d: %.17g vs %.17g", k, D, R.diff_estimate(k));
    ASSERT_TRUE(std::abs(V - R.variance(R.diff, k)) <= 1e-9 * R.variance(R.diff, k));
  }
  for (int k = 0; k < K; k++) {
    int64_t cert = 0, want = 0;
    for (int t = 0; t < T2; t++) cert += R.tile_cert[t * K + k];
    for (int d = 0; d <= p.depth; d++) want += R.certified[d * K + k] << (2 * (p.depth - d));
    ASSERT_EQ(cert, want);
  }
}

TEST(leaf_stats) {
  // Allocation statistics partition the leaves, and their costs add up to the leaf iterations
  auto p = small_params();
  p.m = 16;
  p.leaf_stats = true;
  const auto R = run_tree(p);
  ASSERT_TRUE(p.m == 16 && R.leaves > 100);
  for (int k = 0; k < int(p.ks.size()); k++) {
    // Samples 8..15 below k, from the area sums: Σ_l c2 over classes must match the area count minus pilot
    for (int pi = 0; pi < 3; pi++) {
      int64_t N = 0, W = 0, below = 0;
      for (int c = 0; c <= kPilots[pi]; c++) {
        N += R.alloc[alloc_index(k, pi, c, 0)];
        W += R.alloc[alloc_index(k, pi, c, 2)] + R.alloc[alloc_index(k, pi, c, 3)];
        below += c * R.alloc[alloc_index(k, pi, c, 0)];
      }
      ASSERT_EQ(N, R.leaves);
      if (pi == 2) {
        ASSERT_EQ(W, R.leaf_iters);
        ASSERT_TRUE(below <= R.area[k].s && 2 * below > R.area[k].s / 2) << tfm::format("%d of %d", below, R.area[k].s);
      }
    }
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

TEST(rounded) {
  std::mt19937_64 rng(3);
  std::uniform_real_distribution<double> u(-10, 10);
  for (int i = 0; i < 100000; i++) {
    const double a = u(rng), r = Rounded<30>::round(a);
    // At most 30 significant bits, within half an ulp
    int e;
    const double m = std::frexp(r, &e);
    ASSERT_EQ(std::ldexp(m, 30), std::round(std::ldexp(m, 30))) << tfm::format("%.17g -> %.17g", a, r);
    ASSERT_LE(std::abs(r - a), std::ldexp(std::abs(a), -30)) << tfm::format("%.17g -> %.17g", a, r);
    // 24 bits agrees with float rounding
    ASSERT_EQ(Rounded<24>::round(a), double(float(a))) << tfm::format("%.17g", a);
  }
}

TEST(expansion_orbits) {
  // Ordering on Expansion<2> follows its value, and comparisons with nan are false (as the orbit's block test needs)
  typedef Expansion<2> E;
  const E a(1.0), b(1.0, 0x1p-60, Nonoverlap()), c(1.0, -0x1p-60, Nonoverlap()), nan(double(NAN));
  ASSERT_TRUE(c < a && a < b && c < b && a <= a && b >= a && !(a < a) && !(b < a));
  ASSERT_TRUE(!(nan < a) && !(nan <= a) && !(nan > a) && !(nan >= a) && !(a < nan) && !(a <= nan));
  // Escape classifications over Expansion<2> orbits mostly agree with double
  auto p = small_params();
  p.prec = "comparedd";
  const auto R = run_tree(p);
  print("  comparedd: %d flips of %d samples", R.flips, R.leaves * p.m);
  ASSERT_LE(R.flips, R.leaves * p.m / 1000);
}

TEST(expansion_step) {
  // The fused double-double step is accurate to 2 · 2^-104 (|z|^2 + |c|) per step, against Expansion<3>
  typedef Expansion<2> E;
  typedef Expansion<3> F;
  std::mt19937_64 rng(7);
  std::uniform_real_distribution<double> u(-1, 1);
  double worst = 0;
  for (int i = 0; i < 1000000; i++) {
    const double r = 2 * std::abs(u(rng)), a = r * u(rng), b = r * u(rng), x = 2 * u(rng), y = 2 * u(rng);
    const double al = std::ldexp(u(rng), std::ilogb(a) - 53), bl = std::ldexp(u(rng), std::ilogb(b) - 53);
    E zx(a, al, nonoverlap), zy(b, bl, nonoverlap), zy2 = zy * zy, r2 = zx * zx + zy2;
    const F X(a, al, 0.0, nonoverlap), Y(b, bl, 0.0, nonoverlap), ex = X * X - Y * Y + x, ey = 2 * (X * Y) + y;
    orbit_step(zx, zy, zy2, r2, x, y);
    const F dx = F(zx.x[0], zx.x[1], 0.0, nonoverlap) - ex, dy = F(zy.x[0], zy.x[1], 0.0, nonoverlap) - ey;
    const double scale = std::ldexp(a * a + b * b + std::hypot(x, y), -104);
    worst = std::max(worst, std::max(std::abs(double(dx)), std::abs(double(dy))) / scale);
    ASSERT_LE(std::abs(r2.x[0] - (zx.x[0] * zx.x[0] + zy.x[0] * zy.x[0])), 1e-15 * r2.x[0]);
  }
  print("  worst step error %.3g · 2^-104 (|z|^2 + |c|)", worst);
  ASSERT_LE(worst, 2);
}

TEST(expansion3_step) {
  // The fused triple-double step is accurate to 2^-156 (|z|^2 + |c|) per step, against Expansion<4>, including
  // when the result cancels; and long orbits at an attracting parameter stay that close
  typedef Expansion<3> E;
  typedef Expansion<4> F;
  const auto up = [](const E e) { return F(e.x[0], e.x[1], e.x[2], 0.0, nonoverlap); };
  std::mt19937_64 rng(7);
  std::uniform_real_distribution<double> u(-1, 1);
  double worst[2] = {0, 0};
  int64_t overlaps[2] = {0, 0};
  const int n = 1000000;
  for (int i = 0; i < 2 * n; i++) {
    const bool cancel = i >= n;  // Second half: c ≈ -z^2, so z' is tiny
    const double r = 2 * std::abs(u(rng)), a = r * u(rng), b = r * u(rng);
    const double a1 = std::ldexp(u(rng), std::ilogb(a) - 53), b1 = std::ldexp(u(rng), std::ilogb(b) - 53);
    const double a2 = std::ldexp(u(rng), std::ilogb(a1) - 53), b2 = std::ldexp(u(rng), std::ilogb(b1) - 53);
    const double x = cancel ? -(a * a - b * b) * (1 + std::ldexp(u(rng), -40 - int(rng() % 13))) : 2 * u(rng);
    const double y = cancel ? -(2 * a * b) * (1 + std::ldexp(u(rng), -40 - int(rng() % 13))) : 2 * u(rng);
    E zx(a, a1, a2, nonoverlap), zy(b, b1, b2, nonoverlap), zy2 = zy * zy, r2 = zx * zx + zy2;
    const F X = up(zx), Y = up(zy), ex = X * X - Y * Y + x, ey = 2 * (X * Y) + y;
    orbit_step(zx, zy, zy2, r2, x, y);
    const F dx = up(zx) - ex, dy = up(zy) - ey;
    const double scale = std::ldexp(a * a + b * b + std::hypot(x, y), -155);
    // Whole sums: a difference of expansions that round differently is not normalized
    worst[cancel] = std::max(worst[cancel], std::max(std::abs(double(dx)), std::abs(double(dy))) / scale);
    for (const E z : {zx, zy})
      overlaps[cancel] += std::abs(z.x[1]) > std::ldexp(1, std::ilogb(z.x[0]) - 53) ||
                          std::abs(z.x[2]) > std::ldexp(1, std::ilogb(z.x[1]) - 53);
    ASSERT_LE(std::abs(r2.x[0] - (zx.x[0] * zx.x[0] + zy.x[0] * zy.x[0])), 1e-15 * r2.x[0]);
  }
  print("  worst step error %.3g · 2^-155 (|z|^2 + |c|), %.3g with cancellation; overlapping results %d, %d",
        worst[0], worst[1], overlaps[0], overlaps[1]);
  ASSERT_LE(std::max(worst[0], worst[1]), 0.5);
  ASSERT_EQ(overlaps[0], 0);

  // 10^4 steps at c = 0.2 + 0.5i, in the main cardioid, from z = c: contracting, so errors stay at the step size
  const double x = 0.2, y = 0.5;
  E zx(x), zy(y), zy2 = zy * zy, r2 = zx * zx + zy2;
  F X(x), Y(y);
  for (int i = 0; i < 10000; i++) {
    orbit_step(zx, zy, zy2, r2, x, y);
    const F t = X * X - Y * Y + x;
    Y = 2 * (X * Y) + y;
    X = t;
  }
  const double err = std::max(std::abs(double(up(zx) - X)), std::abs(double(up(zy) - Y)));
  print("  10^4 steps: error %.3g · 2^-155", std::ldexp(err, 155));
  ASSERT_LE(err, std::ldexp(1, -150));
}

TEST(expansion3_orbits) {
  // Orbit<Expansion<3>> (generic Newton and cycle code around the fused step) classifies like double and
  // Expansion<2>, apart from rare precision-sensitive samples
  std::mt19937_64 rng(11);
  std::uniform_real_distribution<double> ux(-2, 0.5), uy(0, 1.2);
  const int n = 4000;
  int differ2 = 0, differ1 = 0, escaped = 0, interior = 0;
  for (int i = 0; i < n; i++) {
    const double x = ux(rng), y = uy(rng);
    Orbit<double> o1;
    Orbit<Expansion<2>> o2;
    Orbit<Expansion<3>> o3;
    o1.start(x, y); o2.start(x, y); o3.start(x, y);
    o1.finish(1 << 14); o2.finish(1 << 14); o3.finish(1 << 14);
    const Escape e1 = o1.result(), e2 = o2.result(), e3 = o3.result();
    differ1 += e1.steps != e3.steps;
    differ2 += e2.steps != e3.steps;
    escaped += e3.steps >= 0;
    interior += e3.period > 0;
  }
  print("  %d samples: %d escaped, %d attracting; Expansion<3> differs from Expansion<2> on %d, from double on %d",
        n, escaped, interior, differ2, differ1);
  ASSERT_LE(n / 4, escaped);
  ASSERT_LE(n / 8, interior);
  ASSERT_LE(differ2, 2);
  ASSERT_LE(differ1, 4);
}

TEST(compare_rounded) {
  // Fewer bits flip more samples; 48 bits flips few
  auto p = small_params();
  int64_t last = -1;
  for (const string prec : {"compare48", "compare36", "compare30"}) {
    p.prec = prec;
    const auto R = run_tree(p);
    print("  %s: %d flips of %d samples", prec, R.flips, R.leaves * p.m);
    ASSERT_LE(last, R.flips);
    last = R.flips;
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
    // The deep queue on the GPU: bit-identical to the GPU without it
    p.prec = "double";
    p.cuda = true;
    const auto g0 = run_tree(p);
    p.deep_from = 1032;
    p.deep_batch = 200;
    const auto g1 = run_tree(p);
    print("cuda deep queue: %d samples suspended, %d deep passes", g1.deep_samples, g1.deep_passes);
    ASSERT_LT(1, g1.deep_passes);
    ASSERT_EQ(g0.leaf_iters, g1.leaf_iters);
    for (int k = 0; k < int(p.ks.size()); k++) {
      ASSERT_EQ(g0.area[k].s, g1.area[k].s);
      ASSERT_EQ(g0.area[k].q, g1.area[k].q);
      ASSERT_EQ(g0.area[k].p, g1.area[k].p);
    }
  })
}

}  // namespace
}  // namespace mandelbrot
