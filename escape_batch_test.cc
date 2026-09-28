// Escape batch tests

#include "escape_batch.h"
#include "escape.h"
#include "tests.h"
#include <random>
#include <vector>
namespace mandelbrot {
namespace {

using std::vector;

// Leaves of size h around random points that take a while to classify, so they straddle the boundary
vector<Leaf> boundary_leaves(const int count, const double h, const uint64_t seed) {
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> ux(-2, 0.5), uy(0, 1.2);
  vector<Leaf> leaves;
  while (int(leaves.size()) < count) {
    const double x = ux(rng), y = uy(rng);
    if (escape(x, y, 1 << 12).iters >= 64) leaves.push_back({x - h / 2, y - h / 2, h, h});
  }
  return leaves;
}

SampleParams params(const int m, const int strata, const int64_t max_iter) {
  SampleParams p{m, strata, 17, max_iter, 4, {}};
  p.ks[0] = 32; p.ks[1] = 128; p.ks[2] = 1024; p.ks[3] = int(max_iter - 8);
  return p;
}

TEST(sample_point) {
  const auto p = params(8, 2, 1 << 12);
  const Leaf l{-1, 0.25, 0.5, 0.25};
  for (int64_t i = 0; i < 4000; i++) {
    double x, y;
    sample_point(l, p, i, x, y);
    // Sample i lies in stratum (i % m) % 4 of the 2 × 2 grid
    const int j = int(i % p.m) % 4;
    const double x0 = l.x + (j % 2) * l.w / 2, y0 = l.y + (j / 2) * l.h / 2;
    ASSERT_TRUE(x0 <= x && x < x0 + l.w / 2 && y0 <= y && y < y0 + l.h / 2) << tfm::format("i %d, x %g, y %g", i, x, y);
    double x2, y2;
    sample_point(l, p, i, x2, y2);
    ASSERT_EQ(x, x2);
    ASSERT_EQ(y, y2);
  }
}

TEST(scramble) {
  for (const int64_t n : {1, 2, 3, 16, 1000, 1024, 65537, 1 << 20}) {
    const int64_t s = scramble_stride(n);
    vector<bool> seen(n);
    for (int64_t j = 0; j < n; j++) {
      const int64_t i = scramble(j, s, n);
      ASSERT_TRUE(0 <= i && i < n && !seen[i]) << tfm::format("n %d, stride %d, j %d -> %d", n, s, j, i);
      seen[i] = true;
    }
    // Consecutive claims should be far apart
    if (n >= 1000) ASSERT_LE(n / 4, std::min(s, n - s)) << tfm::format("n %d, stride %d", n, s);
  }
}

TEST(cpu_double_matches_escape) {
  const auto leaves = boundary_leaves(500, 1e-3, 3);
  const auto p = params(8, 2, 1 << 14);
  vector<uint32_t> bits(leaves.size() * p.m);
  const int64_t iters = sample_leaves_cpu<double>(leaves, p, bits);
  int64_t expected_iters = 0;
  for (size_t i = 0; i < bits.size(); i++) {
    double x, y;
    sample_point(leaves[i / p.m], p, i, x, y);
    const auto e = escape(x, y, p.max_iter);
    expected_iters += e.iters;
    ASSERT_EQ(bits[i], below_bits(e, p)) << tfm::format("sample %d at %.17g %.17g", i, x, y);
  }
  ASSERT_EQ(iters, expected_iters);
  // Leaves were chosen with g ≈ 2^-64 or less, so thresholds 2^-128 and 2^-1024 should split them
  for (int k = 1; k <= 2; k++) {
    int64_t n = 0;
    for (const auto b : bits) n += (b >> k) & 1;
    ASSERT_TRUE(0 < n && n < int64_t(bits.size())) << tfm::format("k %d: %d of %d below", p.ks[k], n, bits.size());
  }
}

TEST(cpu_float_mostly_matches_double) {
  const auto leaves = boundary_leaves(2000, 1e-3, 5);
  const auto p = params(8, 2, 1 << 14);
  vector<uint32_t> bd(leaves.size() * p.m), bf(bd.size());
  sample_leaves_cpu<double>(leaves, p, bd);
  sample_leaves_cpu<float>(leaves, p, bf);
  int64_t flips = 0;
  for (size_t i = 0; i < bd.size(); i++) flips += bd[i] != bf[i];
  // These leaves straddle the boundary, so flips are far more common than for uniform samples
  ASSERT_LE(flips, int64_t(bd.size()) / 20) << tfm::format("%d flips of %d", flips, bd.size());
}

TEST(cuda_matches_cpu) {
  IF_CUDA({
    const auto leaves = boundary_leaves(2000, 1e-3, 7);
    const auto p = params(16, 2, 1 << 16);
    vector<uint32_t> cpu(leaves.size() * p.m), gpu(cpu.size());
    for (const bool single : {false, true}) {
      const int64_t ci = single ? sample_leaves_cpu<float>(leaves, p, cpu) : sample_leaves_cpu<double>(leaves, p, cpu);
      const int64_t gi = single ? sample_leaves_cuda<float>(leaves, p, gpu) : sample_leaves_cuda<double>(leaves, p, gpu);
      int64_t flips = 0;
      for (size_t i = 0; i < cpu.size(); i++) flips += cpu[i] != gpu[i];
      print("cuda vs cpu (%s): %d flips of %d, iterations %d vs %d", single ? "float" : "double", flips, cpu.size(),
            gi, ci);
      ASSERT_LE(flips, int64_t(cpu.size()) / (single ? 1000 : 10000));
    }
  })
}

}  // namespace
}  // namespace mandelbrot
