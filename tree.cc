// Certified adaptive quadtree Monte Carlo, on CPU threads or the GPU (see tree.h)

#include "tree.h"
#include "debug.h"
#include "engine.h"
#include "orbit.h"
#include <cmath>
namespace mandelbrot {
namespace {

// Domain: the upper half plane part of the box, doubled at the end by symmetry
const double X0 = -2, X1 = 0.5, Y0 = 0, Y1 = 1.2;

// Cell (ix, iy) at depth d is [X0 + ix w_d, X0 + (ix + 1) w_d] × [Y0 + iy h_d, ...]
struct Cell { int32_t ix, iy; };

// Cells of one level: explicit, or implicitly the base grid rows [row0, row0 + n / base)
struct Level {
  const Cell* cells;
  int64_t base, row0;
  __host__ __device__ Cell at(const int64_t i) const {
    // Level sizes are below 2^31, so 32-bit division suffices
    const uint32_t i32 = uint32_t(i), b = uint32_t(base);
    return cells ? cells[i] : Cell{int32_t(i32 % b), int32_t(row0 + i32 / b)};
  }
};

// Center classification.  status[i] = kUncertified, or the below-threshold mask of a certified cell.
const uint32_t kUncertified = 0xffffffff;

struct CenterTask {
  typedef OrbitDE State;
  int64_t burst;  // Steps per run call
  Level level;
  double w, h, r;  // Cell size and half-diagonal at this depth
  double safety;
  int64_t max_iter, first_newton;
  int K;
  int ks[32];
  uint32_t* status;

  __host__ __device__ bool start(State& o, const int64_t i) const {
    const Cell c = level.at(i);
    return o.start(X0 + (c.ix + 0.5) * w, Y0 + (c.iy + 0.5) * h, first_newton);
  }
  __host__ __device__ bool run(State& o) const { return o.run(max_iter, burst); }
  __host__ __device__ int64_t iters(const State& o) const { return o.r.e.iters; }
  __host__ __device__ void finish(const State& o, const int64_t i) const {
    const EscapeDE& e = o.r;
    uint32_t s = kUncertified;
    if (e.dist > 0 && r * safety <= e.dist) {
      if (e.e.steps < 0) {
        s = K == 32 ? 0xffffffff : (1u << K) - 1;  // Interior: below every threshold
      } else {
        // Exterior: Harnack on the disk of radius dist gives g ∈ [lo, hi] · g0 on the cell
        const double t = r / e.dist, lo = std::log2((1 - t) / (1 + t)), hi = -lo;
        bool certain = true;
        uint32_t below = 0;
        for (int k = 0; k < K; k++) {
          const bool b = e.e.log2g + hi < -ks[k];
          certain &= b || e.e.log2g + lo >= -ks[k];
          below |= uint32_t(b) << k;
        }
        if (certain) s = below;
      }
    }
    // A certified all-below mask for K = 32 would collide with kUncertified; K ≤ 31 is enforced
    status[i] = s;
  }
};

// Per chunk of statuses: [uncertified, certified, certified below k for k < K]
const int64_t kChunk = 4096;

struct CountChunk {
  const uint32_t* status;
  int64_t n;
  int K;
  int64_t* out;  // [chunk * (K + 2) + ...]
  __host__ __device__ void operator()(const int64_t c) const {
    int64_t* o = out + c * (K + 2);
    for (int k = 0; k < K + 2; k++) o[k] = 0;
    const int64_t hi = (c + 1) * kChunk < n ? (c + 1) * kChunk : n;
    for (int64_t i = c * kChunk; i < hi; i++) {
      const uint32_t s = status[i];
      if (s == kUncertified) { o[0]++; continue; }
      o[1]++;
      for (int k = 0; k < K; k++) o[2 + k] += (s >> k) & 1;
    }
  }
};

// Write uncertified cells' children (or the cells themselves, as leaves) in index order
struct EmitChunk {
  const uint32_t* status;
  Level level;
  int64_t n;
  const int64_t* offsets;  // Per chunk: number of uncertified cells in earlier chunks
  bool children;
  Cell* out;
  __host__ __device__ void operator()(const int64_t c) const {
    int64_t j = offsets[c];
    const int64_t hi = (c + 1) * kChunk < n ? (c + 1) * kChunk : n;
    for (int64_t i = c * kChunk; i < hi; i++) {
      if (status[i] != kUncertified) continue;
      const Cell a = level.at(i);
      if (children) {
        Cell* o = out + 4 * j;
        o[0] = {2 * a.ix, 2 * a.iy};
        o[1] = {2 * a.ix + 1, 2 * a.iy};
        o[2] = {2 * a.ix, 2 * a.iy + 1};
        o[3] = {2 * a.ix + 1, 2 * a.iy + 1};
      } else {
        out[j] = a;
      }
      j++;
    }
  }
};

// Leaf samples: sample i is sample i % m of leaf i / m, placed by counter-based randomness keyed by the
// leaf's coordinates, so that it does not depend on processing order
template<class T> struct SampleTask {
  typedef Orbit<T> State;
  int64_t burst;  // Steps per run call
  const Cell* leaves;
  int m, strata;
  uint64_t seed;
  double w, h;  // Leaf size
  int64_t max_iter, first_newton;
  int max_period, newton_iters;
  int K;
  int ks[32];
  uint32_t* bits;

  __host__ __device__ bool start(State& o, const int64_t i) const {
    // Item counts are below 2^31 (scramble_stride checks), so 32-bit division suffices
    const uint32_t i32 = uint32_t(i), m32 = uint32_t(m);
    const Cell l = leaves[i32 / m32];
    const int s = int(i32 % m32), j = s % (strata * strata), jx = j % strata, jy = j / strata;
    const uint64_t key = mix64(uint64_t(uint32_t(l.ix)) | uint64_t(uint32_t(l.iy)) << 32) + uint64_t(s);
    const double x = X0 + (l.ix + (jx + uniform(seed, key, 0)) / strata) * w,
                 y = Y0 + (l.iy + (jy + uniform(seed, key, 1)) / strata) * h;
    return o.start(x, y, first_newton, max_period, false, newton_iters);
  }
  __host__ __device__ bool run(State& o) const { return o.run(max_iter, burst); }
  __host__ __device__ int64_t iters(const State& o) const { return o.e.iters; }
  __host__ __device__ void finish(const State& o, const int64_t i) const {
    uint32_t b = 0;
    for (int k = 0; k < K; k++) b |= uint32_t(o.e.steps < 0 || escaped_below(o.e.steps, o.e.r2, ks[k])) << k;
    bits[i] = b;
  }
};

// Group sums of a per-sample value over a chunk of leaves.  Values: area (bit k of a), diff (bit k minus bit
// k + 1 of a), delta (bit k of b minus bit k of a), or flips (a != b, into s only).
enum Kind { kArea, kDiff, kDelta, kFlips };
const int64_t kLeafChunk = 1024;

struct ReduceChunk {
  const uint32_t* a;
  const uint32_t* b;
  Kind kind;
  int k, m, ss;
  int64_t leaves;
  int64_t* out;  // [3 * chunk + {s, q, p}]
  __host__ __device__ int value(const int64_t i) const {
    switch (kind) {
      case kArea: return (a[i] >> k) & 1;
      case kDiff: return int((a[i] >> k) & 1) - int((a[i] >> (k + 1)) & 1);
      case kDelta: return int((b[i] >> k) & 1) - int((a[i] >> k) & 1);
      default: return a[i] != b[i];
    }
  }
  __host__ __device__ void operator()(const int64_t c) const {
    int64_t s = 0, q = 0, p = 0;
    const int64_t hi = (c + 1) * kLeafChunk < leaves ? (c + 1) * kLeafChunk : leaves;
    for (int64_t l = c * kLeafChunk; l < hi; l++) {
      int64_t total = 0;
      for (int g = 0; g < m / ss; g++) {
        int64_t cg = 0;
        for (int t = 0; t < ss; t++) cg += value(l * m + g * ss + t);
        total += cg;
        q += cg * cg;
      }
      s += total;
      p += total * total;
    }
    out[3 * c] = s; out[3 * c + 1] = q; out[3 * c + 2] = p;
  }
};

GroupSums reduce(const Mem<uint32_t>& a, const Mem<uint32_t>* b, const Kind kind, const int k,
                 const TreeParams& p, const int64_t leaves) {
  const int64_t chunks = (leaves + kLeafChunk - 1) / kLeafChunk;
  Mem<int64_t> out(3 * chunks, p.cuda);
  for_each(chunks, ReduceChunk{a.p, b ? b->p : nullptr, kind, k, p.m, p.strata * p.strata, leaves, out.p}, p.cuda);
  vector<int64_t> h(3 * chunks);
  out.to_host(h.data(), 3 * chunks);
  GroupSums g;
  for (int64_t c = 0; c < chunks; c++) g += GroupSums{h[3 * c], h[3 * c + 1], h[3 * c + 2]};
  return g;
}

template<class T> int64_t sample(const Mem<Cell>& leaves, const int64_t n_leaves, const TreeParams& p,
                                 const double w, const double h, Mem<uint32_t>& bits, int64_t& overflow) {
  SampleTask<T> task{p.burst, leaves.p, p.m, p.strata, p.seed, w, h, p.max_iter, p.first_newton, p.newton_max_period,
                     p.newton_iters, int(p.ks.size()), {}, bits.p};
  for (size_t k = 0; k < p.ks.size(); k++) task.ks[k] = p.ks[k];
  const auto stats = run_orbits(task, n_leaves * p.m, p.cuda);
  overflow += stats.overflow;
  return stats.iters;
}

double secs_since(const std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

}  // namespace

double TreeResult::cell_area(const int d) const {
  return (X1 - X0) / double(p.base << d) * ((Y1 - Y0) / double(p.base << d));
}

double TreeResult::estimate(const vector<GroupSums>& sums, const int k, const bool with_certified,
                            const int sign) const {
  // Leaves: Σ_l a · mean_l = a Σ c / m.  Certified cells: exact.  Doubled for the lower half plane.
  const int K = p.ks.size();
  double e = cell_area(p.depth) * double(sums[k].s) / p.m;
  if (with_certified)
    for (int d = 0; d <= p.depth; d++) e += sign * cell_area(d) * double(certified[d * K + k]);
  return 2 * e;
}

double TreeResult::diff_estimate(const int k) const {
  const int K = p.ks.size();
  double e = cell_area(p.depth) * double(diff[k].s) / p.m;
  for (int d = 0; d <= p.depth; d++) e += cell_area(d) * double(certified[d * K + k] - certified[d * K + k + 1]);
  return 2 * e;
}

double TreeResult::variance(const vector<GroupSums>& sums, const int k) const {
  // Per leaf with group sums c_g and G = m / ss groups: var = a^2 (Σ c_g^2 - (Σ c_g)^2 / G) / (ss^2 G (G - 1))
  const double ss = p.strata * p.strata, G = p.m / ss, a = cell_area(p.depth);
  const auto& g = sums[k];
  return 4 * a * a * (double(g.q) - double(g.p) / G) / (ss * ss * G * (G - 1));  // 4: doubled estimate
}

TreeResult run_tree(const TreeParams& p) {
  const int K = p.ks.size();
  slow_assert(0 < K && K <= 31, "need 1 to 31 thresholds, got %d", K);
  slow_assert(p.prec == "double" || p.prec == "float" || p.prec == "compare", "bad prec %s", p.prec);
  const int ss = p.strata * p.strata;
  slow_assert(p.strata >= 1 && p.m % ss == 0 && p.m / ss >= 2,
              "need m a multiple of strata^2 with at least 2 groups for variance estimates");
  slow_assert((p.base << p.depth) < (int64_t(1) << 31), "grid too fine for 32-bit cell coordinates");
  const bool compare = p.prec == "compare", single = p.prec == "float";

  TreeResult R;
  R.p = p;
  R.certified.assign((p.depth + 1) * K, 0);
  R.exact.assign(p.depth + 1, 0);
  R.area.resize(K); R.diff.resize(K); R.float_area.resize(K); R.delta.resize(K);
  const auto t0 = std::chrono::steady_clock::now();
  const int64_t rows = p.rows < 0 ? p.base : std::min(p.rows, p.base);

  int64_t rows_per_batch = std::max<int64_t>(1, std::min<int64_t>(rows, 16));
  for (int64_t row0 = 0; row0 < rows;) {
    const int64_t row1 = std::min(rows, row0 + rows_per_batch);
    R.batches++;

    // Tree levels
    const auto t1 = std::chrono::steady_clock::now();
    int64_t n = (row1 - row0) * p.base;
    Mem<Cell> cells(0, p.cuda), leaves(0, p.cuda);
    int64_t n_leaves = 0;
    for (int d = 0; d <= p.depth; d++) {
      const Level level{d ? cells.p : nullptr, p.base, row0};
      const double w = (X1 - X0) / double(p.base << d), h = (Y1 - Y0) / double(p.base << d);
      Mem<uint32_t> status(n, p.cuda);
      CenterTask task{p.burst, level, w, h, 0.5 * std::hypot(w, h), p.safety, std::min(p.max_iter, p.center_max_iter),
                      p.center_first_newton, K, {},
                      status.p};
      for (int k = 0; k < K; k++) task.ks[k] = p.ks[k];
      const auto stats = run_orbits(task, n, p.cuda);
      R.center_kernel_secs += stats.secs;
      R.centers += n;
      R.center_iters += stats.iters;
      R.overflow += stats.overflow;

      // Count per chunk, then emit the next level in index order
      const int64_t chunks = (n + kChunk - 1) / kChunk;
      Mem<int64_t> counts(chunks * (K + 2), p.cuda);
      for_each(chunks, CountChunk{status.p, n, K, counts.p}, p.cuda);
      vector<int64_t> hc(chunks * (K + 2)), offsets(chunks);
      counts.to_host(hc.data(), hc.size());
      int64_t uncertified = 0;
      for (int64_t c = 0; c < chunks; c++) {
        const int64_t* o = hc.data() + c * (K + 2);
        offsets[c] = uncertified;
        uncertified += o[0];
        R.exact[d] += o[1];
        for (int k = 0; k < K; k++) R.certified[d * K + k] += o[2 + k];
      }
      Mem<int64_t> doffsets(chunks, p.cuda);
      doffsets.from_host(offsets.data(), chunks);
      const bool last = d == p.depth;
      Mem<Cell> next(last ? uncertified : 4 * uncertified, p.cuda);
      for_each(chunks, EmitChunk{status.p, level, n, doffsets.p, !last, next.p}, p.cuda);
      if (last) {
        std::swap(leaves.p, next.p); std::swap(leaves.n, next.n);
        n_leaves = uncertified;
      } else {
        std::swap(cells.p, next.p); std::swap(cells.n, next.n);
        n = 4 * uncertified;
      }
    }
    R.tree_secs += secs_since(t1);

    // Leaf samples
    const auto t2 = std::chrono::steady_clock::now();
    const double w = (X1 - X0) / double(p.base << p.depth), h = (Y1 - Y0) / double(p.base << p.depth);
    Mem<uint32_t> bits(n_leaves * p.m, p.cuda), fbits(compare ? n_leaves * p.m : 0, p.cuda);
    R.leaf_iters += single ? sample<float>(leaves, n_leaves, p, w, h, bits, R.overflow)
                           : sample<double>(leaves, n_leaves, p, w, h, bits, R.overflow);
    if (compare) R.leaf_iters += sample<float>(leaves, n_leaves, p, w, h, fbits, R.overflow);
    R.sample_secs += secs_since(t2);

    // Reductions
    const auto t3 = std::chrono::steady_clock::now();
    for (int k = 0; k < K; k++) {
      R.area[k] += reduce(bits, nullptr, kArea, k, p, n_leaves);
      if (k + 1 < K) R.diff[k] += reduce(bits, nullptr, kDiff, k, p, n_leaves);
      if (compare) {
        R.float_area[k] += reduce(fbits, nullptr, kArea, k, p, n_leaves);
        R.delta[k] += reduce(bits, &fbits, kDelta, k, p, n_leaves);
      }
    }
    if (compare) R.flips += reduce(bits, &fbits, kFlips, 0, p, n_leaves).s;
    R.reduce_secs += secs_since(t3);
    R.leaves += n_leaves;

    // Aim for about p.batch leaves per batch
    const double per_row = double(n_leaves) / double(row1 - row0);
    rows_per_batch = std::max<int64_t>(1, int64_t(double(p.batch) / std::max(1.0, per_row)));
    row0 = row1;
  }
  R.secs = secs_since(t0);
  return R;
}

}  // namespace mandelbrot
