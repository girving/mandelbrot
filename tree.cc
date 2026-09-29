// Certified adaptive quadtree Monte Carlo, on CPU threads or the GPU (see tree.h)

#include "tree.h"
#include "debug.h"
#include "engine.h"
#include "orbit.h"
#include "rounded.h"
#include "double_double.h"
#include <cmath>
namespace mandelbrot {
namespace {

// Domain: the box [p.x0, p.x1] × [p.y0, p.y1] (by default the upper half plane part of [-2, 0.5] × [-1.2, 1.2]),
// doubled at the end by conjugate symmetry.  Cell (ix, iy) at depth d is [x0 + ix w_d, x0 + (ix + 1) w_d] ×
// [y0 + iy h_d, ...]
struct Cell { int32_t ix, iy; };

// Cells of one level: explicit, or implicitly base grid cells [cell0, cell0 + n) in row-major order
struct Level {
  const Cell* cells;
  int64_t base, cell0;
  __host__ __device__ Cell at(const int64_t i) const {
    if (cells) return cells[i];
    const int64_t j = cell0 + i;
    return Cell{int32_t(j % base), int32_t(j / base)};
  }
};

// Center classification.  status[i] = kUncertified, or the below-threshold mask of a certified cell.
const uint32_t kUncertified = 0xffffffff;

struct CenterTask {
  typedef OrbitDE State;
  int64_t burst;  // Steps per run call
  int min_blocks;
  Level level;
  double x0, y0;   // Domain corner
  double w, h, r;  // Cell size and half-diagonal at this depth
  double safety;
  int64_t max_iter, first_newton;
  int max_period;
  int K;
  int64_t ks[32];
  uint32_t* status;

  __host__ __device__ bool start(State& o, const int64_t i) const {
    const Cell c = level.at(i);
    return o.start(x0 + (c.ix + 0.5) * w, y0 + (c.iy + 0.5) * h, first_newton, true);
  }
  // Newton, Brent's period recovery, and cardioid/disk distances are deferred, like SampleTask's Newton
  __host__ __device__ bool run(State& o) const { return o.run(max_iter, burst, true, max_period); }
  __host__ __device__ int64_t iters(const State& o) const { return o.iters(); }
  __host__ __device__ int64_t progress(const State& o) const { return o.n; }
  __host__ __device__ bool pending(const State& o) const { return o.status >= 4 && o.status <= 6; }
  __host__ __device__ bool settle(State& o) const { return o.settle(max_iter, max_period); }
  __host__ __device__ void finish(const State& o, const int64_t i) const {
    const EscapeDE e = o.result();
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
  int min_blocks;
  const Cell* leaves;
  int m, strata;
  uint64_t seed;
  double x0, y0;  // Domain corner
  double w, h;    // Leaf size
  int64_t max_iter, first_newton;
  int max_period;
  NewtonOptions newton;
  int K;
  int64_t ks[32];
  uint32_t* bits;
  uint32_t* iters_out;  // Per-sample iterations, or null

  __host__ __device__ bool start(State& o, const int64_t i) const {
    // Item counts are below 2^31 (scramble_stride checks), so 32-bit division suffices
    const uint32_t i32 = uint32_t(i), m32 = uint32_t(m);
    const Cell l = leaves[i32 / m32];
    const int s = int(i32 % m32), j = s % (strata * strata), jx = j % strata, jy = j / strata;
    const uint64_t key = mix64(uint64_t(uint32_t(l.ix)) | uint64_t(uint32_t(l.iy)) << 32) + uint64_t(s);
    const double x = x0 + (l.ix + (jx + uniform(seed, key, 0)) / strata) * w,
                 y = y0 + (l.iy + (jy + uniform(seed, key, 1)) / strata) * h;
    return o.start(x, y, first_newton);
  }
  // Newton is deferred (the GPU engine settles pending orbits together; the CPU settles them at once)
  __host__ __device__ bool run(State& o) const { return o.run(max_iter, burst, max_period); }
  __host__ __device__ int64_t iters(const State& o) const { return o.iters(); }
  __host__ __device__ int64_t progress(const State& o) const { return o.n; }
  __host__ __device__ bool pending(const State& o) const { return o.pending(); }
  __host__ __device__ bool settle(State& o) const { return o.settle(max_iter, max_period, newton); }
  __host__ __device__ void finish(const State& o, const int64_t i) const {
    uint32_t b = 0;
    for (int k = 0; k < K; k++) b |= uint32_t(o.status != 1 || escaped_below(o.n, double(o.cx), ks[k])) << k;
    bits[i] = b;
    if (iters_out) {
      const int64_t n = o.iters();
      iters_out[i] = n < int64_t(0xffffffff) ? uint32_t(n) : 0xffffffffu;
    }
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

// Tile-resolved accounting: tile of cell (ix, iy) at depth d, on a tiles × tiles grid of the box
__host__ __device__ inline int tile_of(const int32_t ix, const int32_t iy, const int64_t grid, const int T) {
  return int(int64_t(iy) * T / grid) * T + int(int64_t(ix) * T / grid);
}

// Certified-below area per tile and threshold over chunks of cells, in finest-cell units
const int64_t kTileChunk = 1 << 16;  // Small chunks keep many GPU threads busy; per-chunk arrays are tiles² K

struct TileCountChunk {
  const uint32_t* status;
  Level level;
  int64_t n, grid, weight;
  int T, K;
  int64_t* out;  // [chunk * T² K + tile * K + k]
  __host__ __device__ void operator()(const int64_t c) const {
    int64_t* o = out + c * T * T * K;
    for (int j = 0; j < T * T * K; j++) o[j] = 0;
    const int64_t hi = (c + 1) * kTileChunk < n ? (c + 1) * kTileChunk : n;
    for (int64_t i = c * kTileChunk; i < hi; i++) {
      const uint32_t s = status[i];
      if (s == kUncertified || !s) continue;
      const Cell a = level.at(i);
      int64_t* ot = o + tile_of(a.ix, a.iy, grid, T) * K;
      for (int k = 0; k < K; k++) ot[k] += weight * ((s >> k) & 1);
    }
  }
};

// Leaf group sums of consecutive differences per tile over chunks of leaves
struct TileReduceChunk {
  const uint32_t* bits;
  const Cell* leaves;
  int64_t n, grid;
  int T, K, m, ss;
  int64_t* out;  // [chunk * T² K 3 + (tile * K + k) * 3 + {s, q, p}]
  __host__ __device__ void operator()(const int64_t c) const {
    int64_t* o = out + c * T * T * K * 3;
    for (int j = 0; j < T * T * K * 3; j++) o[j] = 0;
    const int64_t hi = (c + 1) * kTileChunk < n ? (c + 1) * kTileChunk : n;
    for (int64_t l = c * kTileChunk; l < hi; l++) {
      int64_t* ot = o + tile_of(leaves[l].ix, leaves[l].iy, grid, T) * K * 3;
      for (int k = 0; k + 1 < K; k++) {
        int64_t total = 0, q = 0;
        for (int g = 0; g < m / ss; g++) {
          int64_t cg = 0;
          for (int t = 0; t < ss; t++) {
            const uint32_t b = bits[l * m + g * ss + t];
            cg += int((b >> k) & 1) - int((b >> (k + 1)) & 1);
          }
          total += cg;
          q += cg * cg;
        }
        ot[3 * k] += total; ot[3 * k + 1] += q; ot[3 * k + 2] += total * total;
      }
    }
  }
};

// Sum per-chunk arrays of S int64s over chunks: out[j] = Σ_c in[c S + j]
struct SumChunks {
  const int64_t* in;
  int64_t chunks, S;
  int64_t* out;
  __host__ __device__ void operator()(const int64_t j) const {
    int64_t s = 0;
    for (int64_t c = 0; c < chunks; c++) s += in[c * S + j];
    out[j] = s;
  }
};

vector<int64_t> sum_chunks(const Mem<int64_t>& in, const int64_t chunks, const int64_t S, const bool cuda) {
  Mem<int64_t> out(S, cuda);
  for_each(S, SumChunks{in.p, chunks, S, out.p}, cuda);
  vector<int64_t> h(S);
  out.to_host(h.data(), S);
  return h;
}

// Pilot-allocation statistics over chunks of leaves (see alloc_index in tree.h), with m = 16
const int64_t kStatsChunk = 16384;

struct LeafStatsChunk {
  const uint32_t* bits;
  const uint32_t* iters;
  int K;
  int64_t leaves;
  int64_t* out;  // [chunk * K * 108 + alloc_index(...)]
  __host__ __device__ void operator()(const int64_t c) const {
    int64_t* o = out + c * K * 108;
    for (int j = 0; j < K * 108; j++) o[j] = 0;
    const int64_t hi = (c + 1) * kStatsChunk < leaves ? (c + 1) * kStatsChunk : leaves;
    for (int64_t l = c * kStatsChunk; l < hi; l++) {
      const uint32_t* b = bits + 16 * l;
      const uint32_t* it = iters + 16 * l;
      int64_t w2 = 0;
      for (int t = 8; t < 16; t++) w2 += it[t];
      for (int k = 0; k < K; k++) {
        int c2 = 0;
        for (int t = 8; t < 16; t++) c2 += (b[t] >> k) & 1;
        for (int pi = 0; pi < 3; pi++) {
          const int P = 2 << pi;
          int c1 = 0;
          int64_t wp = 0;
          for (int t = 0; t < P; t++) { c1 += (b[t] >> k) & 1; wp += it[t]; }
          int64_t* q = o + alloc_index(k, pi, c1, 0);
          q[0]++; q[1] += c2 * (8 - c2); q[2] += w2; q[3] += wp;
        }
      }
    }
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

template<class T> int64_t sample(const Cell* leaves, const int64_t n_leaves, const TreeParams& p,
                                 const double w, const double h, Mem<uint32_t>& bits, int64_t& overflow,
                                 uint32_t* iters = nullptr) {
  SampleTask<T> task{p.burst, p.sample_min_blocks, leaves, p.m, p.strata, p.seed, p.x0, p.y0, w, h, p.max_iter, p.first_newton, p.newton_max_period,
                     NewtonOptions{p.newton_iters, p.newton_close2, p.newton_tol < 0 ? -1 : p.newton_tol * p.newton_tol,
                                   p.newton_margin, false},
                     int(p.ks.size()), {}, bits.p, iters};
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
  return (p.x1 - p.x0) / double(p.base << d) * ((p.y1 - p.y0) / double(p.base << d));
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

double TreeResult::tile_diff_estimate(const int tile, const int k) const {
  const int K = p.ks.size();
  const double a = cell_area(p.depth);
  return 2 * a * (double(tile_diff[tile * K + k].s) / p.m + double(tile_cert[tile * K + k] - tile_cert[tile * K + k + 1]));
}

double TreeResult::tile_diff_variance(const int tile, const int k) const {
  const int K = p.ks.size();
  const double ss = p.strata * p.strata, G = p.m / ss, a = cell_area(p.depth);
  const auto& g = tile_diff[tile * K + k];
  return 4 * a * a * (double(g.q) - double(g.p) / G) / (ss * ss * G * (G - 1));
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
  slow_assert(p.max_iter < kOrbitNever, "max_iter %d needs a -DMANDELBROT_ORBIT64 build", p.max_iter);
  slow_assert(p.prec == "double" || p.prec == "float" || p.prec == "compare" || p.prec == "compare30" ||
              p.prec == "compare36" || p.prec == "compare42" || p.prec == "compare48" || p.prec == "comparedd",
              "bad prec %s", p.prec);
  const int ss = p.strata * p.strata;
  slow_assert(p.strata >= 1 && p.m % ss == 0 && p.m / ss >= 2,
              "need m a multiple of strata^2 with at least 2 groups for variance estimates");
  slow_assert((p.base << p.depth) < (int64_t(1) << 31), "grid too fine for 32-bit cell coordinates");
  const bool compare = p.prec.starts_with("compare"), single = p.prec == "float";

  TreeResult R;
  R.p = p;
  R.certified.assign((p.depth + 1) * K, 0);
  R.exact.assign(p.depth + 1, 0);
  R.area.resize(K); R.diff.resize(K); R.float_area.resize(K); R.delta.resize(K);
  slow_assert(!p.leaf_stats || (p.m == 16 && !p.prec.starts_with("compare")), "leaf_stats needs m = 16, no compare");
  if (p.leaf_stats) R.alloc.assign(int64_t(K) * 108, 0);
  slow_assert(p.tiles >= 0 && p.tiles <= 64, "tiles must be in [0, 64]");
  if (p.tiles) { R.tile_diff.resize(int64_t(p.tiles) * p.tiles * K); R.tile_cert.assign(int64_t(p.tiles) * p.tiles * K, 0); }
  const auto t0 = std::chrono::steady_clock::now();
  const int64_t rows = p.rows < 0 ? p.base : std::min(p.rows, p.base);

  // Batches are runs of base cells in row-major order, sized adaptively for about p.batch leaves (and fewer
  // than 2^31 samples)
  const int64_t total = rows * p.base, max_leaves = std::min<int64_t>(p.batch, ((int64_t(1) << 31) - 1) / p.m);
  int64_t cells_per_batch = std::min<int64_t>(total, 64);  // A small first batch, to estimate leaves per cell
  for (int64_t cell0 = 0; cell0 < total;) {
    const int64_t cell1 = std::min(total, cell0 + cells_per_batch);
    R.batches++;

    // Tree levels
    const auto t1 = std::chrono::steady_clock::now();
    int64_t n = cell1 - cell0;
    Mem<Cell> cells(0, p.cuda), leaves(0, p.cuda);
    int64_t n_leaves = 0;
    for (int d = 0; d <= p.depth; d++) {
      slow_assert(n < (int64_t(1) << 31), "level %d of a batch has %d cells; lower --batch", d, n);
      const Level level{d ? cells.p : nullptr, p.base, cell0};
      const double w = (p.x1 - p.x0) / double(p.base << d), h = (p.y1 - p.y0) / double(p.base << d);
      Mem<uint32_t> status(n, p.cuda);
      CenterTask task{p.burst, p.center_min_blocks, level, p.x0, p.y0, w, h, 0.5 * std::hypot(w, h), p.safety, std::min(p.max_iter, p.center_max_iter),
                      p.center_first_newton, p.center_max_period, K, {},
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
      if (p.tiles) {
        const int64_t tc = (n + kTileChunk - 1) / kTileChunk, S = int64_t(p.tiles) * p.tiles * K;
        Mem<int64_t> tout(tc * S, p.cuda);
        for_each(tc, TileCountChunk{status.p, level, n, p.base << d, int64_t(1) << (2 * (p.depth - d)), p.tiles, K,
                                    tout.p}, p.cuda);
        const auto h = sum_chunks(tout, tc, S, p.cuda);
        for (int64_t j = 0; j < S; j++) R.tile_cert[j] += h[j];
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

    // Leaf samples and reductions, in sub-batches of fewer than 2^31 samples
    const double w = (p.x1 - p.x0) / double(p.base << p.depth), h = (p.y1 - p.y0) / double(p.base << p.depth);
    for (int64_t l0 = 0; l0 < n_leaves; l0 += max_leaves) {
      const int64_t nl = std::min(max_leaves, n_leaves - l0);
      const auto t2 = std::chrono::steady_clock::now();
      Mem<uint32_t> bits(nl * p.m, p.cuda), fbits(compare ? nl * p.m : 0, p.cuda),
                    iters(p.leaf_stats ? nl * p.m : 0, p.cuda);
      uint32_t* ip = p.leaf_stats ? iters.p : nullptr;
      R.leaf_iters += single ? sample<float>(leaves.p + l0, nl, p, w, h, bits, R.overflow, ip)
                             : sample<double>(leaves.p + l0, nl, p, w, h, bits, R.overflow, ip);
      if (compare) {
        // The alternative precision: float, or double rounded to fewer bits
        const Cell* lp = leaves.p + l0;
        R.leaf_iters += p.prec == "compare30" ? sample<Rounded<30>>(lp, nl, p, w, h, fbits, R.overflow)
                      : p.prec == "compare36" ? sample<Rounded<36>>(lp, nl, p, w, h, fbits, R.overflow)
                      : p.prec == "compare42" ? sample<Rounded<42>>(lp, nl, p, w, h, fbits, R.overflow)
                      : p.prec == "compare48" ? sample<Rounded<48>>(lp, nl, p, w, h, fbits, R.overflow)
                      : p.prec == "comparedd" ? sample<DoubleDouble>(lp, nl, p, w, h, fbits, R.overflow)
                                              : sample<float>(lp, nl, p, w, h, fbits, R.overflow);
      }
      R.sample_secs += secs_since(t2);

      const auto t3 = std::chrono::steady_clock::now();
      for (int k = 0; k < K; k++) {
        R.area[k] += reduce(bits, nullptr, kArea, k, p, nl);
        if (k + 1 < K) R.diff[k] += reduce(bits, nullptr, kDiff, k, p, nl);
        if (compare) {
          R.float_area[k] += reduce(fbits, nullptr, kArea, k, p, nl);
          R.delta[k] += reduce(bits, &fbits, kDelta, k, p, nl);
        }
      }
      if (compare) R.flips += reduce(bits, &fbits, kFlips, 0, p, nl).s;
      if (p.tiles) {
        const int64_t tc = (nl + kTileChunk - 1) / kTileChunk, S = int64_t(p.tiles) * p.tiles * K * 3;
        Mem<int64_t> tout(tc * S, p.cuda);
        for_each(tc, TileReduceChunk{bits.p, leaves.p + l0, nl, p.base << p.depth, p.tiles, K, p.m, p.strata * p.strata,
                                     tout.p}, p.cuda);
        const auto h = sum_chunks(tout, tc, S, p.cuda);
        for (int64_t j = 0; j < S / 3; j++) R.tile_diff[j] += GroupSums{h[3 * j], h[3 * j + 1], h[3 * j + 2]};
      }
      if (p.leaf_stats) {
        const int64_t chunks = (nl + kStatsChunk - 1) / kStatsChunk, size = int64_t(K) * 108;
        Mem<int64_t> out(chunks * size, p.cuda);
        for_each(chunks, LeafStatsChunk{bits.p, iters.p, K, nl, out.p}, p.cuda);
        vector<int64_t> h(chunks * size);
        out.to_host(h.data(), h.size());
        for (int64_t c = 0; c < chunks; c++)
          for (int64_t j = 0; j < size; j++) R.alloc[j] += h[c * size + j];
      }
      R.reduce_secs += secs_since(t3);
    }
    R.leaves += n_leaves;

    // Aim for about p.batch leaves per batch
    // Leaves per base cell so far (a running average, since single batches can be tiny)
    const double per_cell = double(R.leaves) / double(cell1);
    // At most 4× the last batch, since sparse early cells (common in --box domains) underestimate the density
    cells_per_batch = std::max<int64_t>(1, std::min<int64_t>(4 * (cell1 - cell0),
                                                             int64_t(double(max_leaves) / std::max(1e-3, per_cell))));
    cell0 = cell1;
  }
  R.secs = secs_since(t0);
  return R;
}

}  // namespace mandelbrot
