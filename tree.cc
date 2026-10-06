// Certified adaptive quadtree Monte Carlo, on CPU threads or the GPU (see tree.h)

#include "tree.h"
#include "debug.h"
#include "engine.h"
#include "orbit.h"
#include "rounded.h"
#include "orbit_expansion.h"
#include <cmath>
#include <functional>
#include <mutex>
#include <thread>
namespace mandelbrot {
namespace {

// Domain: the box [p.x0, p.x1] × [p.y0, p.y1] (by default the upper half plane part of [-2, 0.5] × [-1.2, 1.2]),
// doubled at the end by conjugate symmetry.  Cell (ix, iy) at depth d is [x0 + ix w_d, x0 + (ix + 1) w_d] ×
// [y0 + iy h_d, ...]
struct Cell { int32_t ix, iy; };

// Cells of one level: explicit, or implicitly base grid cells [cell0, cell0 + n) in row-major order over this
// shard's rows (row r of the shard is base row r · shards + shard)
struct Level {
  const Cell* cells;
  int64_t base, cell0;
  int shard, shards;
  __host__ __device__ Cell at(const int64_t i) const {
    if (cells) return cells[i];
    const int64_t j = cell0 + i;
    return Cell{int32_t(j % base), int32_t(j / base * shards + shard)};
  }
};

// Center classification.  status[i] = kUncertified, kUncertifiedInterior (the center is interior, but too
// near its component's boundary to certify the cell), or the below-threshold mask of a certified cell.
const uint32_t kUncertified = 0xffffffff, kUncertifiedInterior = 0xfffffffe;
__host__ __device__ static inline bool certified(const uint32_t s) { return s < kUncertifiedInterior; }

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
  const uint8_t* hints;  // Per cell: its parent's center was interior (or null)
  int64_t hint_newton;   // First Newton step for hinted cells (0: first_newton): they are mostly interior

  __host__ __device__ bool start(State& o, const int64_t i) const {
    const Cell c = level.at(i);
    const int64_t first = hint_newton && hints && hints[i] ? hint_newton : first_newton;
    return o.start(x0 + (c.ix + 0.5) * w, y0 + (c.iy + 0.5) * h, first, true);
  }
  // Newton, Brent's period recovery, and cardioid/disk distances are deferred, like SampleTask's Newton
  __host__ __device__ bool run(State& o) const { return o.run(max_iter, burst, true, max_period); }
  __host__ __device__ int64_t iters(const State& o) const { return o.iters(); }
  __host__ __device__ int64_t progress(const State& o) const { return o.n; }
  __host__ __device__ bool pending(const State& o) const { return o.status >= 4 && o.status <= 6; }
  __host__ __device__ bool settle(State& o) const { return o.settle(max_iter, max_period); }
  // Settle work, for sorting settles (engine.h): Newton's period, or the most for Brent's period recovery
  __host__ __device__ int settle_key(const State& o) const {
    return o.status == 4 ? o.atom_candidate(max_period) : o.status == 5 ? 511 : 0;
  }
  __host__ __device__ void finish(const State& o, const int64_t i) const { finish(o.result(), i); }

  // Compact finish input, which the GPU engine buffers so that a warp finishes 32 cells at once (the distance
  // bound and Harnack's logarithms are a few hundred instructions).  Escaped: |z|^2, dz/dc, step, and exponent.
  // Otherwise dexp = -1, with the distance, log2 g, and steps.
  struct Record { double a, b, c; int32_t n, dexp; };
  __host__ __device__ Record record(const State& o) const {
    if (o.status == 7) return {o.cx, o.dx, o.dy, int32_t(o.n), o.dexp};
    return {o.r.dist, o.r.e.log2g, 0, int32_t(o.r.e.steps), -1};
  }
  __host__ __device__ void finish_record(const Record& c, const int64_t i) const {
    if (c.dexp >= 0) {
      finish(escaped(c.n, c.a, c.b, c.c, c.dexp), i);
    } else {
      EscapeDE e;
      e.e.steps = c.n;
      e.e.log2g = c.b;
      e.dist = c.a;
      finish(e, i);
    }
  }

  __host__ __device__ void finish(const EscapeDE& e, const int64_t i) const {
    uint32_t s = e.e.steps < 0 && e.dist > 0 ? kUncertifiedInterior : kUncertified;
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
      if (!certified(s)) { o[0]++; continue; }
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
  uint8_t* hints;  // With children: whether each child's parent center was interior
  __host__ __device__ void operator()(const int64_t c) const {
    int64_t j = offsets[c];
    const int64_t hi = (c + 1) * kChunk < n ? (c + 1) * kChunk : n;
    for (int64_t i = c * kChunk; i < hi; i++) {
      if (certified(status[i])) continue;
      const Cell a = level.at(i);
      if (children) {
        Cell* o = out + 4 * j;
        o[0] = {2 * a.ix, 2 * a.iy};
        o[1] = {2 * a.ix + 1, 2 * a.iy};
        o[2] = {2 * a.ix, 2 * a.iy + 1};
        o[3] = {2 * a.ix + 1, 2 * a.iy + 1};
        const uint8_t hint = status[i] == kUncertifiedInterior;
        for (int q = 0; q < 4; q++) hints[4 * j + q] = hint;
      } else {
        out[j] = a;
      }
      j++;
    }
  }
};

// Roulette weight byte flag: the sample was killed (its exponent is in the low bits)
constexpr uint8_t kKilled = 0x80;

// A sample's outcome for --flip-stats: kind (0 escaped at step n, 1 certified interior, 2 hit max_iter) in the top
// two bits, and n / 4 in the rest (steps below 2^32)
__host__ __device__ static inline uint32_t outcome_code(const int kind, const int64_t n) {
  const int64_t q = n / 4 < (int64_t(1) << 30) - 1 ? n / 4 : (int64_t(1) << 30) - 1;
  return uint32_t(kind) << 30 | uint32_t(q);
}

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
  uint32_t* outcome_out;  // Per-sample outcome (outcome_code), or null
  int64_t park_steps;   // GPU steps before parking: the first Newton step (see run_orbits)
  // Russian roulette (see TreeParams::roulette_from): at Newton steps d ≥ roulette_from, an orbit Newton did not
  // certify goes on with probability 2^-roulette_log2 and stops otherwise (status 9, killed); rexp receives each
  // sample's weight exponent, the roulette decisions it survived (see roulette_exponent), with kKilled set if
  // roulette stopped it
  int64_t roulette_from;
  int roulette_log2, roulette_stride;
  uint8_t* rexp;
  // The deep queue (TreeParams::deep_from): orbits still unsettled at step suspend_at stop there (status 10) and
  // their states go to deep_states[k], deep_items[k] = i, for k from *deep_count; a deep pass then resumes them as
  // items of their own (resume[i]), with iters_offset steps already counted
  int64_t suspend_at;
  State* deep_states;
  int64_t* deep_items;
  uint64_t* deep_count;
  int64_t deep_cap;
  const State* resume;
  int64_t iters_offset;

  // Whether Newton step d (one of the d_j below) is a roulette decision: every roulette_stride-th from the first at
  // or after roulette_from
  __host__ __device__ bool decision(const int64_t d) const {
    if (!roulette_log2 || d < roulette_from) return false;
    int j = 0;
    for (int64_t e = first_newton + 8; e < d; e = 2 * (e - 8) + 8) j += e >= roulette_from;
    return j % roulette_stride == 0;
  }

  // Roulette decisions happen at the Newton steps d_j = first_newton · 2^j + 8 (the first multiple of 8 past each
  // doubling) with d_j ≥ roulette_from.  A sample's exponent counts those it survived, which follows from where it
  // ended: every decision before its last step (and at it, for escapes: an orbit escaping at a Newton step faced
  // the decision there first).
  __host__ __device__ int roulette_exponent(const State& o) const {
    if (!roulette_log2) return 0;
    const bool escaped = o.status == 1 || o.status == 8;
    int e = 0;
    for (int64_t d = first_newton + 8; d <= o.n; d = 2 * (d - 8) + 8)
      if (decision(d) && (d < o.n || escaped)) e++;
    return e * roulette_log2;
  }
  // Whether roulette stops an orbit at Newton step n: a hash of c and n, so deterministic and independent of the
  // orbit's future
  __host__ __device__ bool killed(const State& o) const {
    uint64_t bx, by;
#ifdef __CUDA_ARCH__
    bx = uint64_t(__double_as_longlong(double(o.x)));
    by = uint64_t(__double_as_longlong(double(o.y)));
#else
    const double dx = double(o.x), dy = double(o.y);
    memcpy(&bx, &dx, 8);
    memcpy(&by, &dy, 8);
#endif
    return (mix64(seed ^ mix64(bx ^ mix64(by + uint64_t(o.n)))) & ((uint64_t(1) << roulette_log2) - 1)) != 0;
  }

  __host__ __device__ bool start(State& o, const int64_t i) const {
    if (resume) {
      o = resume[i];
      o.status = 0;
      return false;
    }
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
  __host__ __device__ bool run(State& o) const {
    if (!suspend_at) return o.run(max_iter, burst);
    if (o.n >= suspend_at) { o.status = 10; return true; }  // Unsettled at suspend_at: to the deep queue
    return o.run(max_iter, burst < suspend_at - o.n ? burst : suspend_at - o.n);
  }
  __host__ __device__ int64_t iters(const State& o) const { return o.iters() - iters_offset; }
  __host__ __device__ int64_t progress(const State& o) const { return o.n; }
  // An escaped fast block (status 8, restored to its start n) needs the exact escape step only if some threshold
  // depends on it: escaping at a step in (n, n + 8] puts g below 2^-k for all of them if n - k ≥ 6 and for none
  // if n - k ≤ -3 (escaped_below).  Otherwise finish classifies it as is, skipping the step-by-step redo (its
  // iteration count is then n, up to 7 short).
  __host__ __device__ bool exact_needed(const State& o) const {
    for (int k = 0; k < K; k++)
      if (o.n >= ks[k] - 2 && o.n <= ks[k] + 5) return true;
    return false;
  }
  __host__ __device__ bool pending(const State& o) const {
    return o.pending() && (o.status != 8 || exact_needed(o));
  }
  __host__ __device__ bool immediate(const State& o) const { return o.immediate() && exact_needed(o); }
  __host__ __device__ bool settle(State& o) const {
    const bool newton_step = o.status == 4;
    if (o.settle(max_iter, max_period, newton)) return true;
    if (newton_step && decision(o.n) && killed(o)) {
      o.status = 9;  // Killed by roulette at step n: known not to have escaped by n
      return true;
    }
    return false;
  }
  // Settle work, for sorting settles (engine.h): Newton's candidate period
  __host__ __device__ int settle_key(const State& o) const {
    return o.status != 4 ? 0 : o.candidate && o.n >= max_period ? int(o.candidate) : int(o.atom_candidate(max_period));
  }
  __host__ __device__ void finish(const State& o, const int64_t i) const {
    if (o.status == 10) {
      const uint64_t k = atomic_fetch_add(deep_count, uint64_t(1));
      if (int64_t(k) < deep_cap) { deep_states[k] = o; deep_items[k] = i; }
      bits[i] = 0;  // Filled in by the deep pass
      if (rexp) rexp[i] = 0;
      return;
    }
    uint32_t b = 0;
    // (Status 8: escaped within (n, n + 8]; status 9: killed at n, so escaped after n if ever.  Either way below
    // 2^-k if n - k ≥ 6; roulette requires thresholds clear of its decision steps, so killed samples are unknown
    // only past their step, where their surviving siblings' weights stand in.)
    for (int k = 0; k < K; k++)
      b |= uint32_t(o.status == 8 || o.status == 9 ? o.n - ks[k] >= 6
                                                   : o.status != 1 || escaped_below(o.n, double(o.cx), ks[k])) << k;
    bits[i] = b;
    if (rexp) rexp[i] = uint8_t(roulette_exponent(o) | (o.status == 9 ? kKilled : 0));
    if (iters_out) {
      const int64_t n = o.iters();
      iters_out[i] = n < int64_t(0xffffffff) ? uint32_t(n) : 0xffffffffu;
    }
    if (outcome_out) outcome_out[i] = outcome_code(o.status == 1 || o.status == 8 ? 0 : o.status == 2 ? 1 : 2, o.n);
  }
};

// Group sums of per-sample values over chunks of leaves, for several series at once.  Values: area (bit k of a),
// diff (bit k minus bit k + 1 of a), delta (bit k of b minus bit k of a), or flips (a != b, into s only); with
// swap, a and b trade places (so area of b is a series).  All sums are integers, so exact in any order.
// With roulette, only escapes are reweighted: a sample's value at threshold k past the reference threshold r (the
// last before any roulette decision) is bit_r - W (bit_r - bit_k), with W = 2^e for its e survived decisions, or 0
// if roulette killed it.  Unbiased: an escape at step s between r and k survives roulette with probability 2^-e(s)
// and then counts 2^e(s) times.  A sample that never escapes (interior, or at max_iter) counts exactly bit_r at
// every threshold, so it adds no variance, and differences are weighted sums of escapes alone.
enum Kind { kArea, kDiff, kDelta, kFlips };
struct Series { int8_t kind, k; bool swap; int8_t ref = -1; };  // ref: the roulette reference r, or -1 if none
const int kMaxSeries = 4 * 32 + 1;
// Chunk t (a GPU thread) sums 32 leaves interleaved with its tile's 31 other chunks, so that a warp reads 32
// consecutive leaves at a time
const int64_t kLeafChunk = 32, kLeafTile = 32 * kLeafChunk;

struct alignas(16) Bits4 { uint32_t x, y, z, w; };  // Four samples' bits, one 16-byte load on GPUs

struct ReduceChunk {
  const uint32_t* a;
  const uint32_t* b;
  const uint8_t* rexp;  // Roulette weight exponents, or null
  int m, ss, S;
  int64_t leaves;
  Series series[kMaxSeries];
  int64_t* out;  // [(chunk * S + series) * 3 + {s, q, p}], with ceil(leaves / kLeafTile) * 32 chunks

  template<int kind> __host__ __device__ static int value(const uint32_t x, const uint32_t y, const int k) {
    if constexpr (kind == kArea) return (x >> k) & 1;
    else if constexpr (kind == kDiff) return int((x >> k) & 1) - int((x >> (k + 1)) & 1);
    else if constexpr (kind == kDelta) return int((y >> k) & 1) - int((x >> k) & 1);
    else return x != y;
  }
  // Area at threshold k with roulette (see Series), for bits x and weight byte r
  __host__ __device__ static int64_t roulette_area(const uint32_t x, const uint8_t r, const int k, const int ref) {
    const int bk = (x >> k) & 1;
    if (k <= ref) return bk;
    const int br = (x >> ref) & 1;
    const int64_t W = r & kKilled ? 0 : int64_t(1) << (r & ~kKilled);
    return br - W * (br - bk);
  }

  // Sums of one series over this chunk's leaves, specialized by kind, with a vectorized path for 16 samples per
  // leaf in strata of 4
  template<int kind> __host__ __device__ void sums(const Series e, const int64_t l0, int64_t* o) const {
    const uint32_t* x = e.swap ? b : a;
    const uint32_t* y = e.swap ? a : b;
    const int k = e.k;
    int64_t s = 0, q = 0, p = 0;
    const int64_t hi = l0 + kLeafTile < leaves ? l0 + kLeafTile : leaves;
    if (m == 16 && ss == 4 && !rexp) {
      for (int64_t l = l0; l < hi; l += 32) {
        const Bits4* xv = reinterpret_cast<const Bits4*>(x + 16 * l);
        const Bits4* yv = kind >= kDelta ? reinterpret_cast<const Bits4*>(y + 16 * l) : nullptr;
        int total = 0;
        for (int g = 0; g < 4; g++) {
          const Bits4 u = xv[g], v = kind >= kDelta ? yv[g] : Bits4{0, 0, 0, 0};
          const int cg = value<kind>(u.x, v.x, k) + value<kind>(u.y, v.y, k) + value<kind>(u.z, v.z, k) +
                         value<kind>(u.w, v.w, k);
          total += cg;
          q += cg * cg;
        }
        s += total;
        p += total * total;
      }
    } else {
      for (int64_t l = l0; l < hi; l += 32) {
        int64_t total = 0;
        for (int g = 0; g < m / ss; g++) {
          int64_t cg = 0;
          for (int t = 0; t < ss; t++) {
            const int64_t i = l * m + g * ss + t;
            if (kind <= kDiff && rexp) {
              const int64_t v = roulette_area(x[i], rexp[i], k, e.ref);
              cg += kind == kArea ? v : v - roulette_area(x[i], rexp[i], k + 1, e.ref);
            } else {
              cg += value<kind>(x[i], kind >= kDelta ? y[i] : 0, k);
            }
          }
          total += cg;
          q += cg * cg;
        }
        s += total;
        p += total * total;
      }
    }
    o[0] = s; o[1] = q; o[2] = p;
  }

  __host__ __device__ void operator()(const int64_t c) const {
    const int64_t l0 = c / 32 * kLeafTile + c % 32;
    for (int e = 0; e < S; e++) {
      int64_t* o = out + (c * S + e) * 3;
      switch (series[e].kind) {
        case kArea: sums<kArea>(series[e], l0, o); break;
        case kDiff: sums<kDiff>(series[e], l0, o); break;
        case kDelta: sums<kDelta>(series[e], l0, o); break;
        default: sums<kFlips>(series[e], l0, o); break;
      }
    }
  }
};

// Flip statistics (compare mode, --flip-stats), over chunks of samples whose classifications differ between the
// precisions (a: double, b: the alternative).  Layout (flip_stats_size): [9] flips by outcome pair (kind_a · 3 +
// kind_b, kinds as outcome_code), [9 · K] their net contributions to A_b(k) - A_a(k) (in samples), [33] flips with
// both escaped by the ratio of escape steps (bins of a quarter octave of n_b / n_a, centered, clamped at ±4
// octaves), [34 · 2] flips with both escaped by the octave of n_a and whether b escaped later.
const int64_t kFlipChunk = 1 << 16;
struct FlipChunk {
  const uint32_t *a, *b, *oa, *ob;
  int64_t n;
  int K;
  int64_t* out;
  __host__ __device__ void operator()(const int64_t c) const {
    const int S = flip_stats_size(K);
    int64_t* o = out + c * S;
    for (int j = 0; j < S; j++) o[j] = 0;
    const int64_t hi = (c + 1) * kFlipChunk < n ? (c + 1) * kFlipChunk : n;
    for (int64_t i = c * kFlipChunk; i < hi; i++) {
      if (a[i] == b[i]) continue;
      const int ka = int(oa[i] >> 30), kb = int(ob[i] >> 30), cat = ka * 3 + kb;
      o[cat]++;
      for (int k = 0; k < K; k++) o[9 + cat * K + k] += int((b[i] >> k) & 1) - int((a[i] >> k) & 1);
      if (!ka && !kb) {
        const double na = 4.0 * (oa[i] & 0x3fffffff) + 2, nb = 4.0 * (ob[i] & 0x3fffffff) + 2;
        const int bin = int(std::floor(4 * std::log2(nb / na) + 16.5));
        o[9 + 9 * K + (bin < 0 ? 0 : bin > 32 ? 32 : bin)]++;
        const int oct = int(std::log2(na));
        o[9 + 9 * K + 33 + 2 * (oct < 0 ? 0 : oct > 33 ? 33 : oct) + (nb > na)]++;
      }
    }
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
      if (!certified(s) || !s) continue;
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
// Sums of groups of rows of a rows × S array: out[g * S + j] = Σ in[r * S + j] over rows r of group g
const int64_t kSumGroup = 256;
struct SumGroups {
  const int64_t* in;
  int64_t rows, S;
  int64_t* out;
  __host__ __device__ void operator()(const int64_t t) const {
    const int64_t g = t / S, j = t % S, hi = (g + 1) * kSumGroup < rows ? (g + 1) * kSumGroup : rows;
    int64_t s = 0;
    for (int64_t r = g * kSumGroup; r < hi; r++) s += in[r * S + j];
    out[t] = s;
  }
};

// Column sums of a chunks × S array, in parallel passes over groups of rows (integers, so exact in any order)
vector<int64_t> sum_chunks(const Mem<int64_t>& in, const int64_t chunks, const int64_t S, const bool cuda) {
  if (chunks > kSumGroup) {
    const int64_t groups = (chunks + kSumGroup - 1) / kSumGroup;
    Mem<int64_t> out(groups * S, cuda);
    for_each(groups * S, SumGroups{in.p, chunks, S, out.p}, cuda);
    return sum_chunks(out, groups, S, cuda);
  }
  vector<int64_t> h(chunks * S), sums(S);
  in.to_host(h.data(), chunks * S);
  for (int64_t c = 0; c < chunks; c++)
    for (int64_t j = 0; j < S; j++) sums[j] += h[c * S + j];
  return sums;
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

// Group sums of each series over leaves, in one pass
vector<GroupSums> reduce(const Mem<uint32_t>& a, const Mem<uint32_t>* b, const vector<Series>& series,
                         const TreeParams& p, const int64_t leaves, const uint8_t* rexp = nullptr) {
  const int S = int(series.size());
  slow_assert(S <= kMaxSeries);
  const int64_t chunks = (leaves + kLeafTile - 1) / kLeafTile * 32;
  static const bool timing = env_int("MANDELBROT_REDUCE_TIMING", 0);
  const auto t0 = std::chrono::steady_clock::now();
  Mem<int64_t> out(chunks * S * 3, p.cuda);
  ReduceChunk r{a.p, b ? b->p : nullptr, rexp, p.m, p.strata * p.strata, S, leaves, {}, out.p};
  for (int e = 0; e < S; e++) r.series[e] = series[e];
  for_each(chunks, r, p.cuda);
  if (timing) { int64_t x; out.to_host(&x, 1); }  // Wait for the pass
  const auto t1 = std::chrono::steady_clock::now();
  const auto h = sum_chunks(out, chunks, S * 3, p.cuda);
  if (timing)
    print("      reduce: %d leaves, %d series: chunks %.1f ms, sums %.1f ms", leaves, S,
          1e3 * std::chrono::duration<double>(t1 - t0).count(),
          1e3 * std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count());
  vector<GroupSums> g(S);
  for (int e = 0; e < S; e++) g[e] = GroupSums{h[3 * e], h[3 * e + 1], h[3 * e + 2]};
  return g;
}

template<class T> SampleTask<T> sample_task(const Cell* leaves, const TreeParams& p, const double w, const double h,
                                           uint32_t* bits, uint32_t* iters, uint32_t* outcomes, uint8_t* rexp) {
  SampleTask<T> task{p.burst, p.sample_min_blocks, leaves, p.m, p.strata, p.seed, p.x0, p.y0, w, h, p.max_iter, p.first_newton, p.newton_max_period,
                     NewtonOptions{p.newton_iters, p.newton_close2, p.newton_tol < 0 ? -1 : p.newton_tol * p.newton_tol,
                                   p.newton_margin, false, p.newton_repel2, false},
                     int(p.ks.size()), {}, bits, iters, outcomes, p.first_newton, p.roulette_from,
                     p.roulette_from ? p.roulette_log2 : 0, p.roulette_stride, rexp,
                     0, nullptr, nullptr, nullptr, 0, nullptr, 0};
  for (size_t k = 0; k < p.ks.size(); k++) task.ks[k] = p.ks[k];
  return task;
}

template<class T> int64_t sample(const Cell* leaves, const int64_t n_leaves, const TreeParams& p,
                                 const double w, const double h, Mem<uint32_t>& bits, int64_t& overflow,
                                 uint32_t* iters = nullptr, uint32_t* outcomes = nullptr, uint8_t* rexp = nullptr) {
  const auto task = sample_task<T>(leaves, p, w, h, bits.p, iters, outcomes, rexp);
  const auto stats = run_orbits(task, n_leaves * p.m, p.cuda);
  overflow += stats.overflow;
  return stats.iters;
}

double secs_since(const std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

// The roulette reference threshold (see Series): the last before the first roulette decision step, or -1 if none
int roulette_reference(const TreeParams& p) {
  if (!p.roulette_from) return -1;
  int64_t first = p.first_newton + 8;
  while (first < p.roulette_from) first = 2 * (first - 8) + 8;
  int ref = -1;
  for (int k = 0; k < int(p.ks.size()); k++)
    if (p.ks[k] < first) ref = k;
  return ref;
}

// --dump-cells: append each certified cell of a level as int32s (depth, ix, iy, below-threshold mask), for pictures
void dump_certified(const TreeParams& p, const Level& level, const int d, const int64_t n, const Mem<uint32_t>& status) {
  vector<uint32_t> st(n);
  status.to_host(st.data(), n);
  vector<int32_t> out;
  for (int64_t i = 0; i < n; i++) {
    if (!certified(st[i])) continue;
    const Cell c = level.at(i);
    out.insert(out.end(), {int32_t(d), c.ix, c.iy, int32_t(st[i])});
  }
  static std::mutex mu;
  std::lock_guard<std::mutex> lock(mu);
  FILE* f = fopen(p.dump_cells.c_str(), "ab");
  slow_assert(f, "can't open %s", p.dump_cells);
  slow_assert(fwrite(out.data(), 4, out.size(), f) == out.size(), "short write to %s", p.dump_cells);
  fclose(f);
}

// --dump: append each leaf's cell (ix, iy at the leaf level) and its samples' outcome codes, as int32s
void dump_leaves(const TreeParams& p, const Mem<Cell>& leaves, const int64_t l0, const int64_t nl,
                 const Mem<uint32_t>& outcomes) {
  vector<Cell> all(l0 + nl);
  vector<uint32_t> codes(nl * p.m);
  leaves.to_host(all.data(), l0 + nl);
  const Cell* cells = all.data() + l0;
  outcomes.to_host(codes.data(), codes.size());
  vector<uint32_t> out;
  out.reserve(nl * (2 + p.m));
  for (int64_t i = 0; i < nl; i++) {
    out.push_back(uint32_t(cells[i].ix));
    out.push_back(uint32_t(cells[i].iy));
    for (int j = 0; j < p.m; j++) out.push_back(codes[i * p.m + j]);
  }
  static std::mutex mu;
  std::lock_guard<std::mutex> lock(mu);
  FILE* f = fopen(p.dump.c_str(), "ab");
  slow_assert(f, "can't open %s", p.dump);
  slow_assert(fwrite(out.data(), 4, out.size(), f) == out.size(), "short write to %s", p.dump);
  fclose(f);
}

// The leaf series of a run (areas and consecutive differences, and with compare the alternative's areas, its
// difference from double, and flips), and where each one's sums go in R
void leaf_series(const TreeParams& p, TreeResult& R, vector<Series>& series, vector<GroupSums*>& into,
                 GroupSums* flips) {
  const bool compare = p.prec.starts_with("compare");
  const int K = p.ks.size();
  for (int k = 0; k < K; k++) {
    const int8_t k8 = int8_t(k), ref = int8_t(roulette_reference(p));
    series.push_back({kArea, k8, false, ref}); into.push_back(&R.area[k]);
    if (k + 1 < K) { series.push_back({kDiff, k8, false, ref}); into.push_back(&R.diff[k]); }
    if (compare) {
      series.push_back({kArea, k8, true}); into.push_back(&R.float_area[k]);
      series.push_back({kDelta, k8, false}); into.push_back(&R.delta[k]);
    }
  }
  if (compare) { series.push_back({kFlips, 0, false}); into.push_back(flips); }
}

// Deep queue (TreeParams::deep_from).  A sub-batch's leaves holding a suspended sample are saved here (their m
// result words and roulette weights) and zeroed in the sub-batch, which makes their contribution to every series
// exactly zero there; the suspended states wait here too, until a deep pass runs them all at full GPU width, fills
// in their words, and reduces the saved leaves, adding exactly what the sub-batch would have.
struct DeepQueue {
  std::mutex mu;
  vector<char> states;    // Suspended orbit states, packed
  vector<int64_t> slots;  // Each state's place in bits and rexp below
  vector<uint32_t> bits;  // Saved leaves' result words, m per leaf
  vector<uint8_t> rexp;   // Their roulette weights (with roulette)
};

// Copy leaves[u]'s m words of bits (and rexp) to out, then zero them in bits
struct GatherLeaves {
  uint32_t* bits;
  const uint8_t* rexp;
  const int64_t* leaves;
  int m;
  uint32_t* out_bits;
  uint8_t* out_rexp;
  __host__ __device__ void operator()(const int64_t u) const {
    const int64_t l = leaves[u];
    for (int j = 0; j < m; j++) {
      out_bits[u * m + j] = bits[l * m + j];
      if (rexp) out_rexp[u * m + j] = rexp[l * m + j];
      bits[l * m + j] = 0;
    }
  }
};

// Sample a sub-batch with suspension at p.deep_from, moving suspended samples and their leaves to the queue
template<class T> int64_t sample_suspending(const Cell* leaves, const int64_t nl, const TreeParams& p,
                                            const double w, const double h, Mem<uint32_t>& bits, Mem<uint8_t>& rexp,
                                            int64_t& overflow, DeepQueue& q) {
  typedef Orbit<T> State;
  const int64_t n = nl * p.m;
  // Room for the suspended states: usually a few per million samples, but dense regions (cusps, the real axis)
  // can have many more, so if they overflow, run the sub-batch again with exactly enough room (every sample
  // rewrites its result, so the rerun is exact)
  int64_t cap = p.deep_cap ? p.deep_cap : std::max<int64_t>(int64_t(1) << 16, n / 64);
  Mem<State> states(0, p.cuda);
  Mem<int64_t> items(0, p.cuda);
  Mem<uint64_t> count(1, p.cuda);
  RunStats stats;
  uint64_t c;
  for (;;) {
    Mem<State> s(cap, p.cuda);
    Mem<int64_t> it(cap, p.cuda);
    std::swap(states.p, s.p); std::swap(states.n, s.n);
    std::swap(items.p, it.p); std::swap(items.n, it.n);
    count.zero();
    auto task = sample_task<T>(leaves, p, w, h, bits.p, nullptr, nullptr, p.roulette_from ? rexp.p : nullptr);
    task.suspend_at = p.deep_from;
    task.deep_states = states.p;
    task.deep_items = items.p;
    task.deep_count = count.p;
    task.deep_cap = cap;
    stats = run_orbits(task, n, p.cuda);
    count.to_host(&c, 1);
    if (int64_t(c) <= cap) break;
    print("  deep queue: %d suspended samples in a sub-batch of %d, room for %d: rerunning", c, n, cap);
    cap = c;
  }
  overflow += stats.overflow;
  if (!c) return stats.iters;
  vector<int64_t> hi(c);
  vector<State> hs(c);
  items.to_host(hi.data(), c);
  states.to_host(hs.data(), c);
  // Their leaves, each once
  vector<int64_t> ls(c);
  for (uint64_t k = 0; k < c; k++) ls[k] = hi[k] / p.m;
  std::sort(ls.begin(), ls.end());
  ls.erase(std::unique(ls.begin(), ls.end()), ls.end());
  const int64_t L = ls.size();
  Mem<int64_t> dl(L, p.cuda);
  dl.from_host(ls.data(), L);
  Mem<uint32_t> gb(L * p.m, p.cuda);
  Mem<uint8_t> gr(p.roulette_from ? L * p.m : 0, p.cuda);
  for_each(L, GatherLeaves{bits.p, p.roulette_from ? rexp.p : nullptr, dl.p, p.m, gb.p, gr.p}, p.cuda);
  vector<uint32_t> hb(L * p.m);
  vector<uint8_t> hr(p.roulette_from ? L * p.m : 0);
  gb.to_host(hb.data(), hb.size());
  gr.to_host(hr.data(), hr.size());
  std::lock_guard<std::mutex> lock(q.mu);
  const int64_t base = q.bits.size();
  q.bits.insert(q.bits.end(), hb.begin(), hb.end());
  q.rexp.insert(q.rexp.end(), hr.begin(), hr.end());
  const char* raw = reinterpret_cast<const char*>(hs.data());
  q.states.insert(q.states.end(), raw, raw + c * sizeof(State));
  for (uint64_t k = 0; k < c; k++) {
    const int64_t u = std::lower_bound(ls.begin(), ls.end(), hi[k] / p.m) - ls.begin();
    q.slots.push_back(base + u * p.m + hi[k] % p.m);
  }
  return stats.iters;
}

// A deep pass: run every queued state to completion, fill in its saved leaf, and reduce the saved leaves into R
template<class T> void deep_pass(const TreeParams& p, DeepQueue& q, TreeResult& R) {
  typedef Orbit<T> State;
  vector<char> raw;
  vector<int64_t> slots;
  vector<uint32_t> bits;
  vector<uint8_t> rexp;
  {
    std::lock_guard<std::mutex> lock(q.mu);
    std::swap(raw, q.states); std::swap(slots, q.slots); std::swap(bits, q.bits); std::swap(rexp, q.rexp);
  }
  const int64_t c = slots.size();
  if (!c) return;
  const auto t0 = std::chrono::steady_clock::now();
  Mem<State> st(c, p.cuda);
  st.from_host(reinterpret_cast<const State*>(raw.data()), c);
  Mem<uint32_t> db(c, p.cuda);
  Mem<uint8_t> dr(p.roulette_from ? c : 0, p.cuda);
  auto task = sample_task<T>(nullptr, p, 0, 0, db.p, nullptr, nullptr, p.roulette_from ? dr.p : nullptr);
  task.resume = st.p;
  task.iters_offset = p.deep_from;
  const auto stats = run_orbits(task, c, p.cuda);
  vector<uint32_t> hb(c);
  vector<uint8_t> hr(dr.n);
  db.to_host(hb.data(), c);
  dr.to_host(hr.data(), hr.size());
  for (int64_t k = 0; k < c; k++) {
    bits[slots[k]] = hb[k];
    if (p.roulette_from) rexp[slots[k]] = hr[k];
  }
  // Reduce the saved leaves on the host, with the same integer sums as a sub-batch
  TreeParams ph = p;
  ph.cuda = false;
  const int64_t L = bits.size() / p.m;
  Mem<uint32_t> mb(bits.size(), false);
  mb.from_host(bits.data(), bits.size());
  vector<Series> series;
  vector<GroupSums*> into;
  leaf_series(p, R, series, into, nullptr);
  const auto sums = reduce(mb, nullptr, series, ph, L, p.roulette_from ? rexp.data() : nullptr);
  for (size_t e = 0; e < sums.size(); e++) *into[e] += sums[e];
  R.leaf_iters += stats.iters;
  R.overflow += stats.overflow;
  R.deep_samples += c;
  R.deep_passes++;
  R.sample_secs += secs_since(t0);
}

void deep_pass(const TreeParams& p, DeepQueue& q, TreeResult& R) {
  if (p.prec == "dd") deep_pass<Expansion<2>>(p, q, R);
  else deep_pass<double>(p, q, R);
}

// Cells (and center hints) at TreeParams::split_depth, collected by the first phase of run_tree
struct SplitCells {
  std::mutex mu;
  vector<Cell> cells;
  vector<uint8_t> hints;
};

// One batch: base cells [cell0, cell1), with tree levels, leaf samples and reductions accumulated into Rb.
// Returns false, leaving Rb partial, if the base level itself has p.max_level_cells cells: the caller splits the
// batch.  A deeper level that large is split into pieces of independent subtrees, each continued from that level
// (d0 > 0, with its cells and hints in start_cells and start_hints, which this takes).  With split, the batch
// stops at p.split_depth, appending the cells there to split instead of going on.
bool run_batch(const TreeParams& p, const int64_t cell0, const int64_t cell1, const int64_t max_leaves,
               TreeResult& Rb, DeepQueue& deep, const int d0 = 0, Mem<Cell>* start_cells = nullptr,
               Mem<uint8_t>* start_hints = nullptr, SplitCells* split = nullptr) {
  const int K = p.ks.size();
  const bool compare = p.prec.starts_with("compare"), single = p.prec == "float";

  // Tree levels
  const auto t1 = std::chrono::steady_clock::now();
  int64_t n = cell1 - cell0;
  Mem<Cell> cells(0, p.cuda), leaves(0, p.cuda);
  Mem<uint8_t> hints(0, p.cuda);
  if (d0) {
    std::swap(cells.p, start_cells->p); std::swap(cells.n, start_cells->n);
    std::swap(hints.p, start_hints->p); std::swap(hints.n, start_hints->n);
    n = cells.n;
  }
  int64_t n_leaves = 0;
  for (int d = d0; d <= p.depth; d++) {
    if (split && d == p.split_depth) {
      vector<Cell> hc(n);
      vector<uint8_t> hh(n);
      cells.to_host(hc.data(), n);
      hints.to_host(hh.data(), n);
      std::lock_guard<std::mutex> lock(split->mu);
      split->cells.insert(split->cells.end(), hc.begin(), hc.end());
      split->hints.insert(split->hints.end(), hh.begin(), hh.end());
      Rb.split_cells += n;
      Rb.tree_secs += secs_since(t1);
      return true;
    }
    if (n >= p.max_level_cells) {
      if (!d) return false;
      // Split this level's cells into pieces of independent subtrees, and continue each from here
      vector<Cell> hc(n);
      vector<uint8_t> hh(n);
      cells.to_host(hc.data(), n);
      hints.to_host(hh.data(), n);
      const int64_t piece = std::max<int64_t>(1, p.max_level_cells / 2);
      Rb.tree_secs += secs_since(t1);
      for (int64_t a = 0; a < n; a += piece) {
        const int64_t m = std::min(piece, n - a);
        Mem<Cell> pc(m, p.cuda);
        Mem<uint8_t> ph(m, p.cuda);
        pc.from_host(hc.data() + a, m);
        ph.from_host(hh.data() + a, m);
        run_batch(p, cell0, cell1, max_leaves, Rb, deep, d, &pc, &ph, split);
      }
      return true;
    }
    const Level level{d ? cells.p : nullptr, p.base, cell0, p.shard, p.shards};
    const double w = (p.x1 - p.x0) / double(p.base << d), h = (p.y1 - p.y0) / double(p.base << d);
    Mem<uint32_t> status(n, p.cuda);
    slow_assert(p.center_max_iter < (int64_t(1) << 31), "CenterTask::Record needs 32-bit steps");
    CenterTask task{p.burst, p.center_min_blocks, level, p.x0, p.y0, w, h, 0.5 * std::hypot(w, h), p.safety, std::min(p.max_iter, p.center_max_iter),
                    p.center_first_newton, p.center_max_period, K, {},
                    status.p, d ? hints.p : nullptr, p.center_hint_newton};
    for (int k = 0; k < K; k++) task.ks[k] = p.ks[k];
    const auto stats = run_orbits(task, n, p.cuda);
    if (!p.dump_cells.empty()) dump_certified(p, level, d, n, status);
    Rb.center_kernel_secs += stats.secs;
    Rb.centers += n;
    Rb.center_iters += stats.iters;
    Rb.overflow += stats.overflow;

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
      Rb.exact[d] += o[1];
      for (int k = 0; k < K; k++) Rb.certified[d * K + k] += o[2 + k];
    }
    if (p.tiles) {
      const int64_t tc = (n + kTileChunk - 1) / kTileChunk, S = int64_t(p.tiles) * p.tiles * K;
      Mem<int64_t> tout(tc * S, p.cuda);
      for_each(tc, TileCountChunk{status.p, level, n, p.base << d, int64_t(1) << (2 * (p.depth - d)), p.tiles, K,
                                  tout.p}, p.cuda);
      const auto h = sum_chunks(tout, tc, S, p.cuda);
      for (int64_t j = 0; j < S; j++) Rb.tile_cert[j] += h[j];
    }
    Mem<int64_t> doffsets(chunks, p.cuda);
    doffsets.from_host(offsets.data(), chunks);
    const bool last = d == p.depth;
    Mem<Cell> next(last ? uncertified : 4 * uncertified, p.cuda);
    Mem<uint8_t> next_hints(last ? 0 : 4 * uncertified, p.cuda);
    for_each(chunks, EmitChunk{status.p, level, n, doffsets.p, !last, next.p, next_hints.p}, p.cuda);
    if (last) {
      std::swap(leaves.p, next.p); std::swap(leaves.n, next.n);
      n_leaves = uncertified;
    } else {
      std::swap(cells.p, next.p); std::swap(cells.n, next.n);
      std::swap(hints.p, next_hints.p); std::swap(hints.n, next_hints.n);
      n = 4 * uncertified;
    }
  }
  Rb.tree_secs += secs_since(t1);

  // Leaf samples and reductions, in sub-batches of fewer than 2^31 samples
  const double w = (p.x1 - p.x0) / double(p.base << p.depth), h = (p.y1 - p.y0) / double(p.base << p.depth);
  for (int64_t l0 = 0; l0 < n_leaves; l0 += max_leaves) {
    const int64_t nl = std::min(max_leaves, n_leaves - l0);
    const auto t2 = std::chrono::steady_clock::now();
    Mem<uint32_t> bits(nl * p.m, p.cuda), fbits(compare ? nl * p.m : 0, p.cuda),
                  iters(p.leaf_stats ? nl * p.m : 0, p.cuda),
                  outcomes(p.flip_stats || !p.dump.empty() ? nl * p.m : 0, p.cuda),
                  foutcomes(p.flip_stats ? nl * p.m : 0, p.cuda);
    uint32_t* ip = p.leaf_stats ? iters.p : nullptr;
    uint32_t* op = p.flip_stats || !p.dump.empty() ? outcomes.p : nullptr;
    uint32_t* fop = p.flip_stats ? foutcomes.p : nullptr;
    Mem<uint8_t> rexp(p.roulette_from ? nl * p.m : 0, p.cuda);
    uint8_t* rp = p.roulette_from ? rexp.p : nullptr;
    const Cell* lp = leaves.p + l0;
    if (p.deep_from)
      Rb.leaf_iters += p.prec == "dd" ? sample_suspending<Expansion<2>>(lp, nl, p, w, h, bits, rexp, Rb.overflow, deep)
                                      : sample_suspending<double>(lp, nl, p, w, h, bits, rexp, Rb.overflow, deep);
    else
      Rb.leaf_iters += single ? sample<float>(lp, nl, p, w, h, bits, Rb.overflow, ip)
                     : p.prec == "dd" ? sample<Expansion<2>>(lp, nl, p, w, h, bits, Rb.overflow, ip, op, rp)
                                      : sample<double>(lp, nl, p, w, h, bits, Rb.overflow, ip, op, rp);
    if (compare) {
      // The alternative precision: float, or double rounded to fewer bits
      Rb.leaf_iters += p.prec == "compare24" ? sample<Rounded<24>>(lp, nl, p, w, h, fbits, Rb.overflow, nullptr, fop)
                    : p.prec == "compare27" ? sample<Rounded<27>>(lp, nl, p, w, h, fbits, Rb.overflow, nullptr, fop)
                    : p.prec == "compare30" ? sample<Rounded<30>>(lp, nl, p, w, h, fbits, Rb.overflow, nullptr, fop)
                    : p.prec == "compare36" ? sample<Rounded<36>>(lp, nl, p, w, h, fbits, Rb.overflow, nullptr, fop)
                    : p.prec == "compare42" ? sample<Rounded<42>>(lp, nl, p, w, h, fbits, Rb.overflow, nullptr, fop)
                    : p.prec == "compare48" ? sample<Rounded<48>>(lp, nl, p, w, h, fbits, Rb.overflow, nullptr, fop)
                    : p.prec == "comparedd" ? sample<Expansion<2>>(lp, nl, p, w, h, fbits, Rb.overflow, nullptr, fop)
                                            : sample<float>(lp, nl, p, w, h, fbits, Rb.overflow, nullptr, fop);
    }
    if (!p.dump.empty()) dump_leaves(p, leaves, l0, nl, outcomes);
    if (p.flip_stats) {
      const int64_t samples = nl * p.m, chunks = (samples + kFlipChunk - 1) / kFlipChunk, S = flip_stats_size(K);
      Mem<int64_t> fout(chunks * S, p.cuda);
      for_each(chunks, FlipChunk{bits.p, fbits.p, outcomes.p, foutcomes.p, samples, K, fout.p}, p.cuda);
      const auto h = sum_chunks(fout, chunks, S, p.cuda);
      for (int64_t j = 0; j < S; j++) Rb.flip_stats[j] += h[j];
    }
    Rb.sample_secs += secs_since(t2);

    const auto t3 = std::chrono::steady_clock::now();
    // Every series in one pass: areas and consecutive differences, and with compare, the alternative's areas, its
    // difference from double, and flips
    vector<Series> series;
    vector<GroupSums*> into;
    GroupSums flips;
    leaf_series(p, Rb, series, into, &flips);
    const auto sums = reduce(bits, compare ? &fbits : nullptr, series, p, nl, p.roulette_from ? rexp.p : nullptr);
    for (size_t e = 0; e < sums.size(); e++) *into[e] += sums[e];
    Rb.flips += flips.s;
    if (p.tiles) {
      const int64_t tc = (nl + kTileChunk - 1) / kTileChunk, S = int64_t(p.tiles) * p.tiles * K * 3;
      Mem<int64_t> tout(tc * S, p.cuda);
      for_each(tc, TileReduceChunk{bits.p, leaves.p + l0, nl, p.base << p.depth, p.tiles, K, p.m, p.strata * p.strata,
                                   tout.p}, p.cuda);
      const auto h = sum_chunks(tout, tc, S, p.cuda);
      for (int64_t j = 0; j < S / 3; j++) Rb.tile_diff[j] += GroupSums{h[3 * j], h[3 * j + 1], h[3 * j + 2]};
    }
    if (p.leaf_stats) {
      const int64_t chunks = (nl + kStatsChunk - 1) / kStatsChunk, size = int64_t(K) * 108;
      Mem<int64_t> out(chunks * size, p.cuda);
      for_each(chunks, LeafStatsChunk{bits.p, iters.p, K, nl, out.p}, p.cuda);
      vector<int64_t> h(chunks * size);
      out.to_host(h.data(), h.size());
      for (int64_t c = 0; c < chunks; c++)
        for (int64_t j = 0; j < size; j++) Rb.alloc[j] += h[c * size + j];
    }
    Rb.reduce_secs += secs_since(t3);
  }
  if (p.deep_from) {
    int64_t queued;
    {
      std::lock_guard<std::mutex> lock(deep.mu);
      queued = deep.slots.size();
    }
    if (queued >= p.deep_batch) deep_pass(p, deep, Rb);
  }
  Rb.leaves += n_leaves;
  Rb.batches++;
  return true;
}

// Empty per-threshold sums shaped like R's, for a batch
TreeResult empty_like(const TreeResult& R) {
  TreeResult Rb;
  Rb.p = R.p;
  Rb.certified.assign(R.certified.size(), 0);
  Rb.exact.assign(R.exact.size(), 0);
  Rb.area.resize(R.area.size()); Rb.diff.resize(R.diff.size());
  Rb.float_area.resize(R.float_area.size()); Rb.delta.resize(R.delta.size());
  Rb.alloc.assign(R.alloc.size(), 0);
  Rb.tile_diff.resize(R.tile_diff.size()); Rb.tile_cert.assign(R.tile_cert.size(), 0);
  Rb.flip_stats.assign(R.flip_stats.size(), 0);
  return Rb;
}

}  // namespace

// R += Rb (integer sums, so the order of batches does not matter)
void merge(TreeResult& R, const TreeResult& Rb) {
  for (size_t i = 0; i < R.certified.size(); i++) R.certified[i] += Rb.certified[i];
  for (size_t i = 0; i < R.exact.size(); i++) R.exact[i] += Rb.exact[i];
  for (size_t i = 0; i < R.area.size(); i++) {
    R.area[i] += Rb.area[i]; R.diff[i] += Rb.diff[i]; R.float_area[i] += Rb.float_area[i]; R.delta[i] += Rb.delta[i];
  }
  for (size_t i = 0; i < R.alloc.size(); i++) R.alloc[i] += Rb.alloc[i];
  for (size_t i = 0; i < R.tile_diff.size(); i++) R.tile_diff[i] += Rb.tile_diff[i];
  for (size_t i = 0; i < R.tile_cert.size(); i++) R.tile_cert[i] += Rb.tile_cert[i];
  for (size_t i = 0; i < R.flip_stats.size(); i++) R.flip_stats[i] += Rb.flip_stats[i];
  R.leaves += Rb.leaves; R.centers += Rb.centers; R.center_iters += Rb.center_iters; R.leaf_iters += Rb.leaf_iters;
  R.overflow += Rb.overflow; R.flips += Rb.flips; R.batches += Rb.batches;
  R.deep_samples += Rb.deep_samples; R.deep_passes += Rb.deep_passes; R.split_cells += Rb.split_cells;
  R.tree_secs += Rb.tree_secs; R.center_kernel_secs += Rb.center_kernel_secs; R.sample_secs += Rb.sample_secs;
  R.reduce_secs += Rb.reduce_secs;
}

namespace {

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

// A result with no cells yet, sized for p
TreeResult empty_result(const TreeParams& p) {
  const int K = p.ks.size();
  TreeResult R;
  R.p = p;
  R.certified.assign((p.depth + 1) * K, 0);
  R.exact.assign(p.depth + 1, 0);
  R.area.resize(K); R.diff.resize(K); R.float_area.resize(K); R.delta.resize(K);
  slow_assert(!p.leaf_stats || (p.m == 16 && !p.prec.starts_with("compare")), "leaf_stats needs m = 16, no compare");
  if (p.leaf_stats) R.alloc.assign(int64_t(K) * 108, 0);
  slow_assert(p.tiles >= 0 && p.tiles <= 64, "tiles must be in [0, 64]");
  if (p.tiles) { R.tile_diff.resize(int64_t(p.tiles) * p.tiles * K); R.tile_cert.assign(int64_t(p.tiles) * p.tiles * K, 0); }
  slow_assert(!p.flip_stats || p.prec.starts_with("compare"), "flip_stats needs a compare precision");
  slow_assert(p.dump_cells.empty() || !p.cuda, "dump_cells is CPU only");
  slow_assert(!p.deep_from || ((p.prec == "double" || p.prec == "dd") && !p.tiles && !p.leaf_stats &&
                               p.dump.empty() && p.deep_from % 8 == 0 && p.deep_batch > 0),
              "deep queue needs prec double or dd, no tiles, leaf stats or dump, and deep_from a multiple of 8");
  if (p.roulette_from) {
    slow_assert((p.prec == "double" || p.prec == "dd") && !p.tiles && !p.leaf_stats,
                "roulette needs prec double or dd, no tiles or leaf stats");
    slow_assert(p.first_newton % 8 == 0 && 1 <= p.roulette_log2 && p.roulette_log2 <= 4, "bad roulette settings");
    slow_assert(p.roulette_stride >= 1, "roulette stride must be positive");
    int decisions = 0, j = 0;
    for (int64_t d = p.first_newton + 8; d <= p.max_iter; d = 2 * (d - 8) + 8)
      if (d >= p.roulette_from) decisions += j++ % p.roulette_stride == 0;
    slow_assert(decisions * p.roulette_log2 <= 24, "roulette weights up to 2^%d overflow the sums (at most 2^24)",
                decisions * p.roulette_log2);
    slow_assert(roulette_reference(p) >= 0, "roulette needs a threshold before its first decision step");
    for (const auto k : p.ks)
      for (int64_t d = p.first_newton + 8; d <= p.max_iter; d = 2 * (d - 8) + 8)
        slow_assert(d < p.roulette_from || std::abs(d - k) >= 8,
                    "threshold %d is within 8 steps of roulette decision step %d", k, d);
  }
  if (p.flip_stats) R.flip_stats.assign(flip_stats_size(K), 0);
  return R;
}

TreeResult run_tree(const TreeParams& p) {
  const int K = p.ks.size();
  slow_assert(0 < K && K <= 31, "need 1 to 31 thresholds, got %d", K);
  slow_assert(p.max_iter < kOrbitNever, "max_iter %d needs a -DMANDELBROT_ORBIT64 build", p.max_iter);
  slow_assert(p.prec == "double" || p.prec == "dd" || p.prec == "float" || p.prec == "compare" ||
              p.prec == "compare24" || p.prec == "compare27" || p.prec == "compare30" ||
              p.prec == "compare36" || p.prec == "compare42" || p.prec == "compare48" || p.prec == "comparedd",
              "bad prec %s", p.prec);
  const int ss = p.strata * p.strata;
  slow_assert(p.strata >= 1 && p.m % ss == 0 && p.m / ss >= 2,
              "need m a multiple of strata^2 with at least 2 groups for variance estimates");
  slow_assert((p.base << p.depth) < (int64_t(1) << 31), "grid too fine for 32-bit cell coordinates");

  TreeResult R = empty_result(p);
  const auto t0 = std::chrono::steady_clock::now();
  const int64_t all_rows = p.rows < 0 ? p.base : std::min(p.rows, p.base);
  slow_assert(0 <= p.shard && p.shard < p.shards, "bad shard %d of %d", p.shard, p.shards);
  const int64_t rows = p.shard < all_rows ? (all_rows - p.shard + p.shards - 1) / p.shards : 0;  // This shard's

  // Batches are runs of base cells in row-major order, sized adaptively for about p.batch leaves (and fewer
  // than 2^31 samples).  On the GPU, p.overlap batches run at once in host threads with their own streams, so
  // that one batch's latency-bound work (tree levels, and the long-orbit tails of its rounds) overlaps another's
  // full-GPU work.  With p.split_depth, levels below it are built first for whole runs of base cells (a few
  // large passes instead of every batch building them), collecting the cells at split_depth; batches are then
  // runs of those cells, each continued from split_depth.  Results are the same either way.
  const int64_t total = rows * p.base, max_leaves = std::min<int64_t>(p.batch, ((int64_t(1) << 31) - 1) / p.m);
  std::mutex mu;
  DeepQueue deep;  // Suspended samples (p.deep_from), shared by the workers
  SplitCells split;  // Cells at p.split_depth, from the first phase
  const auto start = std::chrono::steady_clock::now();
  auto last_progress = start;
  const int threads = p.cuda ? p.overlap : 1;
  slow_assert(threads >= 1, "overlap must be at least 1");
  slow_assert(0 <= p.split_depth && p.split_depth <= p.depth, "split depth %d outside [0, depth]", p.split_depth);

  // Run units [0, units) in adaptive batches over the workers: run(a, b, Rb) runs units [a, b) into Rb, returning
  // false if they must be split; size(Rb) is a batch's size, aimed at target per batch
  const auto pool = [&](const char* phase, const int64_t units, const int64_t first,
                        const std::function<bool(int64_t, int64_t, TreeResult&)>& run,
                        const std::function<int64_t(const TreeResult&)>& size, const int64_t target) {
    int64_t next = 0, per_batch = std::min<int64_t>(units, first), done = 0, done_size = 0;
    const auto worker = [&]() {
      for (;;) {
        int64_t u0, u1;
        {
          std::lock_guard<std::mutex> lock(mu);
          if (next >= units) return;
          u0 = next;
          u1 = std::min(units, u0 + per_batch);
          next = u1;
        }
        // Run [a, b), halving it while it must be split
        const std::function<void(int64_t, int64_t)> go = [&](const int64_t a, const int64_t b) {
          TreeResult Rb = empty_like(R);
          if (!run(a, b, Rb)) {
            const int64_t m = a + (b - a) / 2;
            go(a, m);
            go(m, b);
            return;
          }
          std::lock_guard<std::mutex> lock(mu);
          merge(R, Rb);
          // Aim for target per batch, from the size per unit so far, growing at most 4× per batch since sparse
          // early units (common in --box domains) underestimate the density
          done += b - a;
          done_size += size(Rb);
          const double per_unit = double(done_size) / double(done);
          per_batch = std::max<int64_t>(1, std::min<int64_t>(4 * (b - a), int64_t(double(target) / std::max(1e-3, per_unit))));
          if (p.progress > 0 && secs_since(last_progress) >= p.progress) {
            last_progress = std::chrono::steady_clock::now();
            string mem;
            IF_CUDA(if (p.cuda) mem = "; " + gpu_memory();)
            print("progress %.0f s: %s %d / %d, batches %d, leaves %.3g, leaf iterations %.3g, center "
                  "iterations %.3g (sampling %.0f s, tree %.0f s)%s", secs_since(start), phase, done, units, R.batches,
                  double(R.leaves), double(R.leaf_iters), double(R.center_iters), R.sample_secs, R.tree_secs, mem);
          }
        };
        go(u0, u1);
      }
    };
    if (threads == 1) worker();
    else {
      vector<std::thread> workers;
      for (int t = 0; t < threads; t++) workers.emplace_back(worker);
      for (auto& t : workers) t.join();
    }
  };
  // A small first batch, to estimate density: 64 base cells, fewer for deep trees, whose dense base cells have
  // many more leaves
  const int64_t first = std::max<int64_t>(1, 64 >> std::min(6, 2 * std::max(0, p.depth - 8)));
  const auto leaves_of = [](const TreeResult& Rb) { return Rb.leaves; };
  if (!p.split_depth) {
    pool("base cells", total, first,
         [&](const int64_t a, const int64_t b, TreeResult& Rb) { return run_batch(p, a, b, max_leaves, Rb, deep); },
         leaves_of, max_leaves);
  } else {
    // Levels below split_depth for runs of base cells, then batches of the cells there (for the first phase,
    // aiming for a quarter of max_level_cells collected cells per batch)
    pool("base cells", total, first,
         [&](const int64_t a, const int64_t b, TreeResult& Rb) {
           return run_batch(p, a, b, max_leaves, Rb, deep, 0, nullptr, nullptr, &split);
         },
         [](const TreeResult& Rb) { return Rb.split_cells; }, std::max<int64_t>(1, p.max_level_cells / 4));
    const int64_t n = split.cells.size();
    pool("split cells", n, 1,
         [&](const int64_t a, const int64_t b, TreeResult& Rb) {
           Mem<Cell> c(b - a, p.cuda);
           Mem<uint8_t> h(b - a, p.cuda);
           c.from_host(split.cells.data() + a, b - a);
           h.from_host(split.hints.data() + a, b - a);
           return run_batch(p, 0, 0, max_leaves, Rb, deep, p.split_depth, &c, &h);
         },
         leaves_of, max_leaves);
  }
  if (p.deep_from) {
    // The last deep pass, for whatever is still queued
    TreeResult Rb = empty_like(R);
    deep_pass(p, deep, Rb);
    merge(R, Rb);
  }
  R.secs = secs_since(t0);
  return R;
}

namespace {

// The parameters a result depends on (everything but the shard and how the work was scheduled)
string fingerprint(const TreeParams& p) {
  string f = tfm::format("base %d box %.17g %.17g %.17g %.17g depth %d safety %.17g m %d strata %d max_iter %d "
                         "seed %d first_newton %d center_max_iter %d center_first_newton %d center_hint_newton %d "
                         "center_max_period %d newton_max_period %d newton_iters %d newton_close2 %.17g "
                         "newton_repel2 %.17g newton_tol %.17g newton_margin %.17g prec %s rows %d leaf_stats %d "
                         "tiles %d flip_stats %d", p.base, p.x0, p.x1, p.y0, p.y1, p.depth, p.safety, p.m, p.strata, p.max_iter,
                         p.seed, p.first_newton, p.center_max_iter, p.center_first_newton, p.center_hint_newton,
                         p.center_max_period, p.newton_max_period, p.newton_iters, p.newton_close2, p.newton_repel2,
                         p.newton_tol, p.newton_margin, p.prec, p.rows, int(p.leaf_stats), p.tiles,
                         int(p.flip_stats));
  // (Only when on, so results saved before roulette existed still load)
  if (p.roulette_from) f += tfm::format(" roulette %d %d %d", p.roulette_from, p.roulette_log2, p.roulette_stride);
  f += " ks";
  for (const auto k : p.ks) f += tfm::format(" %d", k);
  return f;
}

template<class T> void write_vec(FILE* f, const char* name, const vector<T>& v) {
  fprintf(f, "%s %zu", name, v.size());
  for (const auto& x : v) {
    if constexpr (std::is_same_v<T, GroupSums>) fprintf(f, " %lld %lld %lld", (long long)x.s, (long long)x.q, (long long)x.p);
    else fprintf(f, " %lld", (long long)x);
  }
  fprintf(f, "\n");
}

template<class T> void read_vec(FILE* f, const char* name, vector<T>& v, const string& path) {
  char got[64];
  size_t n;
  slow_assert(fscanf(f, "%63s %zu", got, &n) == 2 && string(got) == name, "%s: expected %s", path, name);
  slow_assert(n == v.size(), "%s: %s has %d entries, expected %d", path, name, n, v.size());
  for (auto& x : v) {
    long long a, b, c;
    if constexpr (std::is_same_v<T, GroupSums>) {
      slow_assert(fscanf(f, "%lld %lld %lld", &a, &b, &c) == 3, "%s: short %s", path, name);
      x = GroupSums{a, b, c};
    } else {
      slow_assert(fscanf(f, "%lld", &a) == 1, "%s: short %s", path, name);
      x = T(a);
    }
  }
}

}  // namespace

void save_result(const TreeResult& R, const string& path) {
  FILE* f = fopen(path.c_str(), "w");
  slow_assert(f, "can't write %s", path);
  fprintf(f, "mandelbrot tree result 1\n%s\nshard %d %d\n", fingerprint(R.p).c_str(), R.p.shard, R.p.shards);
  write_vec(f, "certified", R.certified);
  write_vec(f, "exact", R.exact);
  write_vec(f, "area", R.area);
  write_vec(f, "diff", R.diff);
  write_vec(f, "float_area", R.float_area);
  write_vec(f, "delta", R.delta);
  write_vec(f, "alloc", R.alloc);
  write_vec(f, "tile_diff", R.tile_diff);
  write_vec(f, "tile_cert", R.tile_cert);
  write_vec(f, "flip_stats", R.flip_stats);
  write_vec(f, "counts", vector<int64_t>{R.leaves, R.centers, R.center_iters, R.leaf_iters, R.overflow, R.flips,
                                        R.batches});
  fprintf(f, "secs %.17g %.17g %.17g %.17g %.17g\n", R.tree_secs, R.center_kernel_secs, R.sample_secs,
          R.reduce_secs, R.secs);
  slow_assert(fclose(f) == 0, "error writing %s", path);
}

TreeResult load_result(const string& path, const TreeParams& p) {
  FILE* f = fopen(path.c_str(), "r");
  slow_assert(f, "can't read %s", path);
  char line[4096];
  slow_assert(fgets(line, sizeof(line), f) && string(line) == "mandelbrot tree result 1\n", "%s: not a tree result",
              path);
  slow_assert(fgets(line, sizeof(line), f), "%s: no parameters", path);
  const string want = fingerprint(p) + "\n";
  slow_assert(string(line) == want, "%s: parameters differ:\n  file: %s  run:  %s", path, line, want);
  TreeResult R = empty_result(p);
  slow_assert(fscanf(f, " shard %d %d", &R.p.shard, &R.p.shards) == 2, "%s: no shard", path);
  read_vec(f, "certified", R.certified, path);
  read_vec(f, "exact", R.exact, path);
  read_vec(f, "area", R.area, path);
  read_vec(f, "diff", R.diff, path);
  read_vec(f, "float_area", R.float_area, path);
  read_vec(f, "delta", R.delta, path);
  read_vec(f, "alloc", R.alloc, path);
  read_vec(f, "tile_diff", R.tile_diff, path);
  read_vec(f, "tile_cert", R.tile_cert, path);
  read_vec(f, "flip_stats", R.flip_stats, path);
  vector<int64_t> c(7);
  read_vec(f, "counts", c, path);
  R.leaves = c[0]; R.centers = c[1]; R.center_iters = c[2]; R.leaf_iters = c[3]; R.overflow = c[4]; R.flips = c[5];
  R.batches = c[6];
  slow_assert(fscanf(f, " secs %lf %lf %lf %lf %lf", &R.tree_secs, &R.center_kernel_secs, &R.sample_secs,
                     &R.reduce_secs, &R.secs) == 5, "%s: no secs", path);
  fclose(f);
  return R;
}

}  // namespace mandelbrot
