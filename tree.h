// Certified adaptive quadtree Monte Carlo for the areas of {c : g_M(c) < 2^-k}, on CPU threads or the GPU
//
// Only cells straddling the boundary of these sets contribute variance, and at fine scales they are a tiny
// fraction of all cells.  Starting from a base grid over a box (by default the upper half plane part of [-2, 0.5] × [-1.2, 1.2]),
// each level classifies its cells' centers with distance estimates (OrbitDE).  A cell that fits well inside
// its center's certified disk is decided exactly: interior cells count fully, and exterior cells count per
// threshold using Harnack bounds on g.  Otherwise it splits into 4 children, down to `depth` levels.
// Uncertified cells at the last level are leaves, sampled with m points each, in groups of strata^2 points
// jittered on a strata × strata grid.  Group means are iid, so each leaf gives an unbiased estimate with an
// unbiased variance estimate, and leaves are independent, so the total variance is the sum over leaves.
//
// Everything is computed from exact integer counts, and sample points are keyed by (leaf, sample), not by
// processing order.  So results are deterministic, and CPU and GPU runs agree except where host and device
// libm differ in the last bit.
//
// Besides each area A(k), we report differences D = A(k_i) - A(k_{i+1}) between consecutive thresholds with
// their own variances.  With prec = "compare", each leaf sample is classified with both float and double
// orbits, and we also report the paired difference A_float(k) - A_double(k), whose variance comes only from
// samples that flip, to measure the bias of low precision.
#pragma once

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>
namespace mandelbrot {

using std::string;
using std::vector;

struct TreeParams {
  int64_t base = 1000;       // Base grid size per axis
  double x0 = -2, x1 = 0.5, y0 = 0, y1 = 1.2;  // Domain box; estimates double it (conjugate symmetry), so a
                                               // box with y0 = 0 measures the symmetric region it reflects to
  int depth = 5;             // Refinement levels below the base grid
  double safety = 4;         // Certify a cell if its half-diagonal is at most dist / safety
  int m = 16;                // Samples per leaf
  int strata = 2;            // Each group of strata^2 samples is jittered on a strata × strata grid
  int64_t max_iter = 1 << 20;
  uint64_t seed = 1;
  vector<int64_t> ks;        // Thresholds 2^-k, increasing, at most 31
  int64_t first_newton = 8192;   // First Newton certificate attempt for leaf samples (and the GPU parking step)
  int64_t center_max_iter = 1 << 14;  // Iteration cap for cell centers (long orbits rarely certify a cell)
  int64_t center_first_newton = 8192; // First Newton interior certificate attempt for cell centers
  int64_t center_hint_newton = 0;     // Earlier, for cells whose parent's center was interior (0: no)
  int center_max_period = 256;        // Largest period Newton tries for cell centers
  int newton_max_period = 256;        // Largest period Newton tries for leaf samples
  int newton_iters = 30;             // Newton iterations per certificate attempt for leaf samples
  double newton_close2 = INFINITY;   // Leaf Newton gives up unless |f^p(w) - w|^2 < this after one iteration
  double newton_repel2 = 2;          // Leaf Newton gives up from its second iteration once |λ(w)|^2 exceeds this
  double newton_tol = 1e-10;         // Leaf Newton converges when |step| < this · |w| (-1: 1e-14)
  double newton_margin = 1e-6;       // Leaf Newton certifies when |λ|^2 < 1 - this (-1: 1e-9)
  int64_t burst = 128;               // Orbit steps per run call (finished lanes refill between bursts)
  int center_min_blocks = 2, sample_min_blocks = 3;  // GPU register budgets (H200: centers spill beyond 2)
  string prec = "double";    // Leaf orbits: double, dd (double-double), float, or compare (float), compareNN
                             // (double rounded to NN bits, NN ∈ {24, 27, 30, 36, 42, 48}), comparedd
                             // (double-double), each paired with double
  bool cuda = false;
  int64_t batch = 1 << 22;   // Target leaves per batch
  int64_t rows = -1;         // Only the first `rows` base rows (for tests), or -1 for all
  bool leaf_stats = false;   // Collect pilot-allocation statistics (TreeResult::alloc; needs m = 16)
  int tiles = 0;             // Resolve consecutive differences over a tiles × tiles grid of the box (0: off)
  int overlap = 1;           // Batches in flight at once on the GPU (host threads with their own streams)
  int shard = 0, shards = 1; // Only base rows r with r % shards == shard (for runs split across GPUs)
  bool flip_stats = false;   // With compare: classify the samples whose classifications differ (TreeResult::flip_stats)
  // Russian roulette for long leaf orbits (0: off): at every roulette_stride-th Newton step d ≥ roulette_from
  // (d = first_newton · 2^j + 8), an orbit not yet decided goes on with probability 2^-roulette_log2, weighted by
  // 2^roulette_log2 from then on, and stops otherwise.  Unbiased, with integer weights; thresholds must be at
  // least 8 steps from every decision step.  Keeping 1/2 every other octave (stride 2) balances cost against
  // variance: an octave's escapes fall like 2^-j, so its added variance falls like 2^-j/2, as does its cost.
  int64_t max_level_cells = int64_t(1) << 31;  // Batches split while a tree level reaches this (tests lower it)
  double progress = 0;       // If positive, print progress at most this often (seconds)
  string dump;               // If set, append each leaf's cell and sample outcome codes here (see dump_leaves)
  int64_t roulette_from = 0;
  int roulette_log2 = 1;
  int roulette_stride = 1;
};

// Pilot-allocation statistics: leaves are classed by c1, the count below threshold k among their first P
// samples (P = kPilots[pi]), and each class accumulates, over its leaves, N (leaves), S = Σ c2 (8 - c2) with
// c2 the count among samples 8..15 (S / 56 is unbiased for Σ p (1 - p) for iid samples), W2 = iterations
// of samples 8..15, and Wp = iterations of the first P samples.  Measuring on samples the class did not
// select makes these unbiased for what a second phase would see.
const int kPilots[3] = {2, 4, 8};
constexpr int64_t alloc_index(const int k, const int pi, const int c1, const int field) {
  return ((int64_t(k) * 3 + pi) * 9 + c1) * 4 + field;
}

// Exact sums over leaves l and sample groups g of the group sums c_g of a per-sample value: Σ c_g, Σ c_g^2,
// and Σ_l (Σ_g c_g)^2.  The leaf estimate and its variance follow from these (see TreeResult).
struct GroupSums {
  int64_t s = 0, q = 0, p = 0;
  void operator+=(const GroupSums& o) { s += o.s; q += o.q; p += o.p; }
};

struct TreeResult {
  TreeParams p;
  vector<int64_t> certified;          // [d * K + k]: certified cells at depth d below threshold k
  vector<int64_t> exact;              // Certified cells per depth
  vector<GroupSums> area, diff;       // Leaf sums per threshold for areas, and for consecutive differences
  vector<GroupSums> float_area, delta;  // With compare: float areas, and float - double
  vector<int64_t> alloc;              // With leaf_stats: [alloc_index(k, pi, c1, {N, S, W2, Wp})]
  vector<GroupSums> tile_diff;        // With tiles: [tile * K + k] leaf sums of the k, k + 1 difference
  vector<int64_t> tile_cert;          // With tiles: [tile * K + k] area certified below k, in finest cells
  vector<int64_t> flip_stats;         // With flip_stats: see FlipChunk in tree.cc (flip_stats_size(K) entries)
  int64_t leaves = 0, centers = 0, center_iters = 0, leaf_iters = 0, overflow = 0, flips = 0, batches = 0;
  double tree_secs = 0, center_kernel_secs = 0, sample_secs = 0, reduce_secs = 0, secs = 0;

  // Areas over the whole plane (twice the upper half), consecutive differences, float - double, and
  // variances of each
  double cell_area(int d) const;
  double estimate(const vector<GroupSums>& sums, int k, bool with_certified, int sign = 1) const;
  double variance(const vector<GroupSums>& sums, int k) const;
  double area_estimate(int k) const { return estimate(area, k, true); }
  double diff_estimate(int k) const;
  // Per tile (index ty * tiles + tx): consecutive difference A(k) - A(k + 1) and its variance
  double tile_diff_estimate(int tile, int k) const;
  double tile_diff_variance(int tile, int k) const;
};

TreeResult run_tree(const TreeParams& p);

// Size of TreeResult::flip_stats for K thresholds
constexpr int flip_stats_size(const int K) { return 9 + 9 * K + 33 + 34 * 2; }

// Sharded runs: a shard's result saved to a file (its exact sums, and the parameters that determine them),
// loaded back (checking that the parameters match p, up to the shard), and merged.  Merging every shard of a
// run gives exactly the unsharded result, since all sums are integers and samples are keyed by cell.
void save_result(const TreeResult& R, const string& path);
TreeResult load_result(const string& path, const TreeParams& p);
void merge(TreeResult& R, const TreeResult& S);
TreeResult empty_result(const TreeParams& p);

}  // namespace mandelbrot
