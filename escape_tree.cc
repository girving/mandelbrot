// Certified adaptive quadtree Monte Carlo for the areas of {c : g_M(c) < 2^-k} (see tree.h)

#include "argparse.hpp"
#include "debug.h"
#include "engine.h"
#include "print.h"
#include "tree.h"
#include <cmath>
namespace mandelbrot {
namespace {

void report(const TreeResult& R) {
  const auto& p = R.p;
  const int K = p.ks.size();
  if (p.x0 != -2 || p.x1 != 0.5 || p.y0 != 0 || p.y1 != 1.2)
    print("box [%.17g, %.17g] × [%.17g, %.17g] (doubled)", p.x0, p.x1, p.y0, p.y1);
  print("base %d, depth %d (effective grid %d), safety %g, %d samples/leaf, strata %d, max_iter %d, seed %d, "
        "first Newton %d, center max_iter %d, center first Newton %d, center max period %d, Newton max period %d, Newton iterations %d, Newton close %g, burst %d, prec %s, %s, %d threads: %.1f s (tree %.1f s with centers %.1f s, sampling %.1f s, "
        "reduce %.1f s, %d batches)", p.base, p.depth, p.base << p.depth, p.safety, p.m, p.strata, p.max_iter, p.seed,
        p.first_newton, p.center_max_iter, p.center_first_newton, p.center_max_period, p.newton_max_period, p.newton_iters, std::sqrt(p.newton_close2), p.burst, p.prec, p.cuda ? "cuda" : "cpu", cpu_threads(), R.secs, R.tree_secs, R.center_kernel_secs, R.sample_secs,
        R.reduce_secs, R.batches);
  print("  sampling throughput: %.3g iterations/s", double(R.leaf_iters) / R.sample_secs);
  print("  centers: %.3g cells, %.3g iterations; leaves: %.3g leaves, %.3g samples, %.3g iterations; "
        "%d orbits overflowed", double(R.centers), double(R.center_iters), double(R.leaves),
        double(R.leaves * p.m), double(R.leaf_iters), R.overflow);
  string es = "  certified cells by depth:";
  for (const auto n : R.exact) es += tfm::format(" %.3g", double(n));
  print(es);
  const bool compare = p.prec.starts_with("compare");
  print("  %s:", compare ? "double" : p.prec);
  print("       k      area{g < 2^-k}      std err");
  for (int k = 0; k < K; k++)
    print("    %8d   %.10f   %.2e", p.ks[k], R.area_estimate(k), std::sqrt(R.variance(R.area, k)));
  print("       k → k'           A(k) - A(k')     std err");
  for (int k = 0; k + 1 < K; k++)
    print("    %8d → %-8d  %.6e   %.2e", p.ks[k], p.ks[k + 1], R.diff_estimate(k), std::sqrt(R.variance(R.diff, k)));
  if (compare) {
    print("  %s:", p.prec == "compare" ? "float" : p.prec == "comparedd" ? string("double-double")
                                                 : "rounded to " + p.prec.substr(7) + " bits");
    for (int k = 0; k < K; k++)
      print("    %8d   %.10f   %.2e", p.ks[k], R.estimate(R.float_area, k, true),
            std::sqrt(R.variance(R.float_area, k)));
    print("  %s - double (%d flipped samples of %d):", p.prec == "compare" ? "float" : "rounded", R.flips,
          R.leaves * p.m);
    for (int k = 0; k < K; k++)
      print("    %8d   %+.3e   %.2e", p.ks[k], R.estimate(R.delta, k, false), std::sqrt(R.variance(R.delta, k)));
  }
}

// Two-phase allocation: a pilot of P samples classes each leaf, and the second phase gives class c n_c ≥ 1
// fresh samples per leaf, estimating each leaf from its second-phase samples only (so unbiased).  With class
// variance σ_c^2 and per-sample cost κ_c measured on samples 8..15, minimize Σ N_c n_c κ_c + pilot cost
// subject to Σ N_c σ_c^2 / n_c equal to the variance of uniform sampling with m = 16: the unconstrained
// optimum is n_c = λ σ_c / sqrt(κ_c) (Neyman), clamped below at 1.  Reports sampling cost relative to uniform.
void allocation_report(const TreeResult& R) {
  const int K = R.p.ks.size();
  print("  two-phase allocation (leaf sampling cost at equal variance, relative to uniform m = 16):");
  print("       k    pilot 2   pilot 4   pilot 8   (pilot 8, without pilot cost)");
  for (int k = 0; k < K; k++) {
    string line = tfm::format("    %8d", R.p.ks[k]);
    double free8 = 0;
    for (int pi = 0; pi < 3; pi++) {
      const int P = kPilots[pi];
      double Vu = 0, Bu = 0, Bp = 0;
      vector<double> N, s2, kappa;
      for (int c = 0; c <= P; c++) {
        const auto a = [&](int f) { return double(R.alloc[alloc_index(k, pi, c, f)]); };
        if (!a(0)) continue;
        Vu += a(1) / 56 / 16; Bu += 2 * a(2); Bp += a(3);
        N.push_back(a(0)); s2.push_back(a(1) / (56 * a(0))); kappa.push_back(std::max(1.0, a(2) / (8 * a(0))));
      }
      const auto solve = [&](const double lambda, double& V, double& B) {
        V = 0; B = 0;
        for (size_t c = 0; c < N.size(); c++) {
          const double n = std::max(1.0, lambda * std::sqrt(s2[c] / kappa[c]));
          V += N[c] * s2[c] / n; B += N[c] * n * kappa[c];
        }
      };
      double lo = 1e-12, hi = 1e12, V, B;
      for (int it = 0; it < 200; it++) {
        const double mid = std::sqrt(lo * hi);
        solve(mid, V, B);
        (V > Vu ? lo : hi) = mid;
      }
      solve(hi, V, B);
      line += tfm::format("   %7.3f", (B + Bp) / Bu);
      if (P == 8) free8 = B / Bu;
    }
    print("%s   (%.3f)", line, free8);
  }
  // Class table for the last threshold with an 8-sample pilot
  const int k = K - 1;
  double tN = 0, tS = 0, tW = 0;
  for (int c = 0; c <= 8; c++) {
    tN += R.alloc[alloc_index(k, 2, c, 0)]; tS += R.alloc[alloc_index(k, 2, c, 1)]; tW += R.alloc[alloc_index(k, 2, c, 2)];
  }
  print("    k = %d, pilot 8: class (pilot count below), share of leaves, of variance, of cost", R.p.ks[k]);
  for (int c = 0; c <= 8; c++)
    print("      %d   %6.2f%%   %6.2f%%   %6.2f%%", c, 100 * R.alloc[alloc_index(k, 2, c, 0)] / tN,
          100 * R.alloc[alloc_index(k, 2, c, 1)] / tS, 100 * R.alloc[alloc_index(k, 2, c, 2)] / tW);
}

}  // namespace
}  // namespace mandelbrot

int main(const int argc, const char** argv) {
  using namespace mandelbrot;
  try {
    argparse::ArgumentParser program("escape_tree");
    program.add_argument("ks").help("thresholds g < 2^-k, increasing (options go before these)").remaining()
        .scan<'i', int>();
    program.add_argument("--base").help("base grid size per axis").scan<'i', int64_t>().default_value(int64_t(1000));
    program.add_argument("--depth").help("refinement levels").scan<'i', int>().default_value(5);
    program.add_argument("--safety").help("certify if half-diagonal ≤ dist / safety").scan<'g', double>()
        .default_value(4.0);
    program.add_argument("--m").help("samples per leaf").scan<'i', int>().default_value(16);
    program.add_argument("--strata").help("jitter groups of strata^2 samples").scan<'i', int>().default_value(2);
    program.add_argument("--max-iter").scan<'i', int64_t>().default_value(int64_t(1) << 20);
    program.add_argument("--seed").scan<'i', int64_t>().default_value(int64_t(1));
    program.add_argument("--prec").help("leaf orbit precision: double, float, compare (float vs double), or "
                                        "compareNN (double rounded to NN ∈ {30, 36, 42, 48} bits vs double), or comparedd (double-double vs double)")
        .default_value(string("double"));
    program.add_argument("--cuda").help("run on the GPU").default_value(false).implicit_value(true);
    program.add_argument("--batch").help("target leaves per batch (0: 2^26 on the GPU, 2^22 on the CPU)")
        .scan<'i', int64_t>().default_value(int64_t(0));
    program.add_argument("--first-newton").help("first Newton certificate attempt for leaf samples")
        .scan<'i', int64_t>().default_value(int64_t(16384));
    program.add_argument("--center-max-iter").help("iteration cap for cell centers").scan<'i', int64_t>()
        .default_value(int64_t(1) << 14);
    program.add_argument("--newton-max-period").help("largest period Newton tries for leaf samples")
        .scan<'i', int>().default_value(1024);
    program.add_argument("--newton-iters").help("Newton iterations per certificate attempt").scan<'i', int>()
        .default_value(30);
    program.add_argument("--center-first-newton").help("first Newton attempt for cell centers")
        .scan<'i', int64_t>().default_value(int64_t(8192));
    program.add_argument("--newton-close").help("leaf Newton needs |f^p(w) - w| below this after one iteration")
        .scan<'g', double>().default_value(double(INFINITY));
    program.add_argument("--center-max-period").help("largest period Newton tries for cell centers")
        .scan<'i', int>().default_value(256);
    program.add_argument("--newton-tol").help("leaf Newton converged when |step| < this (relative); -1: 1e-14")
        .scan<'g', double>().default_value(1e-10);
    program.add_argument("--newton-margin").help("leaf Newton certifies when |λ|^2 < 1 - this; -1: 1e-9")
        .scan<'g', double>().default_value(1e-6);
    program.add_argument("--box").help("domain x0 x1 y0 y1 (estimates double it by conjugate symmetry)")
        .nargs(4).scan<'g', double>().default_value(vector<double>{-2, 0.5, 0, 1.2});
    program.add_argument("--leaf-stats").help("report two-phase allocation gains (needs --m 16)")
        .default_value(false).implicit_value(true);
    program.add_argument("--burst").help("orbit steps per run call").scan<'i', int64_t>().default_value(int64_t(64));
    program.parse_args(argc, argv);

    TreeParams p;
    p.base = program.get<int64_t>("--base");
    p.depth = program.get<int>("--depth");
    p.safety = program.get<double>("--safety");
    p.m = program.get<int>("--m");
    p.strata = program.get<int>("--strata");
    p.max_iter = program.get<int64_t>("--max-iter");
    p.seed = uint64_t(program.get<int64_t>("--seed"));
    p.prec = program.get<string>("--prec");
    p.cuda = program.get<bool>("--cuda");
    p.batch = program.get<int64_t>("--batch");
    if (!p.batch) p.batch = int64_t(1) << (p.cuda ? 26 : 22);
    p.first_newton = program.get<int64_t>("--first-newton");
    p.center_max_iter = program.get<int64_t>("--center-max-iter");
    p.newton_max_period = program.get<int>("--newton-max-period");
    p.burst = program.get<int64_t>("--burst");
    p.center_first_newton = program.get<int64_t>("--center-first-newton");
    p.center_max_period = program.get<int>("--center-max-period");
    p.newton_iters = program.get<int>("--newton-iters");
    p.newton_tol = program.get<double>("--newton-tol");
    p.newton_margin = program.get<double>("--newton-margin");
    p.newton_close2 = program.get<double>("--newton-close") * program.get<double>("--newton-close");
    p.ks = program.get<vector<int>>("ks");
    for (const int k : p.ks) slow_assert(k + 8 <= p.max_iter, "need max_iter ≥ k + 8 for k = %d", k);
    for (size_t i = 0; i + 1 < p.ks.size(); i++) slow_assert(p.ks[i] < p.ks[i + 1], "thresholds must increase");
    slow_assert(p.first_newton >= 1, "need --first-newton ≥ 1");
    p.leaf_stats = program.get<bool>("--leaf-stats");
    const auto box = program.get<vector<double>>("--box");
    p.x0 = box[0]; p.x1 = box[1]; p.y0 = box[2]; p.y1 = box[3];
    slow_assert(p.x0 < p.x1 && p.y0 < p.y1, "empty box");
    const auto R = run_tree(p);
    report(R);
    if (p.leaf_stats) allocation_report(R);
    return 0;
  } catch (const std::exception& e) {
    die(e.what());
  }
}
