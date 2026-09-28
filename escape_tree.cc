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
  print("base %d, depth %d (effective grid %d), safety %g, %d samples/leaf, strata %d, max_iter %d, seed %d, "
        "first Newton %d, center max_iter %d, center first Newton %d, Newton max period %d, Newton iterations %d, burst %d, prec %s, %s, %d threads: %.1f s (tree %.1f s with centers %.1f s, sampling %.1f s, "
        "reduce %.1f s, %d batches)", p.base, p.depth, p.base << p.depth, p.safety, p.m, p.strata, p.max_iter, p.seed,
        p.first_newton, p.center_max_iter, p.center_first_newton, p.newton_max_period, p.newton_iters, p.burst, p.prec, p.cuda ? "cuda" : "cpu", cpu_threads(), R.secs, R.tree_secs, R.center_kernel_secs, R.sample_secs,
        R.reduce_secs, R.batches);
  print("  sampling throughput: %.3g iterations/s", double(R.leaf_iters) / R.sample_secs);
  print("  centers: %.3g cells, %.3g iterations; leaves: %.3g leaves, %.3g samples, %.3g iterations; "
        "%d orbits overflowed", double(R.centers), double(R.center_iters), double(R.leaves),
        double(R.leaves * p.m), double(R.leaf_iters), R.overflow);
  string es = "  certified cells by depth:";
  for (const auto n : R.exact) es += tfm::format(" %.3g", double(n));
  print(es);
  const bool compare = p.prec == "compare";
  print("  %s:", compare ? "double" : p.prec);
  print("       k      area{g < 2^-k}      std err");
  for (int k = 0; k < K; k++)
    print("    %8d   %.10f   %.2e", p.ks[k], R.area_estimate(k), std::sqrt(R.variance(R.area, k)));
  print("       k → k'           A(k) - A(k')     std err");
  for (int k = 0; k + 1 < K; k++)
    print("    %8d → %-8d  %.6e   %.2e", p.ks[k], p.ks[k + 1], R.diff_estimate(k), std::sqrt(R.variance(R.diff, k)));
  if (compare) {
    print("  float:");
    for (int k = 0; k < K; k++)
      print("    %8d   %.10f   %.2e", p.ks[k], R.estimate(R.float_area, k, true),
            std::sqrt(R.variance(R.float_area, k)));
    print("  float - double (%d flipped samples of %d):", R.flips, R.leaves * p.m);
    for (int k = 0; k < K; k++)
      print("    %8d   %+.3e   %.2e", p.ks[k], R.estimate(R.delta, k, false), std::sqrt(R.variance(R.delta, k)));
  }
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
    program.add_argument("--prec").help("leaf orbit precision: double, float, or compare")
        .default_value(string("double"));
    program.add_argument("--cuda").help("run on the GPU").default_value(false).implicit_value(true);
    program.add_argument("--batch").help("target leaves per batch (0: 2^26 on the GPU, 2^22 on the CPU)")
        .scan<'i', int64_t>().default_value(int64_t(0));
    program.add_argument("--first-newton").help("first Newton certificate attempt for leaf samples")
        .scan<'i', int64_t>().default_value(int64_t(16384));
    program.add_argument("--center-max-iter").help("iteration cap for cell centers").scan<'i', int64_t>()
        .default_value(int64_t(1) << 14);
    program.add_argument("--newton-max-period").help("largest period Newton tries for leaf samples")
        .scan<'i', int>().default_value(256);
    program.add_argument("--newton-iters").help("Newton iterations per certificate attempt").scan<'i', int>()
        .default_value(30);
    program.add_argument("--center-first-newton").help("first Newton attempt for cell centers")
        .scan<'i', int64_t>().default_value(int64_t(8192));
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
    p.newton_iters = program.get<int>("--newton-iters");
    p.ks = program.get<vector<int>>("ks");
    for (const int k : p.ks) slow_assert(k + 8 <= p.max_iter, "need max_iter ≥ k + 8 for k = %d", k);
    for (size_t i = 0; i + 1 < p.ks.size(); i++) slow_assert(p.ks[i] < p.ks[i + 1], "thresholds must increase");
    slow_assert(p.first_newton >= 1, "need --first-newton ≥ 1");
    report(run_tree(p));
    return 0;
  } catch (const std::exception& e) {
    die(e.what());
  }
}
