// Böttcher octave energy near hyperbolic component roots, by period
//
// Tests whether the whole-circle octave energy S_j is a sum of root-local contributions with one profile
// per root type, S_j ≈ Σ_p w_p Ψ(j/p): for each octave j, the energy within ±W bins (bins of 2^-(j+1) in
// angle, θ and 1-θ folded) of the two landing angles of every root of period p ≤ max_period, in increasing
// period so that each bin counts once, split into satellite and primitive roots.  Only periods p ≤ j - gap
// are used, so that the windows cover a small fraction of the circle (reported as coverage).  Writes lines
//   j p satellite_energy primitive_energy
// to stdout after the summary, with p = 0 for the whole circle S_j.

#include "angles.h"
#include "debug.h"
#include "numpy.h"
#include "octaves.h"
#include "print.h"
#include "wall_time.h"
#include <algorithm>
#include <cmath>
#include <vector>
namespace mandelbrot {
namespace {

using std::max;
using std::min;
using std::vector;

void run(const string& path, const int min_j, const int max_j, const int max_period, const int W, const int gap) {
  auto t0 = wall_time();
  const auto F = read_numpy(path);
  slow_assert(F.shape.size() == 2 && F.shape[1] == 2 && F.shape[0] >= int64_t(2) << max_j, "bad coefficient file");
  const auto roots = lavaurs(max_period);
  print("read %s %s, %d roots of period ≤ %d: %.2f s", path, F.shape, roots.size(), max_period,
        (wall_time() - t0).seconds());

  // Roots by period, with the cusp 1/4 (angle 0) as the period-1 primitive root
  vector<vector<const Root*>> by_p(max_period + 1);
  for (const auto& r : roots) by_p[r.w.lo.q].push_back(&r);
  vector<string> lines;
  for (int j = min_j; j <= max_j; j++) {
    const auto t = wall_time();
    const auto e = octave_energy(F.data, j);
    const int64_t lo = int64_t(1) << j, n = 2 * lo;
    vector<char> used(lo + 1, 0);
    double S = 0;
    for (const double v : e) S += v;
    const auto window = [&](const double theta) {
      const double t = min(theta, 1 - theta);
      const int64_t c = std::llround(t * double(n));
      double s = 0;
      for (int64_t i = max(int64_t(0), c - W); i <= min(lo, c + W); i++)
        if (!used[i]) { used[i] = 1; s += e[i]; }
      return s;
    };
    double claimed = 0;
    lines.push_back(tfm::format("%d 0 %.10e 0", j, S));
    const double cusp = window(0);
    claimed += cusp;
    lines.push_back(tfm::format("%d 1 0 %.10e", j, cusp));
    for (int p = 2; p <= min(max_period, j - gap); p++) {
      double sat = 0, prim = 0;
      for (const Root* r : by_p[p]) {
        const double s = window(r->w.lo.value()) + window(r->w.hi.value());
        (r->satellite ? sat : prim) += s;
      }
      claimed += sat + prim;
      lines.push_back(tfm::format("%d %d %.10e %.10e", j, p, sat, prim));
    }
    int64_t covered = 0;
    for (const char u : used) covered += u;
    print("j %2d: S_j %.6e, roots of period ≤ %d hold %.1f%% on %.1f%% of angles, %.2f s", j, S,
          min(max_period, j - gap), 100 * claimed / S, 100.0 * covered / (lo + 1), (wall_time() - t).seconds());
  }
  print("# j p satellite primitive");
  for (const auto& l : lines) print("%s", l);
}

}  // namespace
}  // namespace mandelbrot

int main(const int argc, const char** argv) {
  using namespace mandelbrot;
  try {
    slow_assert(argc >= 2, "usage: %s f-kK.npy [min_j] [max_j] [max_period] [window] [gap]", argv[0]);
    const int min_j = argc > 2 ? atoi(argv[2]) : 8, max_j = argc > 3 ? atoi(argv[3]) : 26;
    const int max_period = argc > 4 ? atoi(argv[4]) : 20, W = argc > 5 ? atoi(argv[5]) : 4;
    const int gap = argc > 6 ? atoi(argv[6]) : 7;
    run(argv[1], min_j, max_j, max_period, W, gap);
    return 0;
  } catch (const std::exception& e) {
    die(e.what());
  }
}
