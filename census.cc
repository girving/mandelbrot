// Census of Böttcher octave energy near the roots of hyperbolic components
//
// For each octave j, assigns each angle bin within ±W bins of a root ray to that root (lowest period first,
// so no bin counts twice), and tabulates the energy near roots by period and type (satellite/primitive).

#include "angles.h"
#include "arith.h"
#include "debug.h"
#include "numpy.h"
#include "octaves.h"
#include "print.h"
#include "wall_time.h"
#include <algorithm>
#include <cmath>
namespace mandelbrot {
namespace {

using std::max;
using std::min;

void run(const string& path, const int max_j, const int P, const int W) {
  const auto t0 = wall_time();
  const auto F = read_numpy(path);
  slow_assert(F.shape.size() == 2 && F.shape[1] == 2 && F.shape[0] >= int64_t(2) << max_j, "bad coefficient file");
  auto roots = lavaurs(P);
  print("read %s, %d roots of period 2..%d (%.2f s)", path, roots.size(), P, (wall_time() - t0).seconds());

  // near[type][p][j]: energy within ±W bins of period p roots; type 0 = satellite, 1 = primitive.  p = 1 is the cusp.
  vector<vector<vector<double>>> near(2, vector<vector<double>>(P + 1, vector<double>(max_j + 1)));
  vector<double> S(max_j + 1);
  for (int j = 1; j <= max_j; j++) {
    const auto e = octave_energy(F.data, j);
    const int64_t lo = int64_t(1) << j, n = 2*lo;
    for (const double v : e) S[j] += v;
    vector<bool> used(lo + 1);
    const auto take = [&](const double t, double& acc) {
      const int64_t c = std::llround(min(t, 1 - t) * n);
      for (int64_t i = max(int64_t(0), c - W); i <= min(lo, c + W); i++)
        if (!used[i]) { used[i] = true; acc += e[i]; }
    };
    take(0, near[1][1][j]);  // The cardioid's cusp is primitive
    for (const auto& r : roots) {  // lavaurs returns roots in increasing period
      auto& acc = near[r.satellite ? 0 : 1][r.w.lo.q][j];
      take(r.w.lo.value(), acc);
      take(r.w.hi.value(), acc);
    }
  }
  print("octaves done (%.2f s)\n", (wall_time() - t0).seconds());

  // Fraction of each octave's energy near roots, by period and type
  for (int type = 0; type < 2; type++) {
    print("%s roots: percent of S_j within ±%d bins, by root period p (rows) and octave j (columns)",
          type ? "Primitive" : "Satellite", W);
    string h = "    p  ";
    for (int j = 10; j <= max_j; j++) h += tfm::format(" %5d", j);
    print(h);
    for (int p = 1; p <= P; p++) {
      if (!near[type][p][max_j]) continue;
      string row = tfm::format("   %2d  ", p);
      for (int j = 10; j <= max_j; j++) row += tfm::format(" %5.2f", 100 * near[type][p][j] / S[j]);
      print(row);
    }
  }
  string tot = "  all  ";
  for (int j = 10; j <= max_j; j++) {
    double s = 0;
    for (int type = 0; type < 2; type++) for (int p = 1; p <= P; p++) s += near[type][p][j];
    tot += tfm::format(" %5.1f", 100 * s / S[j]);
  }
  print("Total percent of S_j near roots of period <= %d:\n%s", P, tot);

  // Local log-power exponent per period class over j = 18..max_j
  print("\nlog-power exponent gamma of near-root energy, j = 18..%d (resolved roots only: j >= p + 6)", max_j);
  for (int type = 0; type < 2; type++)
    for (int p = 1; p <= P; p++) {
      vector<double> x, y;
      for (int j = max(18, p + 6); j <= max_j; j++)
        if (near[type][p][j] > 0) { x.push_back(std::log(j)); y.push_back(std::log(near[type][p][j])); }
      if (x.size() >= 4)
        print("   %s p = %2d: gamma %5.2f over %d octaves", type ? "primitive" : "satellite", p, -fit_line(x, y).b,
              int(x.size()));
    }
  print("\ntotal %.1f s", (wall_time() - t0).seconds());
}

}  // namespace
}  // namespace mandelbrot

int main(const int argc, const char** argv) {
  using namespace mandelbrot;
  try {
    slow_assert(argc == 5, "usage: %s <f-k27.npy> <max_j> <max_period> <window_bins>", argv[0]);
    run(argv[1], atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
    return 0;
  } catch (const std::exception& e) {
    die(e.what());
  }
}
