// Full-resolution scan of Böttcher octave energy over external angle
//
// For each octave j (coefficients 2^j <= m < 2^(j+1)), computes the energy density over angle at resolution
// 2^-(j+1) and measures:
//   1. Energy within ±W bins of parabolic points (predicted ~ j^-4 at satellite roots, j^-6 at primitive
//      cusps, since ψ has 1/log and 1/log^2 singularities there) and Misiurewicz points (predicted to decay
//      like a power of n instead).
//   2. Inverse tuning: energy on the angles of a tuned copy of M, decoded to M's angles at the resolution
//      octave j supports, compared with the full distribution of earlier octaves j'.

#include "angles.h"
#include "arith.h"
#include "debug.h"
#include "numpy.h"
#include "octaves.h"
#include "print.h"
#include "wall_time.h"
#include <algorithm>
#include <cmath>
#include <functional>
namespace mandelbrot {
namespace {

using std::function;
using std::max;
using std::min;

struct Point { string name; double t; string type; };

// Sum e over bins [lo, hi) into K equal groups (K divides the bin count)
vector<double> coarsen(span<const double> e, const int64_t K) {
  const int64_t r = int64_t(e.size()) / K;
  slow_assert(r * K == int64_t(e.size()), "coarsen: %d bins into %d", e.size(), K);
  vector<double> c(K);
  for (int64_t i = 0; i < int64_t(e.size()); i++) c[i / r] += e[i];
  return c;
}

double tv(span<const double> P, span<const double> Q) {
  double sp = 0, sq = 0, d = 0;
  for (size_t i = 0; i < P.size(); i++) { sp += P[i]; sq += Q[i]; }
  for (size_t i = 0; i < P.size(); i++) d += std::abs(P[i] / sp - Q[i] / sq);
  return d / 2;
}

void run(const string& path, const int max_j, const string& compare) {
  auto t0 = wall_time();
  const auto F = read_numpy(path);
  slow_assert(F.shape.size() == 2 && F.shape[1] == 2 && F.shape[0] >= int64_t(2) << max_j, "bad coefficient file");
  print("read %s %s: %.2f s", path, F.shape, (wall_time() - t0).seconds());

  const vector<Point> points = {
      {"cusp 1/4", 0, "primitive"},       {"period-3 root -1.75", 3./7, "primitive"},
      {"period-4 root in 1/3 limb", 1./5, "primitive"}, {"period-4 root -1.94", 7./15, "primitive"},
      {"root -3/4 (1/2 bulb)", 1./3, "satellite"}, {"1/3 bulb root", 1./7, "satellite"},
      {"1/4 bulb root", 1./15, "satellite"},  {"1/5 bulb root", 1./31, "satellite"},
      {"2/5 bulb root", 9./31, "satellite"},  {"root -5/4 (period 4)", 2./5, "satellite"},
      {"tip -2", 0.5, "Misiurewicz"},          {"c = i", 1./6, "Misiurewicz"},
      {"angle 1/4", 0.25, "Misiurewicz"}};
  const vector<int> windows = {4, 32};
  vector<vector<vector<double>>> near(windows.size(), vector<vector<double>>(points.size(), vector<double>(max_j + 1)));

  struct Copy { string name; Wake w; };
  const vector<Copy> copies = {{"period-2 copy", cardioid_wake(1, 2)}, {"period-3 copy (1/3 bulb)", cardioid_wake(1, 3)}};
  const int Kmax = 1 << 12;             // half-circle bins kept per octave for comparisons
  vector<vector<double>> coarse(max_j + 1);  // coarse[j] = folded map with min(2^j, Kmax) bins
  vector<vector<vector<double>>> decoded(copies.size(), vector<vector<double>>(max_j + 1));
  vector<vector<double>> captured(copies.size(), vector<double>(max_j + 1));
  vector<vector<int>> Ds(copies.size(), vector<int>(max_j + 1));
  vector<double> S(max_j + 1);

  // Optional cross-check against an existing half-circle map (e.g. from numpy)
  Numpy other;
  if (compare.size()) other = read_numpy(compare);

  double parseval = 0;
  print("\n j   time(s)   S_j          ");
  for (int j = 1; j <= max_j; j++) {
    const auto t = wall_time();
    const auto e = octave_energy(F.data, j);
    const int64_t lo = int64_t(1) << j, n = 2*lo;
    double sum = 0, direct = 0;
    for (const double v : e) sum += v;
    for (int64_t m = lo; m < n; m++) direct += (m - 1) * sqr(F.data[2*m] + F.data[2*m+1]);
    S[j] = sum;
    parseval = max(parseval, std::abs(sum / direct - 1));

    // Folded map on 2^j half-circle bins: put the θ = 1/2 entry into the last bin
    vector<double> half(e.begin(), e.end() - 1);
    half.back() += e.back();
    coarse[j] = lo <= Kmax ? half : coarsen(half, Kmax);

    // 1. Parabolic and Misiurewicz windows
    for (size_t w = 0; w < windows.size(); w++)
      for (size_t p = 0; p < points.size(); p++) {
        const double t = min(points[p].t, 1 - points[p].t);
        const int64_t c = std::llround(t * n);
        double s = 0;
        for (int64_t i = max(int64_t(0), c - windows[w]); i <= min(lo, c + windows[w]); i++) s += e[i];
        near[w][p][j] = s;
      }

    // 2. Inverse tuning at the resolution this octave supports
    for (size_t k = 0; k < copies.size(); k++) {
      const auto& w = copies[k].w;
      const int r = w.lo.q, D = max(1, (j + 1) / r - 1);
      Ds[k][j] = D;
      vector<double> P(int64_t(1) << (D - 1));
      double wake = 0, got = 0;
      for (int64_t i = 0; i < lo; i++) {
        if (!w.contains(double(i) / n)) continue;
        wake += e[i];
        uint64_t prefix;
        if (untune(w, uint64_t(i), j + 1, D, prefix) == D) {
          const uint64_t h = prefix < (uint64_t(1) << (D - 1)) ? prefix : (uint64_t(1) << D) - 1 - prefix;
          P[h] += e[i];
          got += e[i];
        }
      }
      decoded[k][j] = P;
      captured[k][j] = got / wake;
    }

    string check;
    if (other.shape.size() == 2 && j < other.shape[0] && lo >= other.shape[1]) {
      const auto mine = coarsen(half, other.shape[1]);
      double d = 0;
      for (int64_t b = 0; b < other.shape[1]; b++) d += std::abs(mine[b] - other.data[j * other.shape[1] + b]);
      check = tfm::format("   vs numpy map: L1 diff / S_j = %.1e", d / sum);
    }
    print("%2d   %6.2f   %.6e%s", j, (wall_time() - t).seconds(), sum, check);
  }
  print("max Parseval relative error %.1e", parseval);

  // 1. Log-power and n-power exponents near each point
  const int j0 = 14, jm = min(20, max_j), j1 = max_j;
  const auto exponents = [&](const vector<double>& E, const int a, const int b, double& gamma, double& alpha) {
    vector<double> lj, x, y;
    for (int j = a; j <= b; j++) { lj.push_back(std::log(j)); x.push_back(j * std::log(2.)); y.push_back(std::log(E[j])); }
    gamma = -fit_line(lj, y).b;
    alpha = -fit_line(x, y).b;
  };
  for (size_t w = 0; w < windows.size(); w++) {
    print("\n1. Energy within ±%d bins (±%d·2^-(j+1)) of each point: fitted exponents", windows[w], windows[w]);
    print("   %-28s %-12s  gamma %2d-%2d  gamma %2d-%2d  gamma %2d-%2d  alpha %2d-%2d  alpha %2d-%2d   share at j=%d",
          "point", "type", j0, j1, j0, jm, jm, j1, j0, jm, jm, j1, j1);
    for (size_t p = 0; p < points.size(); p++) {
      double g, a, g1, a1, g2, a2;
      exponents(near[w][p], j0, j1, g, a);
      exponents(near[w][p], j0, jm, g1, a1);
      exponents(near[w][p], jm, j1, g2, a2);
      print("   %-28s %-12s  %11.2f  %11.2f  %11.2f  %11.3f  %11.3f   %9.2e", points[p].name, points[p].type,
            g, g1, g2, a1, a2, near[w][p][j1] / S[j1]);
    }
  }
  {
    double g, a, g1, a1, g2, a2;
    exponents(S, j0, j1, g, a); exponents(S, j0, jm, g1, a1); exponents(S, jm, j1, g2, a2);
    print("   %-28s %-12s  %11.2f  %11.2f  %11.2f  %11.3f  %11.3f", "whole circle", "", g, g1, g2, a1, a2);
  }

  // 2. Inverse tuning
  for (size_t k = 0; k < copies.size(); k++) {
    const int r = copies[k].w.lo.q;
    print("\n2. %s: decoded octave j vs full octave j' (TV %%, at 2^D angle resolution)", copies[k].name);
    print("    j   D  captured   j'=j   j'=j/%d   best j'  TV(best)   TV curve over j' = D..j", r);
    for (int j = 2*r; j <= max_j; j++) {
      const int D = Ds[k][j], K = 1 << (D - 1);
      string curve;
      int best = -1;
      double bt = 2, at_j = 0, at_r = 0;
      for (int jp = max(1, D - 1); jp <= j; jp++) {
        if (int(coarse[jp].size()) < K) continue;
        const auto Q = coarsen(coarse[jp], K);
        const double d = tv(decoded[k][j], Q);
        if (d < bt) { bt = d; best = jp; }
        if (jp == j) at_j = d;
        if (jp == (j + r/2) / r) at_r = d;
        curve += tfm::format(" %2.0f", 100*d);
      }
      print("   %2d  %2d  %6.1f%%  %5.0f%%  %6.0f%%  %6d  %7.0f%%  %s", j, D, 100*captured[k][j], 100*at_j, 100*at_r, best,
            100*bt, curve);
    }
    // Multiplicative vs additive tracking of the copy's captured energy
    vector<double> Ecopy(max_j + 1);
    for (int j = 1; j <= max_j; j++) {
      double s = 0;
      for (const double v : decoded[k][j]) s += v;
      Ecopy[j] = s;
    }
    const auto logS = [&](const double x) {  // log S at fractional octave, by linear interpolation of log S
      const int a = min(max(1, int(std::floor(x))), max_j - 1);
      const double f = x - a;
      return (1 - f) * std::log(S[a]) + f * std::log(S[a + 1]);
    };
    const auto spread = [&](const function<double(int)>& jp) {
      double s1 = 0, s2 = 0; int m = 0;
      for (int j = j0; j <= j1; j++) {
        const double v = std::log2(Ecopy[j]) - logS(jp(j)) / std::log(2.);
        s1 += v; s2 += v*v; m++;
      }
      return std::sqrt(max(0.0, s2/m - sqr(s1/m)));
    };
    string add;
    for (int s = 0; s <= 8; s += 2) add += tfm::format(" s=%d: %.2f", s, spread([s](int j) { return double(j - s); }));
    print("   std of log2(E_copy(j) / S(j')), j = %d..%d:  multiplicative j' = j/%d: %.2f   additive j' = j-s:%s",
          j0, j1, r, spread([r](int j) { return double(j) / r; }), add);
  }
  print("\ntotal %.1f s", (wall_time() - t0).seconds());
}

}  // namespace
}  // namespace mandelbrot

int main(const int argc, const char** argv) {
  using namespace mandelbrot;
  try {
    slow_assert(argc >= 3, "usage: %s <f-k27.npy> <max_j> [numpy angle map to compare]", argv[0]);
    run(argv[1], atoi(argv[2]), argc > 3 ? argv[3] : "");
    return 0;
  } catch (const std::exception& e) {
    die(e.what());
  }
}
