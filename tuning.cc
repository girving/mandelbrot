// Renormalization structure of hyperbolic component areas, and the tail's Dirichlet series D(s)
//
// Near a tuned copy of period p, g_M(c) ≈ C g_M(χ(c))^p, so the copy's share of the escape-time tail is
// r_W T((k - c_W) / p) with r_W close to its area ratio (notes §5.3).  Over maximal copies (tunings by
// non-renormalizable components W0), T(k) = T_0(k) + Σ r_W T((k - c_W) / p_W), whose Mellin transform gives
// T̂(s) (1 - D(s)) = T̂_0(s) with D(s) = Σ_W r_W p_W^s.  The tail's form at large k follows from the roots of
// D(s) = 1.
//
// Components come from `hyperbolic max_p N 1 dump.txt` (center and area per component); their roots come from
// Lavaurs' algorithm, and tuning from the angles (maximal_tuning).  Each root is matched to its center by
// tracing its lower parameter ray inward to near the root and taking the nearest center of that period; the
// matching must be a bijection.  Weights: r_W = area(W) / area(cardioid).

#include "angles.h"
#include "debug.h"
#include "print.h"
#include "wall_time.h"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <complex>
#include <thread>
#include <vector>
namespace mandelbrot {
namespace {

using std::max;
using std::min;
using std::vector;
typedef std::complex<double> C;

struct Component { int p; C c; double area; };

vector<Component> read_components(const string& path) {
  FILE* f = fopen(path.c_str(), "r");
  slow_assert(f, "can't open %s", path);
  vector<Component> cs;
  int p;
  double x, y, a;
  while (fscanf(f, "%d %lf %lf %lf", &p, &x, &y, &a) == 4) cs.push_back({p, C(x, y), a});
  fclose(f);
  return cs;
}

// Trace the parameter ray of angle θ = k / (2^q - 1) inward through `levels` levels of potential, following
// Φ(c) = ρ e^{2πiθ} via z_n(c) = ρ^{2^(n-1)} e^{2πi 2^(n-1) θ} with S Newton-tracked points per level.  The
// endpoint has potential about log(er) 2^-levels.
C ray_in(const Periodic& theta, const int levels, const int S = 8) {
  const double er = 65536;
  const uint64_t M = (uint64_t(1) << theta.q) - 1;
  uint64_t k = theta.k;  // 2^(n-1) θ = k / M (mod 1)
  C c = std::polar(er, 2 * M_PI * double(k) / double(M));
  for (int n = 1; n <= levels; n++) {
    for (int j = 0; j < S; j++) {
      const double r = std::pow(er, std::pow(0.5, (j + 1.0) / S));
      const C t = std::polar(r, 2 * M_PI * double(k) / double(M));
      for (int it = 0; it < 64; it++) {
        C z = 0, dz = 0;
        for (int m = 0; m < n; m++) { dz = 2.0 * z * dz + 1.0; z = z * z + c; }
        const C step = (z - t) / dz;
        c -= step;
        if (std::abs(step) < 1e-15 * std::abs(c)) break;
      }
    }
    k = 2 * k % M;  // Next level targets z_{n+1} = r² e^{2πi 2θ_n}, where r² = er at the level's start
  }
  return c;
}

void run(const string& path, const int P, const int extra_levels, const double wprim, const double wsat) {
  auto t0 = wall_time();
  const auto comps = read_components(path);
  vector<vector<int>> by_p(P + 1);
  for (size_t i = 0; i < comps.size(); i++) if (comps[i].p <= P) by_p[comps[i].p].push_back(int(i));
  slow_assert(by_p[1].size() == 1, "expected one period-1 component");
  const double a1 = comps[by_p[1][0]].area;
  const auto roots = lavaurs(P);
  const auto parent = maximal_tuning(roots);
  print("%d components, %d roots of period 2..%d: %.2f s", comps.size(), roots.size(), P, (wall_time() - t0).seconds());

  // Match roots to centers: nearest center of the same period to the ray's endpoint, by real-part sweep
  for (int p = 1; p <= P; p++)
    std::sort(by_p[p].begin(), by_p[p].end(), [&](int a, int b) { return comps[a].c.real() < comps[b].c.real(); });
  t0 = wall_time();
  vector<int> match(roots.size(), -1);
  vector<double> margin(roots.size());  // Second-nearest / nearest distance
  std::atomic<int64_t> next(0);
  vector<std::thread> pool;
  for (unsigned t = 0; t < std::thread::hardware_concurrency(); t++)
    pool.emplace_back([&]() {
      for (int64_t i; (i = next.fetch_add(1)) < int64_t(roots.size());) {
        const int p = roots[i].w.lo.q;
        const C c = ray_in(roots[i].w.lo, 2 * p + extra_levels);
        const auto& v = by_p[p];
        const auto lo = std::lower_bound(v.begin(), v.end(), c.real(),
                                         [&](int a, double x) { return comps[a].c.real() < x; });
        double d1 = INFINITY, d2 = INFINITY;
        int best = -1;
        const auto visit = [&](const int j) {
          const double d = std::abs(comps[j].c - c);
          if (d < d1) { d2 = d1; d1 = d; best = j; } else if (d < d2) d2 = d;
        };
        for (auto it = lo; it != v.end() && comps[*it].c.real() - c.real() <= d2; ++it) visit(*it);
        for (auto it = lo; it != v.begin();) { --it; if (c.real() - comps[*it].c.real() > d2) break; visit(*it); }
        match[i] = best;
        margin[i] = d2 / d1;
      }
    });
  for (auto& t : pool) t.join();
  vector<int> hits(comps.size());
  double worst = INFINITY;
  for (size_t i = 0; i < roots.size(); i++) { hits[match[i]]++; worst = min(worst, margin[i]); }
  int64_t dup = 0, missing = 0;
  for (int p = 2; p <= P; p++)
    for (const int j : by_p[p]) { dup += hits[j] > 1; missing += hits[j] == 0; }
  print("ray matching: %.2f s, %d centers hit twice, %d missed, worst second/first distance ratio %.2f",
        (wall_time() - t0).seconds(), dup, missing, worst);
  slow_assert(!dup && !missing, "root-center matching is not a bijection");

  // Shortcut check: classify each center by its closest returns (the n where |z_n| reaches a new minimum,
  // approximating the internal address): renormalizable with period d iff d is a closest return, 1 < d < p,
  // d | p, and every later closest return is a multiple of d
  {
    int64_t agree = 0, false_tuned = 0, false_nr = 0;
    vector<int64_t> bad_p(P + 1);
    for (size_t i = 0; i < roots.size(); i++) {
      const int p = roots[i].w.lo.q;
      const C c = comps[match[i]].c;
      vector<int> returns;
      C z = 0;
      double best = INFINITY;
      for (int n = 1; n <= p; n++) {
        z = z * z + c;
        const double a = std::abs(z);
        if (a < best || n == p) { best = a; returns.push_back(n); }
      }
      bool tuned = false;
      for (size_t a = 0; a < returns.size() && !tuned; a++) {
        const int d = returns[a];
        if (d <= 1 || d >= p || p % d) continue;
        bool all = true;
        for (size_t b = a + 1; b < returns.size(); b++) all &= returns[b] % d == 0;
        tuned = all;
      }
      const bool exact = parent[i] >= 0;
      if (tuned == exact) agree++;
      else { (tuned ? false_tuned : false_nr)++; bad_p[p]++; }
    }
    string bad;
    for (int p = 2; p <= P; p++) if (bad_p[p]) bad += tfm::format(" p%d:%d", p, bad_p[p]);
    print("closest-return shortcut vs angles: %d of %d agree; %d falsely tuned, %d falsely non-renormalizable%s",
          agree, roots.size(), false_tuned, false_nr, bad.size() ? " (" + bad + " )" : "");
  }

  // Per period: counts and areas, all and non-renormalizable (satellite / primitive)
  print("\n  p        N      N_nr   (sat, prim)        a_p          a_nr_p        a_nr_p / a_p   p^3 a_nr_p");
  vector<double> anr(P + 1), anr_sat(P + 1);
  for (int p = 2; p <= P; p++) {
    int64_t N = 0, nr = 0, nr_sat = 0;
    double a = 0;
    for (size_t i = 0; i < roots.size(); i++) {
      if (roots[i].w.lo.q != p) continue;
      const double ar = comps[match[i]].area;
      N++; a += ar;
      if (parent[i] < 0) {
        nr++; anr[p] += ar;
        if (roots[i].satellite) { nr_sat++; anr_sat[p] += ar; }
      }
    }
    print("  %2d %8d %8d   (%d, %d) %*s %.6e   %.6e   %.4f         %.4f", p, N, nr, nr_sat, nr - nr_sat,
          max(0, 12 - int(tfm::format("%d, %d", nr_sat, nr - nr_sat).size())), "", a, anr[p], anr[p] / a,
          p * p * p * anr[p]);
  }

  // D(s) = Σ_{p ≥ 2} (r_p / a1) p^s with r_p the weighted non-renormalizable area, truncated at period P' ≤ P
  vector<double> r(P + 1);
  for (int p = 2; p <= P; p++) r[p] = wprim * (anr[p] - anr_sat[p]) + wsat * anr_sat[p];
  const auto D = [&](const std::complex<double> s, const int Pt) {
    std::complex<double> d = 0;
    for (int p = 2; p <= Pt; p++) d += r[p] / a1 * std::exp(s * std::log(double(p)));
    return d;
  };
  print("\n  satellite share of non-renormalizable area by period:");
  string sh = "   ";
  for (int p = 2; p <= P; p++) sh += tfm::format(" %d:%.3f", p, anr_sat[p] / anr[p]);
  print(sh);
  print("\n  D(s) truncated at period P' (weights r_W = w area / area(cardioid), w = %g primitive, %g satellite):",
        wprim, wsat);
  print("    P'    D(0)       D(1)       D(2)       real root s* of D(s) = 1");
  for (int Pt = 4; Pt <= P; Pt++) {
    double lo = -5, hi = 20;
    for (int it = 0; it < 200; it++) { const double m = (lo + hi) / 2; (D(m, Pt).real() < 1 ? lo : hi) = m; }
    print("    %2d   %.6f   %.6f   %.6f   %.4f", Pt, D(0.0, Pt).real(), D(1.0, Pt).real(), D(2.0, Pt).real(), lo);
  }

  // Complex roots of D(s) = 1 at the full truncation, by Newton from a grid of seeds
  print("\n  complex roots of D(s) = 1 (P' = %d), largest real parts first:", P);
  vector<std::complex<double>> found;
  for (double sr = -2; sr <= 4; sr += 0.5)
    for (double si = 0.5; si <= 40; si += 0.5) {
      std::complex<double> s(sr, si);
      bool ok = false;
      for (int it = 0; it < 100; it++) {
        std::complex<double> f = -1, df = 0;
        for (int p = 2; p <= P; p++) {
          const auto term = r[p] / a1 * std::exp(s * std::log(double(p)));
          f += term; df += term * std::log(double(p));
        }
        const auto step = f / df;
        s -= step;
        if (std::abs(step) < 1e-13) { ok = true; break; }
        if (std::abs(s) > 1e3) break;
      }
      if (!ok || s.imag() <= 1e-6) continue;
      bool dup = false;
      for (const auto& r : found) dup |= std::abs(r - s) < 1e-8;
      if (!dup) found.push_back(s);
    }
  std::sort(found.begin(), found.end(), [](auto a, auto b) { return a.real() > b.real(); });
  for (size_t i = 0; i < min(found.size(), size_t(12)); i++)
    print("    s = %.4f ± %.4fi  (log-period in ln k: %.3f)", found[i].real(), found[i].imag(), 2 * M_PI / found[i].imag());
}

}  // namespace
}  // namespace mandelbrot

int main(const int argc, const char** argv) {
  using namespace mandelbrot;
  try {
    slow_assert(argc >= 2, "usage: %s components.txt [max_period ≤ 16] [extra ray levels] [primitive weight] "
                "[satellite weight]", argv[0]);
    const int P = argc > 2 ? atoi(argv[2]) : 16, extra = argc > 3 ? atoi(argv[3]) : 20;
    // Tail weight per unit area ratio, measured on real copies with --box (notes §5.3): about 0.89 for
    // primitive copies of periods 3–5 and 0.63 for the period-2 satellite copy
    const double wprim = argc > 4 ? atof(argv[4]) : 1, wsat = argc > 5 ? atof(argv[5]) : 1;
    slow_assert(2 <= P && P <= 16, "max_period must be in [2, 16]");
    run(argv[1], P, extra, wprim, wsat);
    return 0;
  } catch (const std::exception& e) {
    die(e.what());
  }
}
