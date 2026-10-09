// Batched Lavaurs-model component areas (lavaurs.h) on CPU threads.
//
//   ./build/release/lavaurs_area [threads] < jobs > out
// reads lines "name r n sigma_re sigma_im" (r transits, excursion n, a guess for the center in the phase σ) and
// prints "name r n center_re center_im area_hi area_lo C conv" with C = (π²/4) area the family constant, or
// "name r n failed".  With $LAVAURS_TUNE = "X:p:c_re:c_im,..." (primitive centers c_X of period p in M), each component
// U is followed by its tunings U*X: r p transits and excursion p n + p - 1 (the cycle (F^n g^r F)^p), from the guess
// center + 2 a_1 c_X (the copy's cardioid is μ/2 - μ²/4), printed as "name*X ..." with the area ratio to U appended.
//
//   ./build/release/lavaurs_area --island [threads] < pairs > out
// reads lines "name n_u u_re u_im n_c c_re c_im" (an island: single-transit center σ_u with excursion n_u; a target:
// single-transit center σ_c with excursion n_c) and finds the two-transit centers Θ(σ) = σ_c + j/2 (j = jmin..1 with
// jmin = $LAVAURS_JMIN or -1; with LAVAURS_R = r > 2 the sources are (r-1)-transit centers and the results r-transit),
// final excursion n_c - j) on both branches, printing "name|side|j 2 n center_re center_im area_hi area_lo C conv
// island_address target_address side own|other island_re island_im" (the combinatorial label: addresses of the island's and the target's preimages
// of the critical point, the branch's side at the critical passage, and j).  With --address, ordinary output lines
// end with the single-transit address of the center.
//
//   ./build/release/lavaurs_area --labels < centers > out
// reads lines "name n sigma_re sigma_im" (single-transit centers) and prints "name n address" (no areas).
//
//   ./build/release/lavaurs_area --locate < "name r sigma_re sigma_im" lines
// prints "name r theta_re theta_im center_re center_im": Θ_r(σ) (an r-transit center maps to a single-transit center,
// shifted by j/2) and the critical point of the horn-map composition reached by Newton (the region's source center).
//
//   ./build/release/lavaurs_area --walk [threads] < walks > out
// reads lines "name r n0 step count s1_re s1_im s2_re s2_im ..." (a class of components whose excursion grows by step:
// seeds for its first members, from M) and walks n = n0, n0 + step, ...: members past the seeds are predicted by
// extrapolating the previous centers in n^{-1/2} (Lagrange through up to 7), and a member is accepted only if its center
// lies within 0.3 of its radius sqrt(area/π) from the prediction (else the walk stops).  Output lines as above, with
// name_<n>.
#include "expansion_arith.h"
#include "lavaurs.h"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <thread>
#include <vector>
using namespace mandelbrot;

static Complex<double> cdiv(const Complex<double> a, const Complex<double> b) {
  const double d = b.r * b.r + b.i * b.i;
  return Complex<double>((a.r * b.r + a.i * b.i) / d, (a.i * b.r - a.r * b.i) / d);
}

static int walk_main(const int threads) {
  struct Walk { std::string name; int r, n0, step, count; std::vector<Complex<double>> seeds; };
  std::vector<Walk> walks;
  char buf[1 << 16];
  while (fgets(buf, sizeof(buf), stdin)) {
    Walk w;
    char name[256];
    int used = 0;
    if (sscanf(buf, "%255s %d %d %d %d%n", name, &w.r, &w.n0, &w.step, &w.count, &used) != 5) continue;
    w.name = name;
    const char* p = buf + used;
    double a, b;
    int k;
    while (sscanf(p, "%lf %lf%n", &a, &b, &k) == 2) { w.seeds.push_back(Complex<double>(a, b)); p += k; }
    walks.push_back(w);
  }
  std::vector<std::string> out(walks.size());
  std::atomic<int64_t> next(0), stopped(0);
  std::vector<std::thread> pool;
  for (int t = 0; t < threads; t++)
    pool.emplace_back([&]() {
      for (int64_t i; (i = next++) < int64_t(walks.size());) {
        const auto& w = walks[i];
        std::vector<Complex<double>> cs;
        std::string text;
        for (int m = 0; m < w.count; m++) {
          const int n = w.n0 + m * w.step;
          Complex<double> pred;
          const int K = std::min(int(cs.size()), 7);
          if (m < int(w.seeds.size())) pred = w.seeds[m];
          else if (K >= 2) {
            // Lagrange through the last K centers in x = n^{-1/2} (a long digit is a pass of ~n/2 turns landing ~n^{-1/2}
            // from the parabolic point, so the centers are smooth in x), evaluated at this n
            const auto x = [&](const int mm) { return 1 / std::sqrt(double(w.n0 + mm * w.step)); };
            pred = Complex<double>(0, 0);
            for (int a = 0; a < K; a++) {
              double l = 1;
              for (int b = 0; b < K; b++)
                if (b != a) l *= (x(m) - x(m - K + b)) / (x(m - K + a) - x(m - K + b));
              pred = pred + Complex<double>(l * cs[cs.size() - K + a].r, l * cs[cs.size() - K + a].i);
            }
          } else break;
          const auto res = lavaurs_area(w.r, n, pred);
          if (!res.ok) { fprintf(stderr, "walk %s: n = %d failed (pred %.12f%+.12fi, %zu seeds)\n", w.name.c_str(), n, pred.r, pred.i, w.seeds.size()); stopped++; break; }
          const double rad = std::sqrt(double(res.area) / M_PI);
          const double move = std::hypot(res.center.r - pred.r, res.center.i - pred.i);
          if (m >= int(w.seeds.size()) && move > 0.3 * rad) { fprintf(stderr, "walk %s: n = %d moved %.2g radii\n", w.name.c_str(), n, move / rad); stopped++; break; }
          cs.push_back(res.center);
          char line[512];
          snprintf(line, sizeof(line), "%s_%d %d %d %.17g %.17g %.17g %.17g %.17g %.1e\n", w.name.c_str(), n, w.r, n,
                   res.center.r, res.center.i, res.area.x[0], res.area.x[1], M_PI * M_PI / 4 * double(res.area),
                   res.conv);
          text += line;
        }
        out[i] = text;
      }
    });
  for (auto& t : pool) t.join();
  for (const auto& s : out) fputs(s.c_str(), stdout);
  fprintf(stderr, "lavaurs_area --walk: %zu walks, %lld stopped early\n", walks.size(), (long long)stopped);
  return 0;
}

static int jmin = -1;  // Shifts j = jmin..1 (env LAVAURS_JMIN)
static int transits = 2;  // r (env LAVAURS_R): sources are (r-1)-transit centers

static int island_main(const int threads) {
  if (getenv("LAVAURS_JMIN")) jmin = atoi(getenv("LAVAURS_JMIN"));
  if (getenv("LAVAURS_R")) transits = atoi(getenv("LAVAURS_R"));
  struct Pair { std::string name; int nu, nc; Complex<double> u, c; };
  std::vector<Pair> pairs;
  char name[256];
  int nu, nc;
  double ur, ui, cr, ci;
  while (scanf("%255s %d %lf %lf %d %lf %lf", name, &nu, &ur, &ui, &nc, &cr, &ci) == 7)
    pairs.push_back({name, nu, nc, Complex<double>(ur, ui), Complex<double>(cr, ci)});
  std::vector<std::string> out(pairs.size());
  std::atomic<int64_t> next(0), found(0);
  std::vector<std::thread> pool;
  for (int t = 0; t < threads; t++)
    pool.emplace_back([&]() {
      for (int64_t i; (i = next++) < int64_t(pairs.size());) {
        const auto& p = pairs[i];
        const std::string ua = lavaurs_address(p.nu, p.u);
        std::string text;
        for (const int branch : {1, -1})
          for (int j = jmin; j <= 1; j++) {
            const int n = p.nc - j;
            if (n < 0) continue;
            Complex<double> s;
            if (!lavaurs_island(p.u, Complex<double>(p.c.r + 0.5 * j, p.c.i), branch, s, transits)) continue;
            const auto res = lavaurs_area(transits, n, s);
            if (!res.ok || std::hypot(res.center.r - s.r, res.center.i - s.i) > 1e-7) continue;
            const char side = lavaurs_island_side(p.nu, res.center);
            // The component's true island (its horn-map critical point), against the island it was tracked from
            Complex<double> isl = res.center;
            const bool own = lavaurs_island_center(isl, transits) && std::hypot(isl.r - p.u.r, isl.i - p.u.i) < 1e-8;
            char line[1024];
            snprintf(line, sizeof(line), "%s|%c|%d %d %d %.17g %.17g %.17g %.17g %.17g %.1e %s %s %c %s %.15g %.15g %.3e\n",
                     p.name.c_str(), side, j, transits, n, res.center.r, res.center.i, res.area.x[0], res.area.x[1],
                     M_PI * M_PI / 4 * double(res.area), res.conv, ua.empty() ? "-" : ua.c_str(),
                     p.nc ? lavaurs_address(p.nc, p.c).c_str() : "-", side, own ? "own" : "other", isl.r, isl.i, res.cusp);
            text += line;
            found++;
          }
        out[i] = text;
      }
    });
  for (auto& t : pool) t.join();
  for (const auto& s : out) fputs(s.c_str(), stdout);
  fprintf(stderr, "lavaurs_area --island: %zu pairs, %lld components\n", pairs.size(), (long long)found);
  return 0;
}

int main(int argc, char** argv) {
  if (argc > 1 && std::string(argv[1]) == "--walk") return walk_main(argc > 2 ? atoi(argv[2]) : 2);
  if (argc > 1 && std::string(argv[1]) == "--island") return island_main(argc > 2 ? atoi(argv[2]) : 2);
  if (argc > 1 && std::string(argv[1]) == "--locate") {  // "name r sigma_re sigma_im" -> Θ_r(σ) and the region's center
    char name[256];
    int r;
    double sr, si;
    while (scanf("%255s %d %lf %lf", name, &r, &sr, &si) == 4) {
      Complex<double> s(sr, si), t, d, dd, c = s;
      if (!lavaurs_theta(s, t, d, dd, r)) { printf("%s %d failed\n", name, r); continue; }
      if (lavaurs_island_center(c, r)) printf("%s %d %.15g %.15g %.15g %.15g\n", name, r, t.r, t.i, c.r, c.i);
      else printf("%s %d %.15g %.15g - -\n", name, r, t.r, t.i);
    }
    return 0;
  }
  if (argc > 1 && std::string(argv[1]) == "--labels") {
    char name[256];
    int n;
    double sr, si;
    while (scanf("%255s %d %lf %lf", name, &n, &sr, &si) == 4)
      printf("%s %d %s\n", name, n, n ? lavaurs_address(n, Complex<double>(sr, si)).c_str() : "-");
    return 0;
  }
  const bool address = argc > 1 && std::string(argv[1]) == "--address";
  if (address) { argc--; argv++; }
  const int threads = argc > 1 ? atoi(argv[1]) : 2;
  struct Job { std::string name; int r, n; double sr, si; };
  std::vector<Job> jobs;
  char name[256];
  int r, n;
  double sr, si;
  while (scanf("%255s %d %d %lf %lf", name, &r, &n, &sr, &si) == 5) jobs.push_back({name, r, n, sr, si});
  struct Tune { std::string name; int p; Complex<double> c; };
  std::vector<Tune> tunes;
  if (const char* e = getenv("LAVAURS_TUNE")) {
    std::string t(e);
    for (size_t i = 0; i < t.size();) {
      size_t j = t.find(',', i);
      if (j == std::string::npos) j = t.size();
      char tn[64];
      int p;
      double cr, ci;
      if (sscanf(t.substr(i, j - i).c_str(), "%63[^:]:%d:%lf:%lf", tn, &p, &cr, &ci) == 4) tunes.push_back({tn, p, {cr, ci}});
      i = j + 1;
    }
  }
  std::vector<std::string> out(jobs.size());
  std::atomic<int64_t> next(0), failed(0);
  std::vector<std::thread> pool;
  for (int t = 0; t < threads; t++)
    pool.emplace_back([&]() {
      for (int64_t i; (i = next++) < int64_t(jobs.size());) {
        const auto& j = jobs[i];
        const auto res = lavaurs_area(j.r, j.n, Complex<double>(j.sr, j.si));
        char line[512];
        if (res.ok)
          snprintf(line, sizeof(line), "%s %d %d %.17g %.17g %.17g %.17g %.17g %.1e %.3e%s%s", j.name.c_str(), j.r, j.n,
                   res.center.r, res.center.i, res.area.x[0], res.area.x[1], M_PI * M_PI / 4 * double(res.area),
                   res.conv, res.cusp, address ? " " : "", address ? (j.n ? lavaurs_address(j.n, res.center).c_str() : "-") : "");
        else {
          snprintf(line, sizeof(line), "%s %d %d failed", j.name.c_str(), j.r, j.n);
          failed++;
        }
        out[i] = line;
        if (!res.ok) continue;
        // The guess interpolates the copy map c -> σ: σ(0) = center, σ'(0) = 2 a_1, and the tunings found so far
        // (σ = center + 2 a_1 c + c² q(c), q Lagrange through them), so list $LAVAURS_TUNE in increasing |c|
        std::vector<std::pair<Complex<double>, Complex<double>>> known;  // (c_X, q value)
        for (const auto& x : tunes) {
          const int r2 = j.r * x.p, n2 = x.p * j.n + x.p - 1;
          Complex<double> q(0);
          for (size_t a = 0; a < known.size(); a++) {
            Complex<double> l = known[a].second;
            for (size_t b = 0; b < known.size(); b++)
              if (b != a) l = cdiv(l * (x.c - known[b].first), known[a].first - known[b].first);
            q = q + l;
          }
          const auto g = res.center + Complex<double>(2, 0) * res.a1 * x.c + x.c * x.c * q;
          auto t = lavaurs_area(r2, n2, g);
          if (!t.ok && known.size()) t = lavaurs_area(r2, n2, res.center + Complex<double>(2, 0) * res.a1 * x.c);
          if (t.ok) known.push_back({x.c, cdiv(t.center - res.center - Complex<double>(2, 0) * res.a1 * x.c, x.c * x.c)});
          if (t.ok)
            snprintf(line, sizeof(line), "\n%s*%s %d %d %.17g %.17g %.17g %.17g %.17g %.1e %.3e %.6e", j.name.c_str(),
                     x.name.c_str(), r2, n2, t.center.r, t.center.i, t.area.x[0], t.area.x[1],
                     M_PI * M_PI / 4 * double(t.area), t.conv, t.cusp, double(t.area) / double(res.area));
          else
            snprintf(line, sizeof(line), "\n%s*%s %d %d failed", j.name.c_str(), x.name.c_str(), r2, n2);
          out[i] += line;
        }
      }
    });
  for (auto& t : pool) t.join();
  for (const auto& s : out) printf("%s\n", s.c_str());
  fprintf(stderr, "lavaurs_area: %zu jobs, %lld failed\n", jobs.size(), (long long)failed);
}
