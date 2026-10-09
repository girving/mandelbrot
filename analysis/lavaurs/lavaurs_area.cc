// Batched Lavaurs-model component areas (lavaurs.h) on CPU threads.
//
//   ./build/release/lavaurs_area [threads] < jobs > out
// reads lines "name r n sigma_re sigma_im" (r transits, excursion n, a guess for the center in the phase σ) and
// prints "name r n center_re center_im area_hi area_lo C conv" with C = (π²/4) area the family constant, or
// "name r n failed".
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

int main(int argc, char** argv) {
  if (argc > 1 && std::string(argv[1]) == "--walk") return walk_main(argc > 2 ? atoi(argv[2]) : 2);
  const int threads = argc > 1 ? atoi(argv[1]) : 2;
  struct Job { std::string name; int r, n; double sr, si; };
  std::vector<Job> jobs;
  char name[256];
  int r, n;
  double sr, si;
  while (scanf("%255s %d %d %lf %lf", name, &r, &n, &sr, &si) == 5) jobs.push_back({name, r, n, sr, si});
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
          snprintf(line, sizeof(line), "%s %d %d %.17g %.17g %.17g %.17g %.17g %.1e", j.name.c_str(), j.r, j.n,
                   res.center.r, res.center.i, res.area.x[0], res.area.x[1], M_PI * M_PI / 4 * double(res.area),
                   res.conv);
        else {
          snprintf(line, sizeof(line), "%s %d %d failed", j.name.c_str(), j.r, j.n);
          failed++;
        }
        out[i] = line;
      }
    });
  for (auto& t : pool) t.join();
  for (const auto& s : out) printf("%s\n", s.c_str());
  fprintf(stderr, "lavaurs_area: %zu jobs, %lld failed\n", jobs.size(), (long long)failed);
}
