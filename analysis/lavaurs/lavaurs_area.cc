// Batched Lavaurs-model component areas (lavaurs.h) on CPU threads.
//
//   ./build/release/lavaurs_area [threads] < jobs > out
// reads lines "name r n sigma_re sigma_im" (r transits, excursion n, a guess for the center in the phase σ) and
// prints "name r n center_re center_im area_hi area_lo C conv" with C = (π²/4) area the family constant, or
// "name r n failed".
#include "expansion_arith.h"
#include "lavaurs.h"
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <thread>
#include <vector>
using namespace mandelbrot;

int main(int argc, char** argv) {
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
