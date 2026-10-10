// The general p/q-root Lavaurs model (glavaurs.h): components, areas, and children.
//
//   glavaurs p q single nmax re0 re1 im0 im1 grid threads   # single-transit centers by grid Newton, with areas
//   glavaurs p q area threads < "name r n re im"              # areas: "name r n cre cim area_σ C conv cusp"
// C = (4π² sin²(πp/q)/q⁴) area_σ is the family constant lim k⁴ area_M (limbs [CF(p/q), k]).
#include "glavaurs.h"
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <string>
#include <thread>
#include <vector>
using namespace mandelbrot;
typedef Complex<double> Cd;

int main(int argc, char** argv) {
  if (argc < 4) { fprintf(stderr, "usage: glavaurs p q mode ...\n"); return 1; }
  const int p = atoi(argv[1]), q = atoi(argv[2]);
  const std::string mode = argv[3];
  const GeneralLavaurs L(p, q, getenv("GL_SIDE") ? atoi(getenv("GL_SIDE")) : 1);
  const double K = 4 * M_PI * M_PI * std::pow(std::sin(M_PI * p / q), 2) / std::pow(q, 4);
  fprintf(stderr, "glavaurs %d/%d: a_{-q} = %.12g%+.12gi, β = %.12g%+.12gi, ζ0 = %.12g%+.12gi, K = %.10g\n", p, q,
          L.a[0].r, L.a[0].i, L.beta.r, L.beta.i, L.zeta0().r, L.zeta0().i, K);
  if (mode == "single") {
    const int nmax = atoi(argv[4]), grid = atoi(argv[9]), threads = atoi(argv[10]);
    const double r0 = atof(argv[5]), r1 = atof(argv[6]), i0 = atof(argv[7]), i1 = atof(argv[8]);
    struct Job { int n; Cd g; };
    std::vector<Job> jobs;
    for (int n = 0; n <= nmax; n++)
      for (int a = 0; a < grid; a++)
        for (int b = 0; b < grid; b++)
          jobs.push_back({n, Cd(r0 + (r1 - r0) * (a + 0.5) / grid, i0 + (i1 - i0) * (b + 0.5) / grid)});
    std::mutex mu;
    std::vector<std::pair<int, Cd>> found;
    std::atomic<int64_t> next(0);
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; t++)
      pool.emplace_back([&]() {
        for (int64_t i; (i = next++) < int64_t(jobs.size());) {
          Cd s = jobs[i].g;
          if (!L.center(1, jobs[i].n, s)) continue;
          if (!(s.r >= r0 - 1 && s.r <= r1 + 1 && s.i >= i0 - 1 && s.i <= i1 + 1)) continue;
          std::lock_guard<std::mutex> lock(mu);
          bool dup = false;
          for (const auto& f : found) dup |= f.first == jobs[i].n && std::hypot(f.second.r - s.r, f.second.i - s.i) < 1e-8;
          if (!dup) found.push_back({jobs[i].n, s});
        }
      });
    for (auto& t : pool) t.join();
    for (const auto& f : found) {
      Cd cen, a1;
      double ar, conv, cusp;
      if (L.area(1, f.first, f.second, cen, ar, conv, cusp, a1))
        printf("S 1 %d %.17g %.17g %.10e %.10e %.1e %.3e\n", f.first, cen.r, cen.i, ar, K * ar, conv, cusp);
      else
        printf("S 1 %d %.17g %.17g area-failed\n", f.first, f.second.r, f.second.i);
    }
    return 0;
  }
  if (mode == "area") {
    const int threads = argc > 4 ? atoi(argv[4]) : 2;
    struct Job { std::string name; int r, n; Cd g; };
    std::vector<Job> jobs;
    char name[256];
    int r, n;
    double sr, si;
    while (scanf("%255s %d %d %lf %lf", name, &r, &n, &sr, &si) == 5) jobs.push_back({name, r, n, Cd(sr, si)});
    std::vector<std::string> out(jobs.size());
    std::atomic<int64_t> next(0);
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; t++)
      pool.emplace_back([&]() {
        for (int64_t i; (i = next++) < int64_t(jobs.size());) {
          const auto& j = jobs[i];
          Cd cen, a1;
          double ar, conv, cusp;
          char line[512];
          if (L.area(j.r, j.n, j.g, cen, ar, conv, cusp, a1))
            snprintf(line, sizeof(line), "%s %d %d %.17g %.17g %.10e %.10e %.1e %.3e", j.name.c_str(), j.r, j.n, cen.r, cen.i, ar, K * ar, conv, cusp);
          else
            snprintf(line, sizeof(line), "%s %d %d failed", j.name.c_str(), j.r, j.n);
          out[i] = line;
        }
      });
    for (auto& t : pool) t.join();
    for (const auto& s : out) printf("%s\n", s.c_str());
    return 0;
  }
  if (mode == "debug") {
    const Cd sig(atof(argv[4]), atof(argv[5]));
    Cd s0, d0, dd0;
    int k;
    const bool ok = L.phi_a(L.v, s0, d0, dd0, k);
    printf("phi_a(v) ok %d = %.12g%+.12gi petal %d (d %.6g%+.6gi)\n", ok, s0.r, s0.i, k, d0.r, d0.i);
    for (int kk = 0; kk < q; kk++) {
      Cd w, d, dd;
      const bool ok2 = L.psi(s0 + sig, kk, w, d, dd);
      if (kk == L.exit_petal(k)) printf("(exit petal) ");
      Cd fw = L.lam * w + w * w;
      printf("psi(ζ0+σ, %d) ok %d = %.12g%+.12gi; f of it %.6g%+.6gi (crit %.6g%+.6gi)\n", kk, ok2, w.r, w.i, fw.r, fw.i, L.crit.r, L.crit.i);
    }
    Cd sg = sig;
    printf("center(1,1) ok %d -> %.12g%+.12gi\n", L.center(1, 1, sg), sg.r, sg.i);
    return 0;
  }
  fprintf(stderr, "unknown mode %s\n", mode.c_str());
  return 1;
}
