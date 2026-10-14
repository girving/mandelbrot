// Parameter-space escape times near roots of the main cardioid, in the multiplier coordinate.
//
// The cardioid is c(λ) = λ/2 - λ^2/4 with λ the multiplier of the fixed point α, and its p/q root (p/q = 0/1: the
// cusp c = 1/4) is λ0 = e^{2πip/q}.  For parameters c(λ) on circles |λ - λ0| = ρ outside M, iterate the critical
// orbit until escape: if the clock near a root is the multiplier, the escape time n scales like 1/ρ (near a
// satellite root λ - λ0 ∝ c - c0, at the cusp ∝ √(c - c0)), and n ρ has a ρ-independent distribution.
//   ./build/release/root_escape p q rho_log2_min rho_log2_max angles max_iter_log2
// prints, per ρ = 2^-k: the fraction of the circle escaping, quantiles of n ρ over escaping points, and the
// fraction neither escaping nor certified interior (cardioid, or a periodic orbit found) by max_iter.
#include <algorithm>
#include <atomic>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <vector>
typedef std::complex<double> C;

int main(int argc, char** argv) {
  if (argc < 7) { fprintf(stderr, "usage: root_escape p q rho_log2_min rho_log2_max angles max_iter_log2\n"); return 1; }
  const int pp = atoi(argv[1]), q = atoi(argv[2]), kmin = atoi(argv[3]), kmax = atoi(argv[4]), A = atoi(argv[5]);
  const int64_t max_iter = int64_t(1) << atoi(argv[6]);
  const C lam0 = std::polar(1.0, 2 * M_PI * pp / q), c0 = lam0 / 2.0 - lam0 * lam0 / 4.0;
  const char* threads_env = getenv("MANDELBROT_THREADS");
  const int T = threads_env ? atoi(threads_env) : int(std::thread::hardware_concurrency());
  printf("root %d/%d: c0 = %.15f %+.15fi, λ0 = %.6f %+.6fi; %d angles per circle, max_iter 2^%s, %d threads\n", pp, q,
         c0.real(), c0.imag(), lam0.real(), lam0.imag(), A, argv[6], T);
  printf("  k (ρ = 2^-k)  escaped  n ρ quantiles: 10%%      25%%      50%%      75%%      90%%     (undecided)\n");
  for (int k = kmin; k <= kmax; k++) {
    const double rho = std::ldexp(1.0, -k);
    std::vector<int64_t> n(A);
    std::atomic<int> next(0);
    std::vector<std::thread> pool;
    for (int t = 0; t < T; t++)
      pool.emplace_back([&]() {
        for (;;) {
          const int a = next.fetch_add(1);
          if (a >= A) return;
          const C lam = lam0 + std::polar(rho, 2 * M_PI * (a + 0.5) / A), c = lam / 2.0 - lam * lam / 4.0;
          if (std::norm(lam) <= 1) { n[a] = -2; continue; }  // In the cardioid
          // Iterate, with a periodicity check against a checkpoint refreshed at powers of two (catches the bulb's
          // attracting cycle, slowly near the root)
          C z = 0, check = 0;
          int64_t i = 0, next_check = 16;
          n[a] = -1;
          for (; i < max_iter; i++) {
            if (std::norm(z) > 4) { n[a] = i; break; }
            z = z * z + c;
            if (std::norm(z - check) < 1e-24) { n[a] = -2; break; }
            if (i == next_check) { check = z; next_check *= 2; }
          }
        }
      });
    for (auto& t : pool) t.join();
    std::vector<double> e;
    int64_t unesc = 0;
    for (const auto x : n) { if (x >= 0) e.push_back(double(x) * rho); else if (x == -1) unesc++; }
    std::sort(e.begin(), e.end());
    const auto qt = [&](const double f) { return e.empty() ? 0.0 : e[std::min(e.size() - 1, size_t(f * e.size()))]; };
    printf("  %2d            %.4f   %8.4f %8.4f %8.4f %8.4f %8.4f   (%.4f)\n", k, double(e.size()) / A, qt(0.1), qt(0.25),
           qt(0.5), qt(0.75), qt(0.9), double(unesc) / A);
    fflush(stdout);
  }
}
