// Escape-octave profiles of parameters near a root, in its parent component's multiplier coordinate.
//
// A component W of period P has multiplier λ_W(c) = (f_c^P)'(z) on its attracting cycle; its p/q root is
// λ_W = e^{2πip/q}.  Sample λ in shells around the root, |λ - λ0| ∈ (ρ0 2^{-b-1}, ρ0 2^{-b}], map to c(λ)
// (explicit for the cardioid and the period 2 disk, else by Newton on (f^P(z) = z, (f^P)'(z) = λ) from the root),
// skip λ inside W (|λ| < 1), and histogram the escape octave of the critical orbit, weighted by λ-area.  If roots
// of one type look alike in multiplier coordinates, profiles at different occurrences agree after shifting time by
// log2 P (one multiplier step is P iterations), in λ-area units.
//   ./build/release/multiplier_profile P p q samples_log2 max_iter_log2 [c_re c_im]
// with (c_re, c_im) a point inside W (its center, say) for P > 2.  CPU threads: $MANDELBROT_THREADS.
#include <algorithm>
#include <atomic>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <thread>
#include <vector>
typedef std::complex<double> C;

// Newton for (z, c) with f_c^P(z) = z, (f_c^P)'(z) = lam, from (z, c)
static bool solve(const int P, const C lam, C& z, C& c) {
  for (int it = 0; it < 50; it++) {
    C x = z, xz = 1, xc = 0, xzz = 0, xzc = 0;
    for (int i = 0; i < P; i++) {
      const C nxzz = 2.0 * (xz * xz + x * xzz), nxzc = 2.0 * (xc * xz + x * xzc);
      xzz = nxzz; xzc = nxzc;
      xc = 2.0 * x * xc + 1.0;
      xz = 2.0 * x * xz;
      x = x * x + c;
    }
    const C F1 = x - z, F2 = xz - lam, a = xz - 1.0, b = xc, d = xzz, e = xzc, det = a * e - b * d;
    const C dz = (F1 * e - b * F2) / det, dc = (a * F2 - d * F1) / det;
    z -= dz; c -= dc;
    if (std::abs(dz) + std::abs(dc) < 1e-14) return true;
  }
  return false;
}

// Whether Newton on f^Q(w) = w from z finds a cycle with |multiplier| < 1 (Q = q P, the child component's period):
// certifies c interior long before a near-parabolic orbit converges
static int Q = 1;
static bool attracting(const C c, C w) {
  for (int it = 0; it < 30; it++) {
    C x = w, d = 1;
    for (int i = 0; i < Q; i++) { d = 2.0 * x * d; x = x * x + c; }
    const C step = (x - w) / (d - 1.0);
    w -= step;
    if (!(std::norm(step) < 1e10)) return false;
    if (std::norm(step) < 1e-28 * (1 + std::norm(w))) {
      C y = w, m = 1;
      for (int i = 0; i < Q; i++) { m = 2.0 * y * m; y = y * y + c; }
      return std::norm(m) < 1 - 1e-9 && std::norm(y - w) < 1e-20;
    }
  }
  return false;
}

int main(int argc, char** argv) {
  if (argc < 6) { fprintf(stderr, "usage: multiplier_profile P p q samples_log2 max_iter_log2 [c_re c_im]\n"); return 1; }
  const int P = atoi(argv[1]), pp = atoi(argv[2]), q = atoi(argv[3]);
  const int64_t S = int64_t(1) << atoi(argv[4]), max_iter = int64_t(1) << atoi(argv[5]);
  const C lam0 = std::polar(1.0, 2 * M_PI * pp / q);
  Q = q * P;
  const double rho0 = 0.125;
  const int bins = 12;  // Down to ρ0 2^-12: escape times there (~2πP/(qρ)) stay well below max_iter
  // The root: explicit for P = 1, 2; else continue from the center (λ = 0, z = 0) to λ0
  C zr = 0, cr = 0;
  if (P == 1) cr = lam0 / 2.0 - lam0 * lam0 / 4.0;
  else if (P == 2) cr = lam0 / 4.0 - 1.0;
  else {
    cr = C(atof(argv[6]), atof(argv[7]));
    zr = 0;
    for (int s = 1; s <= 64; s++)
      if (!solve(P, lam0 * (s / 64.0), zr, cr)) { fprintf(stderr, "continuation failed\n"); return 1; }
  }
  const auto c_of = [&](const C lam, C& c) {
    if (P == 1) { c = lam / 2.0 - lam * lam / 4.0; return true; }
    if (P == 2) { c = lam / 4.0 - 1.0; return true; }
    C z = zr; c = cr;
    return solve(P, lam, z, c);
  };
  const char* te = getenv("MANDELBROT_THREADS");
  const int T = te ? atoi(te) : int(std::thread::hardware_concurrency());
  printf("P %d, root %d/%d: c0 = %.15f %+.15fi; λ shells %g … %g, %lld samples, max_iter 2^%s\n", P, pp, q,
         cr.real(), cr.imag(), rho0, std::ldexp(rho0, -bins), (long long)S, argv[5]);
  std::vector<std::vector<int64_t>> hist(T, std::vector<int64_t>(bins * 66, 0));
  std::atomic<int64_t> fails(0);
  std::vector<std::thread> pool;
  for (int t = 0; t < T; t++)
    pool.emplace_back([&, t]() {
      std::mt19937_64 rng(77 + t);
      std::uniform_real_distribution<double> u(0, 1);
      for (int64_t s = t; s < S; s += T) {
        const int b = int(s % bins);
        const double hi = std::ldexp(rho0, -b), lo = hi / 2;
        const double r = std::sqrt(lo * lo + u(rng) * (hi * hi - lo * lo));
        const C lam = lam0 + std::polar(r, 2 * M_PI * u(rng));
        int out = 64;  // Interior
        C c;
        if (std::norm(lam) >= 1) {
          if (!c_of(lam, c)) { fails++; continue; }
          C z = 0, check = 0;
          int64_t next_check = 16;
          out = 65;
          for (int64_t i = 0; i < max_iter; i++) {
            if (std::norm(z) > 4) { out = 63 - __builtin_clzll(uint64_t(i) | 1); break; }
            z = z * z + c;
            if (std::norm(z - check) < 1e-24) { out = 64; break; }
            if (i == next_check) {
              check = z;
              next_check *= 2;
              if (attracting(c, z)) { out = 64; break; }
            }
          }
        }
        hist[t][b * 66 + out]++;
      }
    });
  for (auto& th : pool) th.join();
  std::vector<double> D(66, 0);
  for (int b = 0; b < bins; b++) {
    const double hi = std::ldexp(rho0, -b), lo = hi / 2, area = M_PI * (hi * hi - lo * lo);
    int64_t nb = 0, cnt[66] = {};
    for (int t = 0; t < T; t++) for (int j = 0; j < 66; j++) { cnt[j] += hist[t][b * 66 + j]; nb += hist[t][b * 66 + j]; }
    for (int j = 0; j < 66; j++) D[j] += nb ? area * cnt[j] / double(nb) : 0;
  }
  printf("Newton failures %lld; undecided λ-area %.3e\n", (long long)fails.load(), D[65]);
  printf("  j   j - log2 P   D_j (λ-area)   n^3-scaled D_j 2^{3(j - log2 P)}\n");
  for (int j = 0; j < 64; j++)
    if (D[j]) printf("  %2d   %6.2f      %.5e     %.5f\n", j, j - std::log2(double(P)), D[j], D[j] * std::pow(2.0, 3 * (j - std::log2(double(P)))));
}
