// Escape-time tails of filled Julia sets at and near parabolic parameters, for any parabolic cycle.
//
// c has a parabolic cycle z_0, ..., z_{P-1} whose multiplier is a primitive q-th root of unity, so near each point
// f^{Pq}(z_k + w) = z_k + w + a_k w^{q+1} + ... (computed by composing the Taylor series of f along the cycle).
// Interior points are certified in the attracting petals Re(-1/(q a_k w^q)) > R of the Fatou coordinate, where
// f^{Pq} acts as u ↦ u + 1 + ....  Two modes:
//   global: sample z uniformly in |z| < 2 and histogram escape steps by octave: T(n) ~ n^-(1 + 2/q);
//   local: sample log-uniformly in |z - z_k| ∈ [r_min, r0] around every cycle point, weighted by area, and
//          estimate the local band constant B in T_local(n) ≈ B n^-(1 + 2/q) per octave.
//   near: c inside the main cardioid (near a root), interior certified in |z - α| < (1 - |λ|)/2, which maps into
//         itself (α the attracting fixed point, λ = 2α): the crossover from n^-(1 + 2/q) to exponential decay.
// The ratio of global to local constants is the total weight W of the cycle's preimages.
//   ./build/release/parabolic_tail global|local|near c_re c_im z_re z_im P q samples_log2 max_iter_log2 threads [r0]
// with (z_re, z_im) a guess for a cycle point (unused by near).
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <thread>
#include <vector>
typedef std::complex<double> C;

// Taylor coefficients (degree ≤ D) of f^m(z_k + w) - z_{k+m} in w, composing f(z_j + w) = z_{j+1} + 2 z_j w + w^2
static std::vector<C> germ(const std::vector<C>& z, const int k, const int m, const int D) {
  std::vector<C> s(D + 1, 0.0);
  s[1] = 1;
  const int P = int(z.size());
  for (int it = 0; it < m; it++) {
    const C zj = z[(k + it) % P];
    std::vector<C> t(D + 1, 0.0);
    for (int i = 0; i <= D; i++) t[i] += 2.0 * zj * s[i];
    for (int i = 0; i <= D; i++) for (int j = 0; i + j <= D; j++) t[i + j] += s[i] * s[j];
    s = t;
  }
  return s;
}

int main(int argc, char** argv) {
  if (argc < 11) { fprintf(stderr, "usage: parabolic_tail global|local c_re c_im z_re z_im P q samples_log2 max_iter_log2 threads [r0]\n"); return 1; }
  const bool local = !strcmp(argv[1], "local"), near = !strcmp(argv[1], "near");
  const C c(atof(argv[2]), atof(argv[3]));
  C z0(atof(argv[4]), atof(argv[5]));
  const int P = atoi(argv[6]), q = atoi(argv[7]), lmi = atoi(argv[9]), T = atoi(argv[10]);
  const int64_t S = int64_t(1) << atoi(argv[8]), max_iter = int64_t(1) << lmi;
  const double r0 = argc > 11 ? atof(argv[11]) : 0.05, s = 1 + 2.0 / q, R = 50;

  const C alpha = (1.0 - std::sqrt(1.0 - 4.0 * c)) / 2.0;
  const double lam = std::abs(2.0 * alpha), rdisk = (1 - lam) / 2;
  if (near) {
    if (lam >= 1) { fprintf(stderr, "near: c is not inside the main cardioid\n"); return 1; }
    printf("near: c = %.15f %+.15fi, |λ| = %.12f, 1 - |λ| = %.3e, certifying |z - α| < %.3e\n", c.real(), c.imag(),
           lam, 1 - lam, rdisk);
  }
  // Refine the cycle point: Newton on f^P(z) - z for q > 1 (a simple fixed point of f^P), and on (f^P)'(z) - 1
  // for q = 1 (where the fixed point is double but the multiplier condition is simple)
  for (int it = 0; it < 100; it++) {
    C x = z0, d = 1, dd = 0;  // f^P, (f^P)', (f^P)''
    for (int i = 0; i < P; i++) { dd = 2.0 * (d * d + x * dd); d = 2.0 * x * d; x = x * x + c; }
    const C step = q > 1 ? (x - z0) / (d - 1.0) : (d - 1.0) / dd;
    z0 -= step;
    if (std::abs(step) < 1e-15) break;
  }
  std::vector<C> z(P);
  z[0] = z0;
  for (int i = 1; i < P; i++) z[i] = z[i - 1] * z[i - 1] + c;
  C mult = 1;
  for (const C x : z) mult *= 2.0 * x;
  printf("c = %.15f %+.15fi, P %d, q %d: cycle multiplier %.12f %+.12fi (|·| - 1 = %.1e, arg q/2π = %.6f)\n",
         c.real(), c.imag(), P, q, mult.real(), mult.imag(), std::abs(mult) - 1, std::arg(mult) * q / (2 * M_PI));
  std::vector<C> a(P), k(P);
  for (int i = 0; i < P; i++) {
    const auto g = germ(z, i, P * q, q + 2);
    for (int j = 2; j <= q; j++)
      if (std::abs(g[j]) > 1e-8 * std::pow(std::abs(g[q + 1]), double(j - 1) / q))
        printf("warning: point %d: w^%d coefficient %.3g\n", i, j, std::abs(g[j]));
    a[i] = g[q + 1];
    k[i] = -1.0 / (double(q) * a[i]);
    printf("  z_%d = %.12f %+.12fi, a = %.6f %+.6fi (|a| %.6g)\n", i, z[i].real(), z[i].imag(), a[i].real(),
           a[i].imag(), std::abs(a[i]));
  }
  double rcert = 1;  // Certification only within this distance of the cycle (the germ's scale)
  for (int i = 0; i < P; i++) rcert = std::min(rcert, 0.3 * std::pow(std::abs(a[i]), -1.0 / q));

  // Escape steps (octave index), 64 for certified interior, 65 for undecided
  const auto classify = [&](C x) {
    for (int64_t n = 0; n < max_iter; n++) {
      if (std::norm(x) > 4) return 63 - __builtin_clzll(uint64_t(n) | 1);
      if (near) {
        if (std::norm(x - alpha) < rdisk * rdisk) return 64;
        x = x * x + c;
        continue;
      }
      for (int i = 0; i < P; i++) {
        const C w = x - z[i];
        if (std::norm(w) < rcert * rcert) {
          C wq = w;
          for (int j = 1; j < q; j++) wq *= w;
          if ((k[i] / wq).real() > R) return 64;
        }
      }
      x = x * x + c;
    }
    return 65;
  };

  // Local mode: scales from r0 down to where the escape time reaches max_iter/4
  double amax = 0;
  for (const C x : a) amax = std::max(amax, std::abs(x));
  const double r_min = std::pow(4.0 * P / (amax * max_iter / q), 1.0 / q), L = std::log(r0 / r_min);
  if (local) printf("local: r in [%.3g, %.3g] around each cycle point\n", r_min, r0);
  std::vector<std::array<double, 66>> acc(T), acc2(T);
  for (int t = 0; t < T; t++) { acc[t].fill(0); acc2[t].fill(0); }
  std::vector<std::thread> pool;
  for (int t = 0; t < T; t++)
    pool.emplace_back([&, t]() {
      std::mt19937_64 rng(4242 + 97 * t);
      std::uniform_real_distribution<double> u(0, 1);
      for (int64_t smp = t; smp < S; smp += T) {
        C x;
        double weight;
        if (local) {
          const int i = int(u(rng) * P) % P;
          const double rr = r0 * std::exp(-L * u(rng));
          x = z[i] + std::polar(rr, 2 * M_PI * u(rng));
          weight = rr * rr * L * 2 * M_PI * P;  // dA = r^2 dt dθ, one of P disks
        } else {
          do { x = C(4 * u(rng) - 2, 4 * u(rng) - 2); } while (std::norm(x) >= 4);
          weight = 4 * M_PI;
        }
        const int out = classify(x);
        acc[t][out] += weight;
        acc2[t][out] += weight * weight;
      }
    });
  for (auto& th : pool) th.join();
  std::array<double, 66> A{}, A2{};
  for (int t = 0; t < T; t++) for (int i = 0; i < 66; i++) { A[i] += acc[t][i]; A2[i] += acc2[t][i]; }
  printf("%s, %lld samples, max_iter 2^%d: in K %.6f, undecided %.2e (area)\n", argv[1],
         (long long)S, lmi, A[64] / S, A[65] / S);
  printf(" j  D_j          ±         B_j = D_j 2^(js) / (1 - 2^-s) ±    D_j / D_(j-1) (pred %.3f)\n", std::pow(2.0, -s));
  for (int j = 0; j < 64; j++) {
    if (!A[j]) continue;
    const double D = A[j] / S, e = std::sqrt(std::max(0.0, A2[j] / S - D * D) / S);
    const double B = D * std::pow(2.0, j * s) / (1 - std::pow(2.0, -s));
    printf("%2d  %.5e  %.1e   %.5f ± %.5f   %s\n", j, D, e, B, e / D * B,
           j && A[j - 1] ? std::to_string(A[j] / A[j - 1]).c_str() : "");
  }
}
