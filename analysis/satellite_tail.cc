// Escape-time tails of K(c) at the satellite roots c of the p/q bulbs: sample z uniformly in |z| < 2, iterate
// until |z| > 2 (escaped) or z is in the attracting petals Re u > R, u = -1/(q a w^q), w = z - α, where
// f^q(α + w) = α + w + a w^{q+1} + ... (u is invariant under w ↦ λw, and f^q acts as u ↦ u + 1 + ...).
// Histogram escape steps by octave.  Prediction from the Fatou coordinate: T(n) ~ n^-(1 + 2/q).
//   ./build/release/satellite_tail p q samples_log2 max_iter_log2 threads
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <thread>
#include <vector>
typedef std::complex<double> C;

// Truncated power series in w (degree ≤ D) composition: g(w) = λ w + w^2, compute g^q
static std::vector<C> compose_q(const C lam, const int q, const int D) {
  std::vector<C> s(D + 1, 0.0);
  s[1] = 1;  // Identity
  for (int it = 0; it < q; it++) {
    // t = λ s + s^2
    std::vector<C> t(D + 1, 0.0);
    for (int i = 0; i <= D; i++) t[i] += lam * s[i];
    for (int i = 0; i <= D; i++) for (int j = 0; i + j <= D; j++) t[i + j] += s[i] * s[j];
    s = t;
  }
  return s;
}

int main(int argc, char** argv) {
  const int pp = atoi(argv[1]), q = atoi(argv[2]);
  const int64_t S = int64_t(1) << atoi(argv[3]), max_iter = int64_t(1) << atoi(argv[4]);
  const int T = atoi(argv[5]);
  const C lam = std::polar(1.0, 2 * M_PI * pp / q), alpha = lam / 2.0, c = alpha - alpha * alpha;
  const auto s = compose_q(lam, q, q + 2);
  for (int i = 2; i <= q; i++) if (std::abs(s[i]) > 1e-12) { printf("unexpected term w^%d: %g\n", i, std::abs(s[i])); return 1; }
  const C a = s[q + 1], k = -1.0 / (double(q) * a);
  const double R = 50;
  printf("p/q %d/%d: c = %.15f %+.15fi, alpha = %.6f %+.6fi, f^q = w + (%.6f %+.6fi) w^%d + ...\n", pp, q, c.real(),
         c.imag(), alpha.real(), alpha.imag(), a.real(), a.imag(), q + 1);
  std::vector<std::array<int64_t, 66>> hist(T);
  for (auto& h : hist) h.fill(0);
  std::vector<std::thread> pool;
  for (int t = 0; t < T; t++)
    pool.emplace_back([&, t]() {
      std::mt19937_64 rng(12345 + 77 * t);
      std::uniform_real_distribution<double> u(-2, 2);
      auto& h = hist[t];
      for (int64_t smp = t; smp < S; smp += T) {
        C z;
        do { z = C(u(rng), u(rng)); } while (std::norm(z) >= 4);
        int out = 65;
        for (int64_t n = 0; n < max_iter; n++) {
          if (std::norm(z) > 4) { out = 63 - __builtin_clzll(uint64_t(n) | 1); break; }
          const C w = z - alpha;
          C wq = w;
          for (int i = 1; i < q; i++) wq *= w;
          if (std::norm(w) < 0.01 && (k / wq).real() > R) { out = 64; break; }
          z = z * z + c;
        }
        h[out]++;
      }
    });
  for (auto& th : pool) th.join();
  std::array<int64_t, 66> H{};
  for (auto& h : hist) for (int i = 0; i < 66; i++) H[i] += h[i];
  const double scale = 4 * M_PI / double(S), pred = 1 + 2.0 / q;
  printf("%lld samples, max_iter 2^%d: in K %.6f, undecided %.3e (area); prediction T(n) ~ n^-%.3f\n",
         (long long)S, atoi(argv[4]), H[64] * scale, H[65] * scale, pred);
  printf(" j  D_j        ±        D_j 2^(j (1+2/q))   D_j / D_(j-1) (pred %.3f)\n", std::pow(2.0, -pred));
  for (int j = 0; j < 64; j++) {
    if (!H[j]) continue;
    const double D = H[j] * scale, e = std::sqrt(double(H[j])) * scale;
    printf("%2d  %.4e  %.1e   %.5f   %s\n", j, D, e, D * std::pow(2.0, j * pred),
           j && H[j - 1] ? std::to_string(double(H[j]) / H[j - 1]).c_str() : "");
  }
}
