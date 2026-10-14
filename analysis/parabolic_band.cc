// Local band constants of parabolic points: the area near α escaping after n steps, T_local(n) ≈ B n^-(1 + 2/q).
// α is the parabolic fixed point at the root of the p/q bulb (p/q = 0/1: the cusp c = 1/4), f^q(α + w) =
// α + w + a w^{q+1} + ....  Sample w log-uniformly in |w| ∈ [r_min, r0] with uniform angle, weight by area
// (|w|^2 log(r0/r_min) 2π per sample), iterate until escape (|z| > 2) or the attracting petals Re(-1/(q a w^q)) > R,
// and estimate B per octave as D_j 2^{j s} / (1 - 2^-s), s = 1 + 2/q.  Also |E0| = B |a|^{2/q} (q + 2)/q, the
// exterior area per unit length of the repelling Écalle cylinder.
//   ./build/release/parabolic_band p q samples_log2 max_iter_log2 threads [r0]
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

static std::vector<C> compose_q(const C lam, const int q, const int D) {
  std::vector<C> s(D + 1, 0.0);
  s[1] = 1;
  for (int it = 0; it < q; it++) {
    std::vector<C> t(D + 1, 0.0);
    for (int i = 0; i <= D; i++) t[i] += lam * s[i];
    for (int i = 0; i <= D; i++) for (int j = 0; i + j <= D; j++) t[i + j] += s[i] * s[j];
    s = t;
  }
  return s;
}

int main(int argc, char** argv) {
  const int pp = atoi(argv[1]), q = atoi(argv[2]), lmi = atoi(argv[4]), T = atoi(argv[5]);
  const int64_t S = int64_t(1) << atoi(argv[3]), max_iter = int64_t(1) << lmi;
  const double r0 = argc > 6 ? atof(argv[6]) : 0.1;
  const C lam = std::polar(1.0, 2 * M_PI * pp / q), alpha = lam / 2.0, c = alpha - alpha * alpha;
  const auto ser = compose_q(lam, q, q + 2);
  const C a = ser[q + 1], k = -1.0 / (double(q) * a);
  const double s = 1 + 2.0 / q, R = 50;
  // Smallest scale: escape time ~ 1/(q |a| |w|^q) steps of f^q, so q/(q |a| |w|^q) = max_iter / 4
  const double r_min = std::pow(4.0 / (std::abs(a) * max_iter), 1.0 / q), L = std::log(r0 / r_min);
  printf("p/q %d/%d: c = %.12f %+.12fi, a = %.6f %+.6fi (|a| %.6f), r in [%.3g, %.3g]\n", pp, q, c.real(), c.imag(),
         a.real(), a.imag(), std::abs(a), r_min, r0);
  std::vector<std::array<double, 66>> acc(T), acc2(T);
  for (int t = 0; t < T; t++) { acc[t].fill(0); acc2[t].fill(0); }
  std::vector<std::thread> pool;
  for (int t = 0; t < T; t++)
    pool.emplace_back([&, t]() {
      std::mt19937_64 rng(999 + 31 * t);
      std::uniform_real_distribution<double> u(0, 1);
      for (int64_t smp = t; smp < S; smp += T) {
        const double rr = r0 * std::exp(-L * u(rng)), th = 2 * M_PI * u(rng);
        const double weight = rr * rr * L * 2 * M_PI;  // dA = r^2 dt dθ
        C z = alpha + std::polar(rr, th);
        int out = 65;
        for (int64_t n = 0; n < max_iter; n++) {
          if (std::norm(z) > 4) { out = 63 - __builtin_clzll(uint64_t(n) | 1); break; }
          const C w = z - alpha;
          if (std::norm(w) < 0.01) {
            C wq = w;
            for (int i = 1; i < q; i++) wq *= w;
            if ((k / wq).real() > R) { out = 64; break; }
          }
          z = z * z + c;
        }
        acc[t][out] += weight;
        acc2[t][out] += weight * weight;
      }
    });
  for (auto& th : pool) th.join();
  std::array<double, 66> A{}, A2{};
  for (int t = 0; t < T; t++) for (int i = 0; i < 66; i++) { A[i] += acc[t][i]; A2[i] += acc2[t][i]; }
  printf("disk area %.6f: in K %.6f, undecided %.2e\n", M_PI * r0 * r0, A[64] / S, A[65] / S);
  printf(" j  D_j          ±        B_j = D_j 2^(js) / (1 - 2^-s)    |E0|_j\n");
  for (int j = 0; j < 64; j++) {
    if (!A[j]) continue;
    const double D = A[j] / S, e = std::sqrt(std::max(0.0, A2[j] / S - D * D) / S);
    const double B = D * std::pow(2.0, j * s) / (1 - std::pow(2.0, -s)), E0 = B * std::pow(std::abs(a), 2.0 / q) * (q + 2) / q;
    printf("%2d  %.5e  %.1e   %.5f ± %.5f   %.5f\n", j, D, e, B, e / D * B, E0);
  }
}
