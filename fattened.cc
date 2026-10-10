// Exact areas of the fattened sets {c : g_M(c) < 2^-k} from the Böttcher coefficients (Grönwall at radius e^δ):
//   F(δ) = π e^{2δ} − π Σ_n n b_n² e^{−2nδ},  δ = 2^-k,
// for comparison with the tree's Monte Carlo estimates (escape_tree with the same k).  With coefficients to
// n = N, the omitted terms weigh at most e^{−2Nδ} times the remaining energy Σ_{n>N} n b_n² < 1; for
// N = 2^27 that is e^{−2^(28−k)}, negligible for k ≤ 22.

#include "debug.h"
#include "numpy.h"
#include "print.h"
#include "wall_time.h"
#include <cmath>
#include <vector>
namespace mandelbrot {
namespace {

void run(const string& path, const int k0, const int k1) {
  const auto t0 = wall_time();
  const auto F = read_numpy(path);
  slow_assert(F.shape.size() == 2 && F.shape[1] == 2, "bad coefficient file");
  const int64_t M = F.shape[0];  // f_m = b_{m-1} for m < M, as unevaluated sums of two doubles
  const int K = k1 - k0 + 1;
  vector<long double> sum(K);
  for (int64_t m = 2; m < M; m++) {
    const long double b = (long double)F.data[2 * m] + F.data[2 * m + 1], n = m - 1, x = n * b * b;
    for (int i = 0; i < K; i++) sum[i] += x * std::exp(-2 * n * std::ldexp(1.0L, -(k0 + i)));
  }
  print("read and summed %d coefficients: %.1f s", M, (wall_time() - t0).seconds());
  print("   k   area{g < 2^-k}     truncation weight e^(-2 N 2^-k)");
  for (int i = 0; i < K; i++) {
    const int k = k0 + i;
    const long double d = std::ldexp(1.0L, -k), pi = 3.141592653589793238462643383279502884L;
    print("  %2d   %.13f   %.1e", k, double(pi * std::exp(2 * d) - pi * sum[i]), double(std::exp(-2.0L * (M - 1) * d)));
  }
}

}  // namespace
}  // namespace mandelbrot

int main(const int argc, const char** argv) {
  using namespace mandelbrot;
  try {
    slow_assert(argc >= 2, "usage: %s f-kK.npy [k0] [k1]", argv[0]);
    run(argv[1], argc > 2 ? atoi(argv[2]) : 8, argc > 3 ? atoi(argv[3]) : 22);
    return 0;
  } catch (const std::exception& e) {
    die(e.what());
  }
}
