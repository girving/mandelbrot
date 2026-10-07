// Escape-time tails of filled Julia sets at and near parabolic parameters (see julia_tail.h), on CPU or GPU.
//   global: T(n) ~ n^-(1 + 2/q) over all of K's complement in |z| ≤ 2
//   local: the band constant B in T_local(n) ≈ B n^-(1 + 2/q) near the parabolic cycle (summed over its points)
//   near: c inside the main cardioid near a root, for the crossover to exponential decay
// The ratio of global to local constants is the total weight W of the cycle's preimages.
//   ./build/release/parabolic_tail global|local|near c_re c_im z_re z_im P q samples_log2 max_iter_log2 [r0] [--cuda]
// with (z_re, z_im) a guess for a cycle point (unused by near).  CPU threads: $MANDELBROT_THREADS.
#include "julia_tail.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
using namespace mandelbrot;

int main(int argc, char** argv) {
  bool cuda = false;
  std::vector<char*> args;
  for (int i = 0; i < argc; i++) {
    if (!strcmp(argv[i], "--cuda")) cuda = true;
    else args.push_back(argv[i]);
  }
  if (args.size() < 10) {
    fprintf(stderr, "usage: parabolic_tail global|local|near c_re c_im z_re z_im P q samples_log2 max_iter_log2 [r0] [--cuda]\n");
    return 1;
  }
  TailParams p;
  p.mode = !strcmp(args[1], "local") ? TailMode::local : !strcmp(args[1], "near") ? TailMode::near : TailMode::global;
  p.c = {atof(args[2]), atof(args[3])};
  p.z_guess = {atof(args[4]), atof(args[5])};
  p.P = atoi(args[6]); p.q = atoi(args[7]);
  p.samples = int64_t(1) << atoi(args[8]);
  p.max_iter = int64_t(1) << atoi(args[9]);
  if (args.size() > 10) p.r0 = atof(args[10]);
  p.cuda = cuda;
  const auto r = julia_tail(p);
  const int q = p.mode == TailMode::near ? p.q : p.q;
  const double s = 1 + 2.0 / q;
  printf("%s c = %.15f %+.15fi, P %d, q %d, %lld samples, max_iter 2^%d, %s: %.1f s, %.3g iterations\n", args[1],
         p.c.real(), p.c.imag(), p.P, q, (long long)p.samples, atoi(args[9]), cuda ? "gpu" : "cpu", r.secs,
         double(r.iters));
  for (size_t i = 0; i < r.cycle.size(); i++)
    printf("  z_%zu = %.12f %+.12fi, a = %.6f %+.6fi (|a| %.6g)\n", i, r.cycle[i].real(), r.cycle[i].imag(),
           r.a[i].real(), r.a[i].imag(), std::abs(r.a[i]));
  if (p.mode == TailMode::local) printf("  shells: %d, half-sides %.3g … %.3g\n", r.bins, r.shell[0], r.shell[r.bins]);
  printf("in K %.6f ± %.1e, undecided %.2e (area)\n", r.area(64), r.error(64), r.area(65));
  printf(" j  D_j          ±         B_j = D_j 2^(js) / (1 - 2^-s) ±    D_j / D_(j-1) (pred %.3f)\n", std::pow(2.0, -s));
  for (int j = 0; j < 64; j++) {
    const double D = r.area(j);
    if (!D) continue;
    const double e = r.error(j), B = D * std::pow(2.0, j * s) / (1 - std::pow(2.0, -s));
    const double prev = j ? r.area(j - 1) : 0;
    printf("%2d  %.5e  %.1e   %.5f ± %.5f   %s\n", j, D, e, B, e / D * B, prev ? std::to_string(D / prev).c_str() : "");
  }
}
