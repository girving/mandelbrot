// Tests for escape-time tails of filled Julia sets

#include "julia_tail.h"
#include "cutil.h"
#include "tests.h"
#include <cmath>
namespace mandelbrot {
namespace {

TEST(global) {
  // c = -3/4: every sample is classified, K's area is ~2.103, and the tail approaches n^-2 per octave
  TailParams p;
  p.c = -0.75; p.z_guess = -0.5; p.P = 1; p.q = 2;
  p.samples = 1 << 26; p.max_iter = 1 << 16;
  const auto r = julia_tail(p);
  int64_t total = 0;
  for (const auto x : r.counts) total += x;
  ASSERT_EQ(total, p.samples);
  const double ratio = r.area(8) / r.area(7);
  print("  c = -3/4: in K %.4f ± %.4f, undecided %.2e, D_8/D_7 %.3f (n^-2: 0.25, 0.21-0.23 this early), %.2f s",
        r.area(64), r.error(64), r.area(65), ratio, r.secs);
  ASSERT_LE(std::abs(r.area(64) - 2.103), 0.005);
  ASSERT_LE(0.17, ratio);
  ASSERT_LE(ratio, 0.28);
}

TEST(local) {
  // The cusp's local band constant B in T_local(n) ≈ B n^-3 is ~1.55
  TailParams p;
  p.mode = TailMode::local;
  p.c = 0.25; p.z_guess = 0.5; p.P = 1; p.q = 1;
  p.samples = 1 << 22; p.max_iter = 1 << 18;
  const auto r = julia_tail(p);
  const double B = r.area(10) * std::pow(2.0, 30) / (1 - 0.125);
  print("  cusp: %d shells, B at octave 10 %.3f ± %.3f, %.2f s", r.bins, B, B * r.error(10) / r.area(10), r.secs);
  ASSERT_LE(std::abs(B - 1.58), 0.25);
}

TEST(near) {
  // Near -3/4 inside the cardioid every sample is classified, and the area is near K(-3/4)'s
  TailParams p;
  p.mode = TailMode::near;
  p.c = -0.748;
  p.samples = 1 << 20; p.max_iter = 1 << 20;
  const auto r = julia_tail(p);
  print("  c = -0.748: in K %.4f, undecided %.2e", r.area(64), r.area(65));
  ASSERT_EQ(r.area(65), 0);
  ASSERT_LE(std::abs(r.area(64) - 2.11), 0.03);
}

TEST(cuda_matches_cpu) {
  IF_CUDA({
    for (const auto mode : {TailMode::global, TailMode::local, TailMode::near}) {
      TailParams p;
      p.mode = mode;
      p.c = mode == TailMode::near ? -0.748 : -0.75; p.z_guess = -0.5; p.P = 1; p.q = 2;
      p.samples = 1 << 22; p.max_iter = 1 << 18;
      const auto cpu = julia_tail(p);
      p.cuda = true;
      const auto gpu = julia_tail(p);
      print("  mode %d: cpu %.2f s, gpu %.2f s", int(mode), cpu.secs, gpu.secs);
      ASSERT_EQ(cpu.counts, gpu.counts);
    }
  })
}

}  // namespace
}  // namespace mandelbrot
