// Bias of Hsing Lo's membership test (github.com/hsingtism/mandelbrot-area, mem-func.c at 3ea0dd8): the
// probability that it calls an escaping point a member, by octave of the escape step.  Samples are uniform in
// leaf cells from escape_tree --dump (base 1000, depth 10, default box); for each, our Orbit<double> finds the
// escape step n (or settles it as interior / undecided, which are skipped: Lo's errors are on escaping points),
// then Lo's test runs on the same c.  A misfire is Lo returning MEMBER for a point that escapes.
//   clang++ -std=c++20 -O3 -mcpu=native -DMANDELBROT_ORBIT64 -I. -Ibuild/release analysis/lo/lo_bias.cc -o lo_bias  (from the repo root)
//   ./lo_bias DUMP... samples_per_leaf max_leaves seed
#include "orbit.h"
#include <atomic>
#include <cmath>
#include <cstdio>
#include <random>
#include <thread>
#include <vector>
using namespace mandelbrot;

enum { MEMBER_CONV, MEMBER_ORBIT, NOT_A_MEMBER, UNDECIDED };

// Lo's membership(), verbatim apart from returning which test fired (and a dwell cap: escaping points end anyway)
static int lo_membership(double re, double im, const uint64_t dwell) {
  if (im < -1.15 || im > 1.15 || re < -2.0 || re > 0.49 || re * re + im * im > 4.0) return NOT_A_MEMBER;
  if (im < 0.25 && re < -0.75 && re > -1.25) {
    const double xp1 = re + 1;
    if (std::signbit(xp1 * xp1 + im * im - 0.0625)) return MEMBER_CONV;
  } else if (re > -0.75 && im < 0.65 && re < 0.375) {
    const double adjx = re - 0.25;
    const double adjx2py2 = adjx * adjx + im * im;
    const double firstterm = adjx2py2 + 2 * 0.25 * adjx;
    if (std::signbit(firstterm * firstterm - 0.25 * adjx2py2)) return MEMBER_CONV;
  }
  const double C = std::ldexp(1.0, -16), O = std::ldexp(1.0, -32);
  const double cRe = re, cIm = im;
  double pRe, pIm, pobRe;
  double obRe = re, obIm = im;
  for (uint64_t i = 0; i < dwell; i++) {
    pRe = re;
    pIm = im;
    re = re * re - im * im + cRe;
    im = 2.0 * pRe * im + cIm;
    if (i % 5 == 1) {
      if (re * re + im * im > 4.0) return NOT_A_MEMBER;
      if (std::fabs(pRe - re) < C && std::fabs(pIm - im) < C) return MEMBER_CONV;
    }
    if (i % 2) {
      pobRe = obRe;
      obRe = obRe * obRe - obIm * obIm + cRe;
      obIm = 2 * pobRe * obIm + cIm;
      if (std::fabs(obRe - re) < O && std::fabs(obIm - im) < O) return MEMBER_ORBIT;
    }
  }
  return UNDECIDED;
}

int main(int argc, char** argv) {
  const int nd = argc - 4;
  const int S = atoi(argv[argc - 3]);
  const int64_t max_leaves = atoll(argv[argc - 2]);
  const uint64_t seed = atoll(argv[argc - 1]);
  std::vector<std::pair<int32_t, int32_t>> leaves;
  for (int f = 0; f < nd; f++) {
    FILE* fp = fopen(argv[1 + f], "rb");
    uint32_t buf[18];
    while (fread(buf, 4, 18, fp) == 18) leaves.push_back({int32_t(buf[0]), int32_t(buf[1])});
    fclose(fp);
  }
  std::mt19937_64 pick(seed);
  std::shuffle(leaves.begin(), leaves.end(), pick);
  if (int64_t(leaves.size()) > max_leaves) leaves.resize(max_leaves);
  const int64_t L = leaves.size();
  const double w = 2.5 / (1000.0 * 1024), h = 1.2 / (1000.0 * 1024);
  const int64_t max_iter = (int64_t(1) << 28) + 8;
  NewtonOptions nw;
  nw.repel2 = 2;
  nw.period = false;
  constexpr int OCT = 40;
  std::atomic<int64_t> next(0);
  std::vector<std::array<int64_t, 4 * OCT>> acc(std::thread::hardware_concurrency());
  for (auto& a : acc) a.fill(0);
  std::atomic<int64_t> interior(0), undecided(0), total(0);
  auto worker = [&](const int t) {
    auto& a = acc[t];
    for (;;) {
      const int64_t i = next++;
      if (i >= L) return;
      std::mt19937_64 rng(seed * 7919 + i);
      std::uniform_real_distribution<double> u(0, 1);
      for (int s = 0; s < S; s++) {
        const double x = -2 + (leaves[i].first + u(rng)) * w, y = (leaves[i].second + u(rng)) * h;
        Orbit<double> o;
        total++;
        if (o.start(x, y, 8192)) { interior++; continue; }
        o.finish(max_iter, 256, nw);
        if (o.status == 2) { interior++; continue; }
        if (o.status != 1) { undecided++; continue; }
        const int oct = int(std::floor(std::log2(double(std::max<int64_t>(o.n, 1)))));
        const int r = lo_membership(x, y, uint64_t(1) << 31);  // Lo's orbit escapes at its own (re-randomized) step
        a[4 * oct + 0]++;                                 // escaping samples in this octave
        if (r == MEMBER_CONV) a[4 * oct + 1]++;           // Lo misfire: convergence test
        if (r == MEMBER_ORBIT) a[4 * oct + 2]++;          // Lo misfire: orbit test
        if (r == UNDECIDED) a[4 * oct + 3]++;             // should not happen for escaping points
      }
    }
  };
  std::vector<std::thread> pool;
  for (unsigned t = 0; t < acc.size(); t++) pool.emplace_back(worker, t);
  for (auto& t : pool) t.join();
  printf("%lld leaves, %lld samples: %lld interior, %lld undecided at 2^28\n", (long long)L, (long long)total.load(),
         (long long)interior.load(), (long long)undecided.load());
  printf("octave  escaping  conv_misfire  orbit_misfire  undecided\n");
  for (int oct = 0; oct < OCT; oct++) {
    int64_t e = 0, c = 0, ob = 0, un = 0;
    for (auto& a : acc) { e += a[4 * oct]; c += a[4 * oct + 1]; ob += a[4 * oct + 2]; un += a[4 * oct + 3]; }
    if (e) printf("%2d %12lld %10lld %10lld %6lld\n", oct, (long long)e, (long long)c, (long long)ob, (long long)un);
  }
}
