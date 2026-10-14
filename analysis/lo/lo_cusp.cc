// Lo's test against ours near the cusp c = 1/4: samples uniform in [1/4 - R, 1/4 + R] x [0, R] for several R
//   clang++ -std=c++20 -O3 -mcpu=native -DMANDELBROT_ORBIT64 -I. -Ibuild/release analysis/lo/lo_cusp.cc -o lo_cusp  (from the repo root)
//   ./lo_cusp R N
#include "orbit.h"
#include <atomic>
#include <cmath>
#include <cstdio>
#include <random>
#include <thread>
#include <vector>
using namespace mandelbrot;
enum { MEMBER_CONV, MEMBER_ORBIT, NOT_A_MEMBER, UNDECIDED };
static int lo_membership(double re, double im, const uint64_t dwell) {
  if (im < -1.15 || im > 1.15 || re < -2.0 || re > 0.49 || re * re + im * im > 4.0) return NOT_A_MEMBER;
  if (im < 0.25 && re < -0.75 && re > -1.25) {
    const double xp1 = re + 1;
    if (std::signbit(xp1 * xp1 + im * im - 0.0625)) return MEMBER_CONV;
  } else if (re > -0.75 && im < 0.65 && re < 0.375) {
    const double adjx = re - 0.25, adjx2py2 = adjx * adjx + im * im, firstterm = adjx2py2 + 2 * 0.25 * adjx;
    if (std::signbit(firstterm * firstterm - 0.25 * adjx2py2)) return 10;  // Cardioid shortcut (a true member)
  }
  const double C = std::ldexp(1.0, -16), O = std::ldexp(1.0, -32);
  const double cRe = re, cIm = im;
  double pRe, pIm, pobRe, obRe = re, obIm = im;
  for (uint64_t i = 0; i < dwell; i++) {
    pRe = re; pIm = im;
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
  const double R = atof(argv[1]);
  const int64_t N = atoll(argv[2]);
  const int64_t max_iter = (int64_t(1) << 30) + 8;
  NewtonOptions nw; nw.repel2 = 2; nw.period = false;
  std::atomic<int64_t> next(0), esc(0), mis_conv(0), mis_orb(0), lo_und(0), ours_int(0), ours_und(0), lo_mem_int(0);
  auto worker = [&](int t) {
    std::mt19937_64 rng(1234 + t);
    std::uniform_real_distribution<double> u(0, 1);
    for (;;) {
      if (next++ >= N) return;
      const double x = 0.25 - R + 2 * R * u(rng), y = R * u(rng);
      Orbit<double> o;
      if (!o.start(x, y, 8192)) o.finish(max_iter, 256, nw);
      if (o.status == 2) { ours_int++; continue; }
      if (o.status != 1) { ours_und++; continue; }
      esc++;
      const int r = lo_membership(x, y, uint64_t(1) << 32);
      if (r == MEMBER_CONV) mis_conv++;
      if (r == MEMBER_ORBIT) mis_orb++;
      if (r == 10) lo_mem_int++;
      if (r == UNDECIDED) lo_und++;
    }
  };
  std::vector<std::thread> pool;
  for (unsigned t = 0; t < std::thread::hardware_concurrency(); t++) pool.emplace_back(worker, t);
  for (auto& t : pool) t.join();
  const double box = 2 * R * R, scale = 2 * box / double(N);  // Doubled: Lo's estimate counts both half planes
  printf("R %.1e: %lld samples, ours interior %lld, undecided %lld, escaping %lld; Lo misfires: conv %lld, orbit %lld, "
         "cardioid-shortcut %lld, undecided %lld -> false-member area %.3e (doubled)\n", R, (long long)N,
         (long long)ours_int.load(), (long long)ours_und.load(), (long long)esc.load(), (long long)mis_conv.load(),
         (long long)mis_orb.load(), (long long)lo_mem_int.load(), (long long)lo_und.load(),
         scale * double(mis_conv + mis_orb + lo_mem_int));
}
