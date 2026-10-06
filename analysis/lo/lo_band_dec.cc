// As lo_band.cc, with η log-uniform in (ETA_LO, ETA), reporting misfire area by decade of η
//   clang++ -std=c++20 -O3 -mcpu=native -DMANDELBROT_ORBIT64 -I. -Ibuild/release analysis/lo/lo_band_dec.cc -o lo_band_dec  (from the repo root)
//   ./lo_band_dec KIND ETA_LO ETA N  (KIND 1: cardioid, 2: period-2 disk)
#include "orbit.h"
#include <atomic>
#include <cmath>
#include <cstdio>
#include <random>
#include <thread>
#include <vector>
#include <array>
#include <complex>
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
  // Exterior bands along the main cardioid (c = λ/2 - λ²/4) or the period-2 disk (c = -1 + λ/4), |λ| = 1 + η with
  // η log-uniform in (eta_lo, eta), θ uniform in (0, π) (the upper half plane): each sample weighs its area element
  // |c'(λ)|² |λ| dη dθ, so weighted misfire sums are areas (doubled for both half planes, as Lo's estimate)
  const int kind = atoi(argv[1]);  // 1: cardioid, 2: period-2 disk
  const double eta_lo = atof(argv[2]), eta = atof(argv[3]);
  const int64_t N = atoll(argv[4]);
  const int64_t max_iter = (int64_t(1) << 26) + 8;
  NewtonOptions nw; nw.repel2 = 2; nw.period = false;
  std::atomic<int64_t> next(0);
  const int T = std::thread::hardware_concurrency();
  std::vector<std::array<double, 6>> acc(T);
  std::vector<std::array<double, 10>> dec(T);
  for (auto& d : dec) d.fill(0);
  for (auto& a : acc) a.fill(0);
  auto worker = [&](int t) {
    std::mt19937_64 rng(4321 + 7 * t + kind);
    std::uniform_real_distribution<double> u(0, 1);
    auto& a = acc[t];
    for (;;) {
      if (next++ >= N) return;
      const double le = std::log(eta_lo) + u(rng) * (std::log(eta) - std::log(eta_lo)), e = std::exp(le), r = 1 + e, th = M_PI * u(rng);
      const std::complex<double> l = std::polar(r, th);
      std::complex<double> c, dc;
      if (kind == 1) { c = l / 2.0 - l * l / 4.0; dc = 0.5 - l / 2.0; }
      else { c = -1.0 + l / 4.0; dc = 0.25; }
      const double wgt = std::norm(dc) * r * e * (std::log(eta) - std::log(eta_lo)) * M_PI;  // Area element / sampling density (log-uniform η in (eta_lo, eta))
      const double x = c.real(), y = c.imag();
      a[0] += wgt;
      const int res = lo_membership(x, y, uint64_t(1) << 30);
      if (res == NOT_A_MEMBER) { a[3] += wgt; continue; }  // Escapes in Lo's own orbit: correct
      if (res == UNDECIDED) { a[2] += wgt; continue; }
      // Lo says member: is it?  Our reference decides (interior, escapes, or unknown by max_iter)
      Orbit<double> o;
      if (!o.start(x, y, 8192)) o.finish(max_iter, 256, nw);
      if (o.status == 2) { a[1] += wgt; continue; }
      if (o.status != 1) { a[2] += wgt; continue; }
      a[4] += wgt; a[5] += 1; dec[t][std::min(9, int(-std::log10(e)))] += wgt;
    }
  };
  std::vector<std::thread> pool;
  for (int t = 0; t < T; t++) pool.emplace_back(worker, t);
  for (auto& t : pool) t.join();
  std::array<double, 6> s{};
  for (auto& a : acc) for (int i = 0; i < 6; i++) s[i] += a[i];
  const double n = double(N);
  printf("%s band eta %.0e..%.0e: band area %.3e, Lo-member interior %.3e, unknown %.3e, Lo-escaping %.3e; Lo misfires %.0f samples, "
         "false-member area %.3e (doubled)\n", kind == 1 ? "cardioid" : "disk", eta_lo, eta, 2 * s[0] / n, 2 * s[1] / n,
         2 * s[2] / n, 2 * s[3] / n, s[5], 2 * s[4] / n);
  printf("  misfire area by decade of eta (10^-j):");
  for (int j = 0; j < 10; j++) { double d = 0; for (auto& x : dec) d += x[j]; if (d) printf(" j=%d %.2e", j, 2 * d / n); }
  printf("\n");
}
