// Escape-time tails of filled Julia sets at and near parabolic parameters

#include "julia_tail.h"
#include "debug.h"
#include "engine.h"
#include <chrono>
#include <cmath>
namespace mandelbrot {

using std::complex;

namespace {

constexpr int kMaxP = 8, kMaxBins = 48, kOut = 66, kCopies = 1024;

__host__ __device__ inline int octave(const int64_t n) {
#ifdef __CUDA_ARCH__
  return 63 - __clzll(n | 1);
#else
  return 63 - __builtin_clzll(uint64_t(n | 1));
#endif
}

struct TailState {
  double x, y;
  int64_t n;
  int32_t out;  // -1 while running
};

struct TailTask {
  typedef TailState State;
  double cx, cy;
  int mode, P, q, bins;
  double zx[kMaxP], zy[kMaxP], kx[kMaxP], ky[kMaxP];  // Cycle, and k = -1/(q a) per point
  double rcert2, R, ax, ay, rdisk2;
  double shell[kMaxBins + 1];
  int64_t max_iter, key0;
  uint64_t seed;
  uint64_t* hist;  // kCopies × bins × kOut
  int64_t burst = 256;
  int min_blocks = 2;

  __host__ __device__ int bin(const int64_t key) const { return mode == 1 ? int(key % bins) : 0; }

  __host__ __device__ bool start(State& o, const int64_t i) const {
    const int64_t key = key0 + i;
    o.n = 0;
    o.out = -1;
    if (mode == 1) {
      // A square shell around one cycle point, by rejection from its outer square
      const int b = bin(key), pt = int((key / bins) % P);
      const double s = shell[b], s1 = shell[b + 1];
      double x = 0, y = 0;
      for (int j = 0; j < 128; j += 2) {
        x = s * (2 * uniform(seed, key, j) - 1);
        y = s * (2 * uniform(seed, key, j + 1) - 1);
        if (fmax(fabs(x), fabs(y)) > s1) break;
      }
      o.x = zx[pt] + x;
      o.y = zy[pt] + y;
    } else {
      o.x = 4 * uniform(seed, key, 0) - 2;
      o.y = 4 * uniform(seed, key, 1) - 2;
    }
    return false;
  }

  __host__ __device__ bool run(State& o) const {
    double x = o.x, y = o.y;
    int64_t n = o.n;
    const int64_t end = n + burst < max_iter ? n + burst : max_iter;
    for (; n < end; n++) {
      const double x2 = x * x, y2 = y * y;
      if (x2 + y2 > 4) { o.out = octave(n); break; }
      if (mode == 2) {
        const double dx = x - ax, dy = y - ay;
        if (dx * dx + dy * dy < rdisk2) { o.out = 64; break; }
      } else {
        bool in = false;
        for (int k = 0; k < P; k++) {
          const double wx = x - zx[k], wy = y - zy[k];
          if (wx * wx + wy * wy < rcert2) {
            // Re(k / w^q) > R  ⟺  Re(k conj(w^q)) > R |w^q|^2
            double px = wx, py = wy;
            for (int j = 1; j < q; j++) {
              const double t = px * wx - py * wy;
              py = px * wy + py * wx;
              px = t;
            }
            if (kx[k] * px + ky[k] * py > R * (px * px + py * py)) in = true;
          }
        }
        if (in) { o.out = 64; break; }
      }
      const double t = x2 - y2 + cx;
      y = 2 * x * y + cy;
      x = t;
    }
    o.x = x; o.y = y; o.n = n;
    if (o.out < 0 && n >= max_iter) o.out = 65;
    return o.out >= 0;
  }

  __host__ __device__ void finish(const State& o, const int64_t i) const {
    const int64_t key = key0 + i;
    atomic_add(hist + (int64_t(i % kCopies) * bins + bin(key)) * kOut + o.out, 1);
  }
  __host__ __device__ int64_t iters(const State& o) const { return o.n; }
  __host__ __device__ int64_t progress(const State& o) const { return o.n; }
  __host__ __device__ bool pending(const State&) const { return false; }
  __host__ __device__ bool settle(State&) const { return true; }
};

// Taylor coefficients (degree ≤ D) of f^m(z_k + w) - z_{k+m} in w, composing f(z_j + w) = z_{j+1} + 2 z_j w + w^2
vector<complex<double>> germ(const vector<complex<double>>& z, const int k, const int m, const int D) {
  vector<complex<double>> s(D + 1, 0.0);
  s[1] = 1;
  const int P = int(z.size());
  for (int it = 0; it < m; it++) {
    const complex<double> zj = z[(k + it) % P];
    vector<complex<double>> t(D + 1, 0.0);
    for (int i = 0; i <= D; i++) t[i] += 2.0 * zj * s[i];
    for (int i = 0; i <= D; i++) for (int j = 0; i + j <= D; j++) t[i + j] += s[i] * s[j];
    s = t;
  }
  return s;
}

}  // namespace

double TailResult::area(const int j) const {
  double sum = 0;
  for (int b = 0; b < bins; b++) sum += counts[b * kOut + j] * weight[b];
  return sum;
}

double TailResult::error(const int j) const {
  // Per bin, counts are binomial out of the bin's samples: var = w^2 c (1 - c/n) ≈ w^2 c
  double var = 0;
  for (int b = 0; b < bins; b++) var += double(counts[b * kOut + j]) * weight[b] * weight[b];
  return std::sqrt(var);
}

TailResult julia_tail(const TailParams& p) {
  const auto t0 = std::chrono::steady_clock::now();
  TailResult res;
  TailTask task{};
  const complex<double> c = p.c;
  task.cx = c.real(); task.cy = c.imag();
  task.mode = int(p.mode);
  task.max_iter = p.max_iter;
  task.seed = p.seed;
  task.R = p.R;
  if (p.mode == TailMode::near) {
    const complex<double> alpha = (1.0 - std::sqrt(1.0 - 4.0 * c)) / 2.0;
    const double lam = std::abs(2.0 * alpha);
    slow_assert(lam < 1, "near: c is not inside the main cardioid");
    task.ax = alpha.real(); task.ay = alpha.imag();
    task.rdisk2 = (1 - lam) * (1 - lam) / 4;
    task.P = 0;
  } else {
    // Refine the cycle: Newton on f^P(z) - z for q > 1 (a simple fixed point of f^P), and on (f^P)'(z) - 1 for
    // q = 1 (where the fixed point is double but the multiplier condition is simple)
    const int P = p.P, q = p.q;
    slow_assert(1 <= P && P <= kMaxP && q >= 1, "need 1 ≤ P ≤ %d", kMaxP);
    complex<double> z0 = p.z_guess;
    for (int it = 0; it < 100; it++) {
      complex<double> x = z0, d = 1, dd = 0;
      for (int i = 0; i < P; i++) { dd = 2.0 * (d * d + x * dd); d = 2.0 * x * d; x = x * x + c; }
      const complex<double> step = q > 1 ? (x - z0) / (d - 1.0) : (d - 1.0) / dd;
      z0 -= step;
      if (std::abs(step) < 1e-15) break;
    }
    res.cycle.resize(P);
    res.cycle[0] = z0;
    for (int i = 1; i < P; i++) res.cycle[i] = res.cycle[i - 1] * res.cycle[i - 1] + c;
    double scale = 1;
    for (int i = 0; i < P; i++) {
      const auto g = germ(res.cycle, i, P * q, q + 2);
      res.a.push_back(g[q + 1]);
      const complex<double> k = -1.0 / (double(q) * g[q + 1]);
      task.zx[i] = res.cycle[i].real(); task.zy[i] = res.cycle[i].imag();
      task.kx[i] = k.real(); task.ky[i] = k.imag();
      scale = std::min(scale, std::pow(std::abs(g[q + 1]), -1.0 / q));
    }
    task.P = P; task.q = q;
    task.rcert2 = 0.3 * scale * 0.3 * scale;
    if (p.mode == TailMode::local) {
      // Shells from r0 down to where the escape time ~P / (|a| r^q) reaches max_iter / 4
      double amax = 0;
      for (const auto x : res.a) amax = std::max(amax, std::abs(x));
      const double r0 = p.r0 ? p.r0 : 0.3 * scale;
      const double r_min = std::pow(4.0 * P / (amax * double(p.max_iter)), 1.0 / q);
      slow_assert(r_min < r0, "local: r_min %g ≥ r0 %g; raise max_iter or r0", r_min, r0);
      res.bins = std::min(kMaxBins, int(std::ceil(std::log2(r0 / r_min))));
      for (int b = 0; b <= res.bins; b++) res.shell.push_back(std::ldexp(r0, -b));
      for (int b = 0; b <= res.bins; b++) task.shell[b] = res.shell[b];
    }
  }
  task.bins = res.bins;

  // Histograms in kCopies replicas, so that atomics rarely collide
  const int64_t hn = int64_t(kCopies) * res.bins * kOut;
  Mem<uint64_t> hist(hn, p.cuda);
  hist.zero();
  task.hist = hist.p;
  const int64_t batch = int64_t(1) << 30;
  for (int64_t done = 0; done < p.samples; done += batch) {
    task.key0 = done;
    res.iters += run_orbits(task, std::min(batch, p.samples - done), p.cuda).iters;
  }
  vector<uint64_t> h(hn);
  hist.to_host(h.data(), hn);
  res.counts.assign(res.bins * kOut, 0);
  for (int64_t r = 0; r < kCopies; r++)
    for (int64_t j = 0; j < res.bins * kOut; j++) res.counts[j] += int64_t(h[r * res.bins * kOut + j]);

  // Area per sample: 16 / samples over the square; per shell, its area 3 s_b^2 (2s outer side, s inner) times
  // P points over the samples landing in it (keys ≡ b mod bins)
  res.weight.resize(res.bins);
  if (p.mode == TailMode::local)
    for (int b = 0; b < res.bins; b++) {
      const int64_t nb = p.samples / res.bins + (b < p.samples % res.bins);
      const double s = res.shell[b], s1 = res.shell[b + 1];
      res.weight[b] = 4 * (s * s - s1 * s1) * p.P / double(nb);
    }
  else
    res.weight[0] = 16.0 / double(p.samples);
  res.secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return res;
}

}  // namespace mandelbrot
