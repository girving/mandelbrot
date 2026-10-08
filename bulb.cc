// Areas of satellite bulbs, batched on CPU threads or the GPU

#include "bulb.h"
#include "cutil.h"
#include "debug.h"
#include "engine.h"
#include "expansion_arith.h"
#include "nearest.h"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <map>
#include <numeric>
#include <thread>
namespace mandelbrot {

namespace {

typedef Complex<double> Cd;
typedef Complex<E2> Ce;

template<class S> __host__ __device__ static inline Complex<S> cdiv(const Complex<S> a, const Complex<S> b) {
  const S d = sqr(b.r) + sqr(b.i);
  const Complex<S> n = a * conj(b);
  return Complex<S>(n.r / d, n.i / d);
}
// sqrt is correctly rounded on host and device (hypot is not), so CPU and GPU take identical paths
__host__ __device__ static inline double cabs(const Cd z) { return sqrt(z.r * z.r + z.i * z.i); }
__host__ __device__ static inline Cd scale(const double a, const Cd z) { return Cd(a * z.r, a * z.i); }
__host__ __device__ static inline Ce to_e(const Cd z) { return Ce(E2(z.r), E2(z.i)); }
__host__ __device__ static inline Cd to_d(const Ce z) { return Cd(double(z.r), double(z.i)); }

// One Newton step for (z, c): f^n(z) = z, (f^n)'(z) = mu with f(x) = x² + s x + c.  Returns |dz| + |dc| (double
// norms) and sets dcdmu = c'(μ).
template<class S> __host__ __device__ static double newton_step(const int n, const int shift, const Complex<S> mu,
                                                                Complex<S>& z, Complex<S>& c, Complex<S>& dcdmu) {
  typedef Complex<S> C;
  const C s(shift), one(1);
  C x = z, xz(1), xc(0), xzz(0), xzc(0);
  for (int i = 0; i < n; i++) {
    const C df = twice(x) + s;
    const C nxzz = twice(xz * xz) + df * xzz, nxzc = twice(xc * xz) + df * xzc;
    xzz = nxzz; xzc = nxzc;
    xc = df * xc + one;
    xz = df * xz;
    x = sqr(x) + s * x + c;
  }
  const C F1 = x - z, F2 = xz - mu, a = xz - one, b = xc, d = xzz, e = xzc, det = a * e - b * d;
  const C dz = cdiv(F1 * e - b * F2, det), dc = cdiv(a * F2 - d * F1, det);
  z -= dz; c -= dc;
  dcdmu = cdiv(a, det);
  const Cd dzd(double(dz.r), double(dz.i)), dcd(double(dc.r), double(dc.i));
  return cabs(dzd) + cabs(dcd);
}

// Newton in double to convergence; false if it diverges or stalls above accept
__host__ __device__ static bool solve(const int n, const int shift, const Cd mu, Cd& z, Cd& c, Cd& dcdmu,
                                      const double accept) {
  double last = INFINITY;
  for (int it = 0; it < 40; it++) {
    last = newton_step(n, shift, mu, z, c, dcdmu);
    if (!(last < 1)) return false;
    if (last < 1e-15 * (1 + cabs(c))) return true;
  }
  return last < accept;
}

// Continue radially from μ = 0 (z at the critical point, c at the center) to mu
__host__ __device__ static bool radial(const int n, const int shift, const Cd mu, Cd& z, Cd& c, Cd& dcdmu,
                                       const int steps, const double accept) {
  for (int s = 1; s <= steps; s++)
    if (!solve(n, shift, scale(double(s) / steps, mu), z, c, dcdmu, accept)) return false;
  return true;
}

__host__ __device__ static BulbResult bulb_one(const BulbJob& job, const Ce* mu_e, const E2 pi, const BulbParams p) {
  BulbResult r;
  r.status = bulb_parent;
  r.conv = 0; r.w = 0;
  const int shift = job.shift, n = job.q * job.P;
  const Cd crit(-0.5 * shift, 0.0), l0 = to_d(job.lam0);
  // Parent root and c_W'(λ0)
  Cd cr, dW;
  if (shift && job.P == 1) {  // Cardioid: c = λ/2 - λ²/4, δ = c - 1/4
    cr = scale(0.5, l0) - scale(0.25, l0 * l0) - Cd(0.25, 0);
    dW = scale(0.5, Cd(1, 0) - l0);
  } else {
    Cd zr = crit;
    cr = job.center;
    if (!radial(job.P, shift, l0, zr, cr, dW, 128, p.accept)) return r;
  }
  // Child center: Newton on f^n(crit) = crit
  r.status = bulb_center;
  const double qq = double(job.q) * job.q;
  Cd c = cr + l0 * Cd(dW.r / qq, dW.i / qq);
  bool conv = false;
  for (int it = 0; it < 200 && !conv; it++) {
    Cd x = crit, dx(0);
    for (int k = 0; k < n; k++) { dx = (twice(x) + Cd(shift)) * dx + Cd(1); x = sqr(x) + Cd(shift) * x + c; }
    const Cd step = cdiv(x - crit, dx);
    c -= step;
    const double st = cabs(step);
    conv = st < 1e-16 * (1 + cabs(c));
    if (!(st < 1)) break;
  }
  if (!conv) return r;
  r.status = bulb_period;
  {
    Cd x = crit;
    for (int k = 1; k < n; k++) {
      x = sqr(x) + Cd(shift) * x + c;
      if (n % k == 0 && cabs(x - crit) < 1e-8) return r;
    }
  }
  // Polish the center in E
  Ce ce = to_e(c);
  {
    const Ce crit_e = to_e(crit), s(shift), one(1);
    for (int it = 0; it < 3; it++) {
      Ce x = crit_e, dx(0);
      for (int k = 0; k < n; k++) { dx = (twice(x) + s) * dx + one; x = sqr(x) + s * x + ce; }
      ce -= cdiv(x - crit_e, dx);
    }
  }
  // Boundary: continuation in double, E polish at each point, Green's sum in E
  r.status = bulb_area;
  const int N = p.N;
  Cd z = crit, cb = to_d(ce), dc;
  if (!radial(n, shift, to_d(mu_e[0]), z, cb, dc, p.radial_steps, p.accept)) return r;
  E2 sum(0.0), half_sum(0.0);
  for (int j = 0; j < N; j++) {
    const Cd mu = to_d(mu_e[j]);
    if (j) {
      const Cd prev = to_d(mu_e[j - 1]);
      // Substeps along the chord from prev to mu, projected to the circle
      for (int s = 1; s <= p.substeps; s++) {
        const double t = double(s) / p.substeps;
        Cd m = prev + scale(t, mu - prev);
        m = scale(1 / cabs(m), m);
        if (!solve(n, shift, m, z, cb, dc, p.accept)) return r;
      }
    } else {
      if (!solve(n, shift, mu, z, cb, dc, p.accept)) return r;
    }
    Ce ze = to_e(z), cee = to_e(cb), dce;
    for (int it = 0; it < p.polish; it++) newton_step(n, shift, mu_e[j], ze, cee, dce);
    const E2 t = (conj(cee - ce) * mu_e[j] * dce).r;
    sum = sum + t;
    if (j % 2 == 0) half_sum = half_sum + t;
  }
  const E2 A = pi * sum / E2(int64_t(N)), A2 = pi * half_sum / E2(int64_t(N / 2));
  r.conv = double((A - A2) / A);
  // |c_W'(λ0)|²: exact for the cardioid, else from the double continuation
  E2 w(dW.r * dW.r + dW.i * dW.i);
  if (shift && job.P == 1) {
    const Ce d = Ce(1) - job.lam0;
    w = (sqr(d.r) + sqr(d.i)) / E2(int64_t(4));
  }
  const E2 q4(int64_t(job.q) * job.q * job.q * job.q);
  r.area = A;
  r.F = A * q4 / (pi * w);
  r.w = double(w);
  r.center = to_d(ce) + Cd(0.25 * shift, 0);
  r.status = bulb_ok;
  return r;
}

#ifdef __CUDACC__
__global__ static void bulb_kernel(const int n, const BulbJob* jobs, const Ce* mu_e, const E2 pi, const BulbParams p,
                                   BulbResult* out) {
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += blockDim.x * gridDim.x)
    out[i] = bulb_one(jobs[i], mu_e, pi, p);
}
#endif

}  // namespace

BulbJob bulb_job(const int P, const Complex<double> center, const int p, const int q) {
  static std::map<std::pair<int, int>, Ce> twiddles;
  BulbJob j;
  j.P = P; j.p = p; j.q = q;
  j.shift = P == 1;
  j.center = j.shift ? center - Cd(0.25, 0) : center;
  const int g = std::gcd(p, q);
  const auto key = std::make_pair(p / g, q / g);
  auto it = twiddles.find(key);
  if (it == twiddles.end()) it = twiddles.emplace(key, nearest_twiddle<E2>(key.first, key.second)).first;
  j.lam0 = it->second;
  return j;
}

vector<BulbResult> bulb_areas(const vector<BulbJob>& jobs, const BulbParams& params) {
  const int N = params.N;
  slow_assert(N >= 2 && N % 2 == 0, "bulb_areas: N must be even");
  vector<Ce> mu(N);
  for (int j = 0; j < N; j++) mu[j] = nearest_twiddle<E2>(2 * j + 1, 2 * N);
  const E2 pi = nearest_pi<E2>();
  // Sort by period so that neighboring threads do similar work
  vector<int64_t> order(jobs.size());
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](int64_t a, int64_t b) {
    return int64_t(jobs[a].q) * jobs[a].P < int64_t(jobs[b].q) * jobs[b].P; });
  vector<BulbJob> sorted(jobs.size());
  for (size_t i = 0; i < order.size(); i++) sorted[i] = jobs[order[i]];
  vector<BulbResult> out_sorted(jobs.size());
  const int64_t n = jobs.size();
  if (params.cuda) {
#ifdef __CUDACC__
    Mem<BulbJob> dj(n, true);
    Mem<Ce> dmu(N, true);
    Mem<BulbResult> dout(n, true);
    dj.from_host(sorted.data(), n);
    dmu.from_host(mu.data(), N);
    bulb_kernel<<<32 * num_sms(), 128, 0, stream()>>>(int(n), dj.p, dmu.p, pi, params, dout.p);
    cuda_check(cudaGetLastError());
    dout.to_host(out_sorted.data(), n);
#else
    die("bulb_areas: no CUDA");
#endif
  } else {
    const char* te = getenv("MANDELBROT_THREADS");
    const int T = te ? atoi(te) : int(std::thread::hardware_concurrency());
    std::atomic<int64_t> next(0);
    vector<std::thread> pool;
    for (int t = 0; t < T; t++)
      pool.emplace_back([&]() {
        for (int64_t i; (i = next++) < n;) out_sorted[i] = bulb_one(sorted[i], mu.data(), pi, params);
      });
    for (auto& th : pool) th.join();
  }
  vector<BulbResult> out(n);
  for (int64_t i = 0; i < n; i++) out[order[i]] = out_sorted[i];
  return out;
}

}  // namespace mandelbrot
