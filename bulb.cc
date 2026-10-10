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

// The parameter: c itself, or (local) c = hi + lo + L Δ for deep components, where (hi, lo) is a double-double base
// point and L a length scale, so that the tracked Δ is O(1) on the component and the double phase keeps relative
// precision.  Non-local jobs take exactly the arithmetic they always did.
struct Par {
  bool local;
  Cd hi, lo;
  double L;
};
__host__ __device__ static inline Cd add_c(const Cd v, const Par& P, const Cd c) {
  return P.local ? (v + P.hi) + (P.lo + scale(P.L, c)) : v + c;
}
__host__ __device__ static inline Ce add_c(const Ce v, const Par& P, const Ce c) {
  return P.local ? v + ((to_e(P.hi) + to_e(P.lo)) + E2(P.L) * c) : v + c;
}
__host__ __device__ static inline Cd dc_one(const Par& P) { return P.local ? Cd(P.L) : Cd(1); }
template<class S> __host__ __device__ static inline Complex<S> dc_one_t(const Par& P);
template<> __host__ __device__ inline Cd dc_one_t<double>(const Par& P) { return dc_one(P); }
template<> __host__ __device__ inline Ce dc_one_t<E2>(const Par& P) { return to_e(dc_one(P)); }
__host__ __device__ static inline double cabs(const Ce z) { return cabs(to_d(z)); }
__host__ __device__ static inline Ce scale(const double a, const Ce z) { return Ce(E2(a) * z.r, E2(a) * z.i); }

// One Newton step for (z, c): f^n(z) = z, (f^n)'(z) = mu with f(x) = x² + s x + c.  Returns |dz| + |dc| (double
// norms) and sets dzdmu = z'(μ), dcdmu = c'(μ) (at the input point: J (z', c') = (0, 1) gives c' = a/det,
// z' = -b/det).
template<class S> __host__ __device__ static double newton_step(const int n, const int shift, const Par& P,
                                                                const Complex<S> mu, Complex<S>& z, Complex<S>& c,
                                                                Complex<S>& dzdmu, Complex<S>& dcdmu) {
  typedef Complex<S> C;
  const C s(shift), one(1), lc = dc_one_t<S>(P);
  C x = z, xz(1), xc(0), xzz(0), xzc(0);
  for (int i = 0; i < n; i++) {
    const C df = twice(x) + s;
    const C nxzz = twice(xz * xz) + df * xzz, nxzc = twice(xc * xz) + df * xzc;
    xzz = nxzz; xzc = nxzc;
    xc = df * xc + lc;
    xz = df * xz;
    x = add_c(sqr(x) + s * x, P, c);
  }
  const C F1 = x - z, F2 = xz - mu, a = xz - one, b = xc, d = xzz, e = xzc, det = a * e - b * d;
  const C dz = cdiv(F1 * e - b * F2, det), dc = cdiv(a * F2 - d * F1, det);
  z -= dz; c -= dc;
  dcdmu = cdiv(a, det);
  dzdmu = -cdiv(b, det);
  const Cd dzd(double(dz.r), double(dz.i)), dcd(double(dc.r), double(dc.i));
  return cabs(dzd) + cabs(dcd);
}

// Newton in double to convergence (at most iters steps); false if it diverges or stalls above accept
// (Expansion<2> instances, for local jobs, are done at 1e-28 instead of 1e-15.)
template<class S> __host__ __device__ static bool solve(const int n, const int shift, const Par& P,
                                                        const Complex<S> mu, Complex<S>& z, Complex<S>& c,
                                                        Complex<S>& dzdmu, Complex<S>& dcdmu, const double accept,
                                                        const int iters = 40) {
  const double done = sizeof(S) == sizeof(double) ? 1e-15 : BULB_E == 2 ? 1e-28 : 1e-44;
  double last = INFINITY;
  for (int it = 0; it < iters; it++) {
    last = newton_step(n, shift, P, mu, z, c, dzdmu, dcdmu);
    if (!(last < 1)) return false;
    if (last < done * (1 + cabs(c))) return true;
  }
  return last < accept;
}

// Jacobian of (f^n(z) - z, (f^n)'(z) - mu) in double at (z, c): a = ∂/∂z, b = ∂/∂c of the first, d, e of the second
__host__ __device__ static void jacobian(const int n, const int shift, const Par& P, const Cd z, const Cd c, Cd& a,
                                         Cd& b, Cd& d, Cd& e) {
  const Cd s(shift), one(1), lc = dc_one(P);
  Cd x = z, xz(1), xc(0), xzz(0), xzc(0);
  for (int i = 0; i < n; i++) {
    const Cd df = twice(x) + s;
    const Cd nxzz = twice(xz * xz) + df * xzz, nxzc = twice(xc * xz) + df * xzc;
    xzz = nxzz; xzc = nxzc;
    xc = df * xc + lc;
    xz = df * xz;
    x = add_c(sqr(x) + s * x, P, c);
  }
  a = xz - one; b = xc; d = xzz; e = xzc;
}

// Simplified Newton step in Expansion<2>: the residual in E, the Jacobian (a, b, d, e) from double.  The Jacobian's
// relative error ~1e-16 makes one step take a 1e-15 point to ~1e-30, with only f and f' carried in E.
__host__ __device__ static double polish_e(const int n, const int shift, const Par& P, const Ce mu, Ce& z, Ce& c,
                                           const Cd a, const Cd b, const Cd d, const Cd e) {
  const Ce s(shift);
  Ce x = z, xz(1);
  for (int i = 0; i < n; i++) {
    xz = (twice(x) + s) * xz;
    x = add_c(sqr(x) + s * x, P, c);
  }
  const Ce F1 = x - z, F2 = xz - mu;
  const Cd det = a * e - b * d;
  const Ce ae = to_e(cdiv(a, det)), be = to_e(cdiv(b, det)), de = to_e(cdiv(d, det)), ee = to_e(cdiv(e, det));
  const Ce dz = F1 * ee - be * F2, dc = ae * F2 - de * F1;
  z -= dz; c -= dc;
  return cabs(to_d(dz)) + cabs(to_d(dc));
}

// Predictor-corrector continuation of (z, c) from mu0 to mu1: an Euler predictor along (z', c'), then at most 8
// Newton steps; on failure, split the step (along the chord, projected to the circle if both ends are on it).
// (z', c') are kept current at the final point.
template<class S> __host__ __device__ static bool track(const int n, const int shift, const Par& P,
                                                        const Complex<S> mu0, const Complex<S> mu1, Complex<S>& z,
                                                        Complex<S>& c, Complex<S>& dz, Complex<S>& dc,
                                                        const double accept) {
  typedef Complex<S> Cd;
  struct Seg { Cd a, b; };
  Seg stack[24];
  int top = 0;
  stack[top++] = {mu0, mu1};
  const bool circle = fabs(cabs(mu0) - 1) < 1e-12 && fabs(cabs(mu1) - 1) < 1e-12;
  int work = 0;
  while (top) {
    const Seg sg = stack[--top];
    const Cd d = sg.b - sg.a;
    const Cd z0 = z, c0 = c, dz0 = dz, dc0 = dc;
    z = z + dz * d; c = c + dc * d;
    if (solve(n, shift, P, sg.b, z, c, dz, dc, accept, 8)) continue;
    z = z0; c = c0; dz = dz0; dc = dc0;
    if (top + 2 > 24 || ++work > 4096) return false;
    Cd mid = scale(0.5, sg.a + sg.b);
    if (circle) mid = scale(1 / cabs(mid), mid);
    stack[top++] = {mid, sg.b};
    stack[top++] = {sg.a, mid};
  }
  return true;
}

// Continue from μ = 0 (z at the critical point, c at the center) to mu, in `steps` initial pieces
template<class S> __host__ __device__ static bool radial(const int n, const int shift, const Par& P,
                                                         const Complex<S> mu, Complex<S>& z, Complex<S>& c,
                                                         Complex<S>& dz, Complex<S>& dc, const int steps,
                                                         const double accept) {
  typedef Complex<S> Cd;
  if (!solve(n, shift, P, Cd(0), z, c, dz, dc, accept)) return false;  // Derivatives at the center
  for (int s = 1; s <= steps; s++)
    if (!track(n, shift, P, scale(double(s - 1) / steps, mu), scale(double(s) / steps, mu), z, c, dz, dc, accept))
      return false;
  return true;
}

constexpr int kMaxN = 128;

// tw: the 2N twiddles e^{iπm/N}, m = 0..2N-1 (boundary points μ_j = tw_{2j+1})
__host__ __device__ static BulbResult bulb_one(const BulbJob& job, const Ce* tw, const E2 pi, const BulbParams p) {
  BulbResult r;
  r.status = bulb_parent;
  r.conv = 0; r.w = 0;
  const int shift = job.shift, n = job.P ? job.q * job.P : job.q;  // P = 0: the component itself, period q
  const Cd crit(-0.5 * shift, 0.0), l0 = to_d(job.lam0);
  Par par{false, Cd(0), Cd(0), 1.0};
  // Parent root and c_W'(λ0) (none for P = 0: job.center is the component's own center, to be refined)
  Cd cr, dW(1);
  if (job.P == 0) {
    cr = job.center;
  } else if (shift && job.P == 1) {  // Cardioid: c = λ/2 - λ²/4, δ = c - 1/4
    cr = scale(0.5, l0) - scale(0.25, l0 * l0) - Cd(0.25, 0);
    dW = scale(0.5, Cd(1, 0) - l0);
  } else {
    Cd zr = crit, dzr;
    cr = job.center;
    if (!radial(job.P, shift, par, l0, zr, cr, dzr, dW, 16, p.accept)) return r;
  }
  // Child center: Newton on f^n(crit) = crit
  r.status = bulb_center;
  const double qq = double(job.q) * job.q;
  Cd c = job.P ? cr + l0 * Cd(dW.r / qq, dW.i / qq) : cr;
  bool conv = false;
  if (job.local) {
    // Deep component: base point = the given double-double center, L = the size estimate |1/(β Λ²)| from its orbit
    par = Par{true, job.center, job.center_lo, 1.0};
    Cd x = crit, prod(1), beta(0);
    for (int k = 1; k < n; k++) { x = add_c(sqr(x), par, Cd(0)); prod = prod * twice(x); beta = beta + cdiv(Cd(1), prod); }
    par.L = 1 / (cabs(beta) * cabs(prod) * cabs(prod));
    c = Cd(0);
    conv = true;
  }
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
      x = add_c(sqr(x) + Cd(shift) * x, par, c);
      if (n % k == 0 && cabs(x - crit) < 1e-8) return r;
    }
  }
  // Polish the center in E
  Ce ce = to_e(c);
  {
    const Ce crit_e = to_e(crit), s(shift), one = to_e(dc_one(par));
    for (int it = 0; it < (par.local ? 6 : 3); it++) {
      Ce x = crit_e, dx(0);
      for (int k = 0; k < n; k++) { dx = (twice(x) + s) * dx + one; x = add_c(sqr(x) + s * x, par, ce); }
      ce -= cdiv(x - crit_e, dx);
    }
  }
  // Boundary: continuation in double, simplified-Newton E polish at each point μ_j = tw_{2j+1}, then the area from
  // the Taylor coefficients of c(μ) - c0 = Σ a_k μ^k: area = π Σ k |a_k|², with a_k by DFT of the boundary values
  r.status = bulb_area;
  const int N = p.N;
  Ce cj[kMaxN];
  if (par.local) {
    // Deep components: the whole continuation in Expansion<2> (a double periodic point would carry the error the
    // local coordinates remove from c: ~1e-16 n |Λ|² relative to the component)
    Ce z = to_e(crit), cb = ce, dz, dc;
    if (!radial(n, shift, par, tw[1], z, cb, dz, dc, p.radial_steps, p.accept)) return r;
    for (int j = 0; j < N; j++) {
      if (j) {
        const Ce mup = tw[2 * j - 1], muj = tw[2 * j + 1];
        for (int s = 1; s <= p.substeps; s++) {
          Ce a = mup + scale(double(s - 1) / p.substeps, muj - mup);
          Ce b = mup + scale(double(s) / p.substeps, muj - mup);
          a = scale(1 / cabs(a), a); b = scale(1 / cabs(b), b);
          if (!track(n, shift, par, a, b, z, cb, dz, dc, p.accept)) return r;
        }
      }
      cj[j] = cb - ce;
    }
  }
  Cd z = crit, cb = to_d(ce), dz, dc;
  if (!par.local && !radial(n, shift, par, to_d(tw[1]), z, cb, dz, dc, p.radial_steps, p.accept)) return r;
  for (int j = 0; !par.local && j < N; j++) {
    const Cd muj = to_d(tw[2 * j + 1]);
    if (j) {
      // Along the circle from the previous point, in `substeps` initial pieces
      const Cd mup = to_d(tw[2 * j - 1]);
      for (int s = 1; s <= p.substeps; s++) {
        Cd a = mup + scale(double(s - 1) / p.substeps, muj - mup);
        Cd b = mup + scale(double(s) / p.substeps, muj - mup);
        a = scale(1 / cabs(a), a); b = scale(1 / cabs(b), b);
        if (!track(n, shift, par, a, b, z, cb, dz, dc, p.accept)) return r;
      }
    }
    // One simplified-Newton step, and more (up to p.polish) only while it still moves (near-parabolic bulbs, where
    // the double solution stalls short of full precision)
    Cd ja, jb, jd, je;
    jacobian(n, shift, par, z, cb, ja, jb, jd, je);
    Ce ze = to_e(z), cee = to_e(cb);
    for (int it = 0; it < p.polish; it++)
      if (polish_e(n, shift, par, tw[2 * j + 1], ze, cee, ja, jb, jd, je) < (BULB_E == 2 ? 1e-24 : 1e-40) * (1 + cabs(cb))) break;
    cj[j] = cee - ce;
  }
  // a_k = (1/N) Σ_j c_j μ_j^-k, μ_j^-k = conj(tw_{k(2j+1) mod 2N}); the N/2 subrule uses even j only
  E2 sum(0.0), half_sum(0.0);
  for (int k = 1; k < N; k++) {
    Ce ak(0), hk(0);
    for (int j = 0; j < N; j++) {
      const Ce t = cj[j] * conj(tw[(int64_t(k) * (2 * j + 1)) % (2 * N)]);
      ak += t;
      if (j % 2 == 0 && k < N / 2) hk += t;
    }
    sum = sum + E2(int64_t(k)) * (sqr(ak.r) + sqr(ak.i));
    if (k < N / 2) half_sum = half_sum + E2(int64_t(k)) * (sqr(hk.r) + sqr(hk.i));
  }
  const E2 NN(int64_t(N) * N), hN(int64_t(N / 2) * (N / 2));
  const E2 L2 = par.local ? E2(par.L) * E2(par.L) : E2(1.0);  // Local jobs: c - c0 = L (Δ - Δ0)
  const E2 A = L2 * (pi * sum / NN), A2 = L2 * (pi * half_sum / hN);
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
  r.center = par.local ? job.center + (job.center_lo + scale(par.L, to_d(ce))) : to_d(ce) + Cd(0.25 * shift, 0);
  r.status = bulb_ok;
  return r;
}

#ifdef __CUDACC__
__global__ static void bulb_kernel(const int n, const BulbJob* jobs, const Ce* tw, const E2 pi, const BulbParams p,
                                   BulbResult* out) {
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += blockDim.x * gridDim.x)
    out[i] = bulb_one(jobs[i], tw, pi, p);
}
#endif

}  // namespace

BulbJob bulb_job(const int P, const Complex<double> center, const int p, const int q) {
  static std::map<std::pair<int, int>, Ce> twiddles;
  BulbJob j;
  j.P = P; j.p = p; j.q = q;
  j.shift = P == 1;
  j.local = false;
  j.center_lo = Cd(0);
  if (P == 0) {  // The component of period q with center near `center` (any type, e.g. primitive)
    j.center = center;
    j.lam0 = Ce(1);
    return j;
  }
  j.center = j.shift ? center - Cd(0.25, 0) : center;
  const int g = std::gcd(p, q);
  const auto key = std::make_pair(p / g, q / g);
  auto it = twiddles.find(key);
  if (it == twiddles.end()) it = twiddles.emplace(key, nearest_twiddle<E2>(key.first, key.second)).first;
  j.lam0 = it->second;
  return j;
}

BulbJob bulb_job_local(const Complex<double> center, const Complex<double> center_lo, const int q) {
  BulbJob j = bulb_job(0, center, 1, q);
  j.local = true;
  j.center_lo = center_lo;
  return j;
}

vector<BulbResult> bulb_areas(const vector<BulbJob>& jobs, const BulbParams& params) {
  const int N = params.N;
  slow_assert(N >= 4 && N % 4 == 0 && N <= kMaxN, "bulb_areas: need 4 | N ≤ %d", kMaxN);
  vector<Ce> mu(2 * N);
  for (int m = 0; m < 2 * N; m++) mu[m] = nearest_twiddle<E2>(m, 2 * N);
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
    Mem<Ce> dmu(2 * N, true);
    Mem<BulbResult> dout(n, true);
    dj.from_host(sorted.data(), n);
    dmu.from_host(mu.data(), 2 * N);
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
