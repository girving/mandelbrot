// GPU children for the general-root Lavaurs model (glavaurs children, glavaurs_core.h)
//
// One thread per (job, shift j, start): the start (the two roots of the source's local quadratic, then `starts`
// points on each of three circles) and Newton on Θ_r(σ) = σ_c + j/q exactly as the CPU mode; the host deduplicates
// per (job, j) in start order; a second kernel verifies each new solution as a center (excursion n_c - j) and computes
// its weight w = |Θ_r' Π H'|^-2.
#include "glavaurs_core.h"
#include <algorithm>
#include <cstdio>
#ifdef __CUDACC__
#include <cuda_runtime.h>
#endif
namespace mandelbrot {

#ifndef __CUDACC__

bool glavaurs_gpu_children(const GLCore&, const std::vector<GLChildJob>&, int, std::vector<GLChild>&) { return false; }

#else  // __CUDACC__

namespace {

typedef Complex<double> Cd;
using glcore::cabs;
using glcore::carg;
using glcore::cdiv;
using glcore::cpolar;

__host__ __device__ static inline Cd csqrt(const Cd z) { return cpolar(::sqrt(cabs(z)), carg(z) / 2); }

struct JobDev {
  int r, nc, jlo, nst;
  Cd u, c, t0, d0, dd0;
  long long first;   // flat index of the job's first task
};

__device__ static Cd task_start(const GLCore& core, const JobDev& J, const int j, const int si, const int m) {
  const Cd y = Cd(J.c.r + double(j) / core.q, J.c.i) + core.target_shift(j);
  const Cd a = Cd(0.5) * J.dd0, b = J.d0, cc = J.t0 - y;
  const Cd sq = csqrt(b * b - Cd(4) * a * cc);
  const Cd r0 = J.u + cdiv(-b + sq, Cd(2) * a);
  if (si == 0) return r0;
  if (si == 1) return J.u + cdiv(-b - sq, Cd(2) * a);
  const double rho = cabs(r0 - J.u);
  const int k = si - 2, fi = k / m, kk = k % m;
  const double f = fi == 0 ? 0.5 : fi == 1 ? 1.0 : 2.0;
  return J.u + cpolar(f * rho, 2 * M_PI * (kk + 0.5) / m);
}

__global__ void newton_kernel(const GLCore core, const JobDev* jobs, const int njobs, const int m, const long long base,
                              const long long count, Cd* sol, unsigned char* ok) {
  const long long i = blockIdx.x * (long long)blockDim.x + threadIdx.x;
  if (i >= count) return;
  const long long t = base + i;
  int lo = 0, hi = njobs - 1;
  while (lo < hi) { const int mid = (lo + hi + 1) / 2; if (jobs[mid].first <= t) lo = mid; else hi = mid - 1; }
  const JobDev& J = jobs[lo];
  const long long off = t - J.first;
  const int j = J.jlo + int(off / J.nst), si = int(off % J.nst);
  ok[i] = 0;
  if (J.nc - j < 0) return;
  const Cd y = Cd(J.c.r + double(j) / core.q, J.c.i) + core.target_shift(j);
  Cd s = task_start(core, J, j, si, m);
  for (int it = 0; it < 60; it++) {
    Cd th, d, dd, hp;
    if (!core.theta(s, J.r, th, d, dd, hp)) break;
    const Cd step = cdiv(th - y, d);
    s = s - step;
    if (!(cabs(step) < 10)) break;
    if (cabs(step) < 1e-11 * (1 + cabs(s))) { ok[i] = 1; break; }
  }
  sol[i] = s;
}

struct Cand { int r, n; Cd s; };

__global__ void verify_kernel(const GLCore core, const Cand* cands, const long long count, unsigned char* valid,
                              double* w) {
  const long long i = blockIdx.x * (long long)blockDim.x + threadIdx.x;
  if (i >= count) return;
  const Cand c = cands[i];
  valid[i] = 0;
  Cd wv = core.crit, cen = c.s;
  if (!(core.newton(c.r, c.n, wv, cen, Cd(0), 1e-14) < 1e-10) || cabs(cen - c.s) > 1e-8) return;
  Cd th, d, dd, hp;
  if (!core.theta(c.s, c.r, th, d, dd, hp)) return;
  valid[i] = 1;
  w[i] = 1 / sqr_abs(d * hp);
}

#define GL_CUDA(e) do { const cudaError_t _e = (e); if (_e != cudaSuccess) { \
  fprintf(stderr, "glavaurs gpu: %s: %s\n", #e, cudaGetErrorString(_e)); return false; } } while (0)

}  // namespace

bool glavaurs_gpu_children(const GLCore& core, const std::vector<GLChildJob>& jobs, const int m,
                           std::vector<GLChild>& out) {
  int devices = 0;
  if (cudaGetDeviceCount(&devices) != cudaSuccess || devices == 0) return false;
  // Per-job source data (Θ_r and its derivatives at the source) on the host
  std::vector<JobDev> jd;
  std::vector<int> jix;   // the job each JobDev came from
  const int nst = 2 + 3 * m;
  long long total = 0;
  for (int i = 0; i < int(jobs.size()); i++) {
    const auto& J = jobs[i];
    Cd t0, d0, dd0, hp0;
    if (!core.theta(J.u, J.r, t0, d0, dd0, hp0)) continue;
    if (J.jhi < J.jlo) continue;
    jd.push_back(JobDev{J.r, J.nc, J.jlo, nst, J.u, J.c, t0, d0, dd0, total});
    jix.push_back(i);
    total += (long long)(J.jhi - J.jlo + 1) * nst;
  }
  out.clear();
  if (jd.empty()) return true;
  JobDev* djobs;
  GL_CUDA(cudaMalloc(&djobs, jd.size() * sizeof(JobDev)));
  GL_CUDA(cudaMemcpy(djobs, jd.data(), jd.size() * sizeof(JobDev), cudaMemcpyHostToDevice));
  const long long chunk = std::min(1ll << 24, total);
  Cd* dsol;
  unsigned char* dok;
  GL_CUDA(cudaMalloc(&dsol, chunk * sizeof(Cd)));
  GL_CUDA(cudaMalloc(&dok, chunk));
  std::vector<Cd> sol(total);
  std::vector<unsigned char> ok(total);
  for (long long base = 0; base < total; base += chunk) {
    const long long count = std::min(chunk, total - base);
    newton_kernel<<<(count + 255) / 256, 256>>>(core, djobs, int(jd.size()), m, base, count, dsol, dok);
    GL_CUDA(cudaGetLastError());
    GL_CUDA(cudaMemcpy(sol.data() + base, dsol, count * sizeof(Cd), cudaMemcpyDeviceToHost));
    GL_CUDA(cudaMemcpy(ok.data() + base, dok, count, cudaMemcpyDeviceToHost));
    fprintf(stderr, "glavaurs gpu: newton %lld / %lld tasks\n", base + count, total);
  }
  cudaFree(dsol);
  cudaFree(dok);
  cudaFree(djobs);
  // Deduplicate per (job, j) in start order, as the CPU mode does (converged solutions, verified or not)
  std::vector<Cand> cands;
  std::vector<GLChild> pend;
  for (size_t ji = 0; ji < jd.size(); ji++) {
    const auto& J = jd[ji];
    const auto& G = jobs[jix[ji]];
    for (int j = G.jlo; j <= G.jhi; j++) {
      std::vector<Cd> found;
      const long long b = J.first + (long long)(j - G.jlo) * nst;
      for (int si = 0; si < nst; si++) {
        if (!ok[b + si]) continue;
        const Cd s = sol[b + si];
        bool dup = false;
        for (const auto& f : found) dup |= cabs(f - s) < 1e-9;
        if (dup) continue;
        found.push_back(s);
        cands.push_back(Cand{G.r, G.nc - j, s});
        pend.push_back(GLChild{jix[ji], j, si, s, 0, cabs(s - G.u) < 4 * G.rad});
      }
    }
  }
  fprintf(stderr, "glavaurs gpu: %zu distinct solutions to verify\n", cands.size());
  if (cands.empty()) return true;
  Cand* dc;
  unsigned char* dv;
  double* dw;
  GL_CUDA(cudaMalloc(&dc, cands.size() * sizeof(Cand)));
  GL_CUDA(cudaMalloc(&dv, cands.size()));
  GL_CUDA(cudaMalloc(&dw, cands.size() * sizeof(double)));
  GL_CUDA(cudaMemcpy(dc, cands.data(), cands.size() * sizeof(Cand), cudaMemcpyHostToDevice));
  verify_kernel<<<(cands.size() + 255) / 256, 256>>>(core, dc, (long long)cands.size(), dv, dw);
  GL_CUDA(cudaGetLastError());
  std::vector<unsigned char> valid(cands.size());
  std::vector<double> w(cands.size());
  GL_CUDA(cudaMemcpy(valid.data(), dv, cands.size(), cudaMemcpyDeviceToHost));
  GL_CUDA(cudaMemcpy(w.data(), dw, cands.size() * sizeof(double), cudaMemcpyDeviceToHost));
  cudaFree(dc);
  cudaFree(dv);
  cudaFree(dw);
  for (size_t i = 0; i < pend.size(); i++)
    if (valid[i]) { pend[i].w = w[i]; out.push_back(pend[i]); }
  return true;
}

#endif  // __CUDACC__

}  // namespace mandelbrot
