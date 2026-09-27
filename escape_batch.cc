// Batched escape-time sampling of leaf cells, on CPU threads or the GPU

#include "escape_batch.h"
#include "debug.h"
#include <algorithm>
#include <atomic>
#include <thread>
#include <vector>
#ifdef __CUDACC__
#include "array.h"
#endif
#ifndef BATCH_LANES
#define BATCH_LANES 4   // Orbits in flight per CPU thread
#endif
#ifndef BATCH_BURST
#define BATCH_BURST 16  // Iterations per orbit per turn
#endif
namespace mandelbrot {

using std::span;

namespace {

// Each CPU thread keeps L orbits in flight, advancing each by short bursts in turn, so that out-of-order
// execution overlaps their independent dependency chains
template<class T, int L> int64_t cpu_worker(span<const Leaf> leaves, const SampleParams& p, span<uint32_t> bits,
                                            std::atomic<int64_t>& next) {
  const int64_t n = int64_t(leaves.size()) * p.m, chunk = 1024, burst = BATCH_BURST;
  Orbit<T> o[L];
  int64_t idx[L], lo = 0, hi = 0, iters = 0;
  // Give lane l its next undecided sample, writing any samples decided at start.  Returns false when out.
  const auto load = [&](const int l) {
    for (;;) {
      if (lo == hi) {
        lo = next.fetch_add(chunk);
        hi = std::min(lo + chunk, n);
        if (lo >= n) { lo = hi = n; idx[l] = -1; return false; }
      }
      const int64_t i = lo++;
      double x, y;
      sample_point(leaves[i / p.m], p, i, x, y);
      if (!o[l].start(x, y)) { idx[l] = i; return true; }
      bits[i] = below_bits(o[l].e, p);
    }
  };
  int active = 0;
  for (int l = 0; l < L; l++) active += load(l);
  while (active) {
    for (int l = 0; l < L; l++) {
      if (idx[l] < 0 || !o[l].run(p.max_iter, burst)) continue;
      bits[idx[l]] = below_bits(o[l].e, p);
      iters += o[l].e.iters;
      if (!load(l)) active--;
    }
  }
  return iters;
}

}  // namespace

template<class T> int64_t sample_leaves_cpu(span<const Leaf> leaves, const SampleParams& p, span<uint32_t> bits) {
  slow_assert(bits.size() == leaves.size() * size_t(p.m));
  slow_assert(0 < p.K && p.K <= 32);
  std::atomic<int64_t> next(0), iters(0);
  std::vector<std::thread> pool;
  for (int t = 0; t < int(std::thread::hardware_concurrency()); t++)
    pool.emplace_back([&]() { iters += cpu_worker<T, BATCH_LANES>(leaves, p, bits, next); });
  for (auto& t : pool) t.join();
  return iters;
}

#ifdef __CUDACC__

namespace {

// Persistent threads: each thread claims chunks of samples from a global counter and refills its orbit as
// soon as it finishes, so a slow orbit only delays its own thread, not its warp's next samples
template<class T> __global__ void sample_kernel(const Leaf* leaves, const SampleParams p, const int64_t n,
                                                uint32_t* bits, unsigned long long* next,
                                                unsigned long long* iters) {
  const int64_t chunk = 16, burst = 16;
  Orbit<T> o;
  int64_t i = -1, end = -1;
  unsigned long long it = 0;
  bool done = true;
  for (;;) {
    if (done) {
      if (i >= 0) { bits[i] = below_bits(o.e, p); it += o.e.iters; }
      if (++i >= end) {
        i = int64_t(atomicAdd(next, (unsigned long long)chunk));
        if (i >= n) break;
        end = min(i + chunk, n);
      }
      double x, y;
      sample_point(leaves[i / p.m], p, i, x, y);
      done = o.start(x, y);
      if (done) continue;
    }
    done = o.run(p.max_iter, burst);
  }
  atomicAdd(iters, it);
}

}  // namespace

template<class T> int64_t sample_leaves_cuda(span<const Leaf> leaves, const SampleParams& p, span<uint32_t> bits) {
  slow_assert(bits.size() == leaves.size() * size_t(p.m));
  slow_assert(0 < p.K && p.K <= 32);
  const int64_t n = int64_t(bits.size());
  if (!n) return 0;
  Array<Device<Leaf>> dleaves(leaves.size());
  Array<Device<uint32_t>> dbits(n);
  Array<Device<unsigned long long>> counters(2);
  host_to_device<Leaf>(dleaves, leaves);
  cuda_check(cudaMemsetAsync(device_get(counters), 0, 2 * sizeof(unsigned long long), stream()));
  sample_kernel<T><<<8 * num_sms(), 256, 0, stream()>>>(device_get(dleaves), p, n, device_get(dbits),
                                                        device_get(counters), device_get(counters) + 1);
  cuda_check(cudaGetLastError());
  device_to_host<uint32_t>(bits, dbits);
  unsigned long long h[2];
  device_to_host<unsigned long long>(span<unsigned long long>(h, 2), counters);
  return int64_t(h[1]);
}

#else  // !__CUDACC__

template<class T> int64_t sample_leaves_cuda(span<const Leaf>, const SampleParams&, span<uint32_t>) {
  die("sample_leaves_cuda: built without CUDA");
}

#endif  // __CUDACC__

#define INSTANTIATE(T) \
  template int64_t sample_leaves_cpu<T>(span<const Leaf>, const SampleParams&, span<uint32_t>); \
  template int64_t sample_leaves_cuda<T>(span<const Leaf>, const SampleParams&, span<uint32_t>);
INSTANTIATE(float)
INSTANTIATE(double)

}  // namespace mandelbrot
