// Parallel execution of resumable orbits and simple loops, on CPU threads or the GPU
//
// The same task code runs in both places, so CPU runs are exact references for GPU runs.  A task describes
// n independent items, each an orbit (Orbit, OrbitDE, ...) that is started, advanced in short bursts, and
// finished.  Workers claim items in scrambled order (slow orbits cluster spatially, so consecutive claims
// should land far apart) and refill each lane or GPU thread as soon as its orbit finishes.  On the GPU,
// orbits that outlive a step budget are parked in an overflow buffer and finished in a second pass, so one
// slow orbit does not hold up the whole launch while the GPU is full.
//
// Task interface (all __host__ __device__, and the task itself trivially copyable):
//   typedef ... State;                       // Orbit type
//   bool start(State& o, int64_t i) const;   // Start item i; true if already done
//   bool run(State& o) const;                // Advance a burst; true when done
//   void finish(const State& o, int64_t i) const;
//   int64_t iters(const State& o) const;     // Iterations performed, for accounting
//   int64_t progress(const State& o) const;  // Steps so far of a running orbit (for timing only)
//   bool pending(const State& o) const;      // run stopped with work to do together with other lanes
//   bool settle(State& o) const;             // Do that work; true if the orbit is done
//   int64_t burst;                           // Steps per run call
//   int min_blocks;                          // GPU: resident 256-thread blocks per SM (1 to 4) to budget registers for
#pragma once

#include "cutil.h"
#include "debug.h"
#include "noncopyable.h"
#include "print.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <thread>
#include <type_traits>
#include <vector>
namespace mandelbrot {

// splitmix64 finalizer
__host__ __device__ static inline uint64_t mix64(uint64_t z) {
  z += 0x9e3779b97f4a7c15;
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
  z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
  return z ^ (z >> 31);
}

// Uniform double in [0, 1) from (seed, key, j)
__host__ __device__ static inline double uniform(const uint64_t seed, const uint64_t key, const int j) {
  return double(mix64(seed ^ mix64(2 * key + j)) >> 11) * 0x1p-53;
}

// Claim order j → j * stride mod n.  scramble_stride returns a stride near n / φ coprime to n (requires
// n < 2^31, so products fit in 64 bits).
int64_t scramble_stride(int64_t n);
__host__ __device__ static inline int64_t scramble(const int64_t j, const int64_t stride, const int64_t n) {
  return int64_t(uint64_t(j) * uint64_t(stride) % uint64_t(n));
}

// CPU threads to use: $MANDELBROT_THREADS if set (say, a container's CPU request), else all hardware threads
int cpu_threads();

// Integer environment variable, or fallback if unset
int env_int(const char* name, int fallback);

// Relaxed atomic add, on host or device
__host__ __device__ static inline void atomic_add(uint64_t* p, const uint64_t v) {
#ifdef __CUDA_ARCH__
  atomicAdd(reinterpret_cast<unsigned long long*>(p), static_cast<unsigned long long>(v));
#else
  __atomic_fetch_add(p, v, __ATOMIC_RELAXED);
#endif
}
__host__ __device__ static inline uint64_t atomic_fetch_add(uint64_t* p, const uint64_t v) {
#ifdef __CUDA_ARCH__
  return atomicAdd(reinterpret_cast<unsigned long long*>(p), static_cast<unsigned long long>(v));
#else
  return __atomic_fetch_add(p, v, __ATOMIC_RELAXED);
#endif
}

// A buffer of trivially copyable T on the host or the device
template<class T> struct Mem : public Noncopyable {
  T* p = nullptr;
  int64_t n = 0;
  bool cuda = false;

  Mem(const int64_t n, const bool cuda) : n(n), cuda(cuda) {
    static_assert(std::is_trivially_copyable_v<T>);
    if (!n) return;
    if (cuda) {
      IF_CUDA(cuda_check(cudaMallocAsync(&p, n * sizeof(T), stream())));
      CUDA_OR_DIE();
    } else {
      p = static_cast<T*>(malloc(n * sizeof(T)));
      slow_assert(p, "Mem: out of memory for %d elements", n);
    }
  }
  ~Mem() {
    if (!p) return;
    if (cuda) IF_CUDA(cudaFreeAsync(p, stream()));
    else free(p);
  }

  void zero() {
    if (!n) return;
    if (cuda) IF_CUDA(cuda_check(cudaMemsetAsync(p, 0, n * sizeof(T), stream())));
    else memset(p, 0, n * sizeof(T));
  }
  void to_host(T* dst, const int64_t count) const {
    slow_assert(count <= n);
    if (!count) return;
    if (cuda) IF_CUDA(cuda_check(cudaMemcpy(dst, p, count * sizeof(T), cudaMemcpyDeviceToHost)));
    else memcpy(dst, p, count * sizeof(T));
  }
  void from_host(const T* src, const int64_t count) {
    slow_assert(count <= n);
    if (!count) return;
    if (cuda) IF_CUDA(cuda_check(cudaMemcpy(p, src, count * sizeof(T), cudaMemcpyHostToDevice)));
    else memcpy(p, src, count * sizeof(T));
  }
  T get(const int64_t i) const { T x; slow_assert(i < n); if (cuda) IF_CUDA(cuda_check(cudaMemcpy(&x, p + i,
      sizeof(T), cudaMemcpyDeviceToHost))); else x = p[i]; return x; }
};

// Statistics of one run_orbits call
struct RunStats {
  int64_t iters = 0;     // Total iterations
  int64_t overflow = 0;  // Orbits deferred to the second pass (GPU only)
  double secs = 0;
};

namespace engine_detail {

// CPU worker: L orbits in flight, advanced by bursts in turn so that out-of-order execution overlaps them
template<class Task, int L> int64_t cpu_worker(const Task& task, const int64_t n, const int64_t stride,
                                               std::atomic<int64_t>& next) {
  const int64_t chunk = 1024;
  typename Task::State o[L];
  int64_t idx[L], lo = 0, hi = 0, iters = 0;
  // Give lane l its next unfinished item, finishing any items done at start.  Returns false when out.
  const auto load = [&](const int l) {
    for (;;) {
      if (lo == hi) {
        lo = next.fetch_add(chunk);
        hi = std::min(lo + chunk, n);
        if (lo >= n) { lo = hi = n; idx[l] = -1; return false; }
      }
      const int64_t i = scramble(lo++, stride, n);
      if (!task.start(o[l], i)) { idx[l] = i; return true; }
      task.finish(o[l], i);
      iters += task.iters(o[l]);
    }
  };
  int active = 0;
  for (int l = 0; l < L; l++) active += load(l);
  while (active) {
    for (int l = 0; l < L; l++) {
      if (idx[l] < 0 || !task.run(o[l])) continue;
      if (task.pending(o[l]) && !task.settle(o[l])) continue;
      task.finish(o[l], idx[l]);
      iters += task.iters(o[l]);
      if (!load(l)) active--;
    }
  }
  return iters;
}

#ifdef __CUDACC__

// Counters: [0: next claim, 1: iterations, 2: parked, 3: max iterations of one thread, 4-7 with timing: cycles
// inside run, total cycles, active lanes at run calls, warp lanes at run calls, 8: next resume claim,
// 9: iterations in rounds, 10-11 with timing: lane steps and warp steps, 12: parked in this round].

// Persistent threads.  Orbits that are pending or have run `budget` bursts are parked (if room) and finished in
// rounds of settle_kernel and resume_kernel.
template<class Task, bool timing, int min_blocks> __global__ void __launch_bounds__(256, min_blocks)
orbit_kernel(const Task task, const int64_t n, const int64_t stride, const int32_t budget,
             typename Task::State* overflow, int64_t* overflow_items, const int64_t overflow_cap,
             unsigned long long* counters, const int32_t chunk) {
  // Lanes stay in the loop until their whole warp is out of work, and reconverge before each burst, so that
  // lanes refilling at different times do not split the warp into groups that each step half empty.
  // Item indices fit in 32 bits (n < 2^31), which saves registers.
  typename Task::State o;
  int32_t end = 0, i = -1, bursts = 0;  // i = current item or -1
  int64_t j = 0, pos = 0;  // pos = scramble(j - 1).  64 bits: pos + stride can exceed 2^31.
  unsigned long long it = 0, run_cycles = 0, active = 0, slots = 0, lane_steps = 0, warp_steps = 0;
  const long long t0 = timing ? clock64() : 0;
  bool done = true, out = false;
  for (;;) {
    // Finish and refill until this lane has a running orbit or runs out of work
    while (done && !out) {
      if (i >= 0) { task.finish(o, i); it += task.iters(o); }
      if (j == end) {
        j = int64_t(atomicAdd(counters, (unsigned long long)chunk));
        if (j >= n) { out = true; i = -1; break; }
        end = int32_t(min(j + chunk, n));
        pos = scramble(j, stride, n);
      } else {
        pos += stride;  // scramble(j) from scramble(j - 1), without a 64-bit modulus
        if (pos >= n) pos -= n;
      }
      i = int32_t(pos);
      j++;
      bursts = 0;
      done = task.start(o, i);
    }
    if (__all_sync(0xffffffff, out)) break;
    __syncwarp();
    if (!out) {
      if constexpr (timing) {
        const unsigned mask = __activemask();
        const bool leader = (threadIdx.x & 31) == __ffs(mask) - 1;
        if (leader) { active += __popc(mask); slots += 32; }
        const int64_t p0 = task.progress(o);
        const long long r0 = clock64();
        done = task.run(o);
        run_cycles += clock64() - r0;
        // Steps this lane took, against the warp's longest: idle lanes within bursts
        const unsigned steps = unsigned(task.progress(o) - p0), longest = __reduce_max_sync(mask, steps);
        lane_steps += steps;
        if (leader) warp_steps += 32ull * longest;
      } else {
        done = task.run(o);
      }
      // Park orbits that are pending (so that lanes settle together later) or have used their step budget.
      // Parked orbits' iterations are counted when they finish.
      const bool pending = done && task.pending(o);
      if (pending || (!done && ++bursts == budget)) {
        const int64_t k = int64_t(atomicAdd(counters + 2, 1ull));
        if (k < overflow_cap) {
          overflow[k] = o;
          overflow_items[k] = i;
          i = -1;
          done = true;
        } else if (pending) {
          done = task.settle(o);  // No room: settle alone
        }
      }
    }
  }
  atomicAdd(counters + 1, it);
  atomicMax(counters + 3, it);
  if constexpr (timing) {
    atomicAdd(counters + 4, run_cycles);
    atomicAdd(counters + 5, (unsigned long long)(clock64() - t0));
    atomicAdd(counters + 6, active);
    atomicAdd(counters + 7, slots);
    atomicAdd(counters + 10, lane_steps);
    atomicAdd(counters + 11, warp_steps);
  }
}

// Launch orbit_kernel with register pressure chosen at run time: min_blocks resident 256-thread blocks per SM
template<class Task, bool timing> void launch_orbit_kernel(const int min_blocks, const int blocks, const Task& task,
                                                           const int64_t n, const int64_t stride, const int32_t budget,
                                                           typename Task::State* overflow, int64_t* items,
                                                           const int64_t cap, unsigned long long* counters,
                                                           const int32_t chunk) {
#define LAUNCH(b) orbit_kernel<Task, timing, b><<<blocks, 256, 0, stream()>>>(task, n, stride, budget, overflow, \
                                                                             items, cap, counters, chunk)
  switch (min_blocks) {
    case 1: LAUNCH(1); break;
    case 2: LAUNCH(2); break;
    case 3: LAUNCH(3); break;
    case 4: LAUNCH(4); break;
    default: die("MANDELBROT_CUDA_MIN_BLOCKS must be 1 to 4, got %d", min_blocks);
  }
#undef LAUNCH
}

// Settle all pending parked orbits at once, finishing those that are done (marked by item -1)
template<class Task> __global__ void settle_kernel(const Task task, typename Task::State* parked, int64_t* items,
                                                   const int64_t count, unsigned long long* counters) {
  unsigned long long it = 0;
  for (int64_t k = int64_t(blockIdx.x) * blockDim.x + threadIdx.x; k < count; k += int64_t(blockDim.x) * gridDim.x) {
    typename Task::State o = parked[k];
    if (!task.pending(o)) continue;
    if (task.settle(o)) {
      task.finish(o, items[k]);
      it += task.iters(o);
      items[k] = -1;
    } else {
      parked[k] = o;
    }
  }
  atomicAdd(counters + 1, it);
  atomicAdd(counters + 9, it);
}

// Resume parked orbits (skipping finished ones) with persistent warp-synchronous threads like orbit_kernel's,
// until done or pending again, when they park into next
template<class Task> __global__ void resume_kernel(const Task task, const typename Task::State* parked,
                                                   const int64_t* items, const int64_t count,
                                                   typename Task::State* next, int64_t* next_items, const int64_t cap,
                                                   unsigned long long* counters) {
  typename Task::State o;
  int64_t i = -1;
  unsigned long long it = 0;
  bool done = true, out = false;
  for (;;) {
    while (done && !out) {
      if (i >= 0) { task.finish(o, i); it += task.iters(o); }
      i = -1;
      const int64_t k = int64_t(atomicAdd(counters + 8, 1ull));
      if (k >= count) { out = true; break; }
      if (items[k] < 0) continue;
      o = parked[k];
      i = items[k];
      done = false;
    }
    if (__all_sync(0xffffffff, out)) break;
    __syncwarp();
    if (!out) {
      done = task.run(o);
      if (done && task.pending(o)) {
        const int64_t k = int64_t(atomicAdd(counters + 12, 1ull));
        if (k < cap) {
          next[k] = o;
          next_items[k] = i;
          i = -1;
        } else {
          done = task.settle(o);  // No room: settle alone
        }
      }
    }
  }
  atomicAdd(counters + 1, it);
  atomicAdd(counters + 9, it);
}

template<class F> __global__ void for_each_kernel(const int64_t n, const F f) {
  for (int64_t i = int64_t(blockIdx.x) * blockDim.x + threadIdx.x; i < n; i += int64_t(blockDim.x) * gridDim.x)
    f(i);
}

#endif  // __CUDACC__

}  // namespace engine_detail

// Run all n items of task
template<class Task> RunStats run_orbits(const Task& task, const int64_t n, const bool cuda) {
  RunStats stats;
  if (!n) return stats;
  const auto t0 = std::chrono::steady_clock::now();
  const int64_t stride = scramble_stride(n);
  if (!cuda) {
    std::atomic<int64_t> next(0), iters(0);
    std::vector<std::thread> pool;
    for (int t = 0; t < cpu_threads(); t++)
      pool.emplace_back([&]() { iters += engine_detail::cpu_worker<Task, 4>(task, n, stride, next); });
    for (auto& t : pool) t.join();
    stats.iters = iters;
  } else {
#ifdef __CUDACC__
    typedef typename Task::State O;
    static const int blocks_per_sm = env_int("MANDELBROT_CUDA_BLOCKS_PER_SM", 8),
                     min_blocks_env = env_int("MANDELBROT_CUDA_MIN_BLOCKS", 0),  // Override task.min_blocks
                     block = 256,
                     // Steps before parking.  At the first leaf Newton step, so that Newton runs in the overflow
                     // pass, with all lanes of a warp at once, instead of one lane stalling the rest.
                     budget_steps = env_int("MANDELBROT_CUDA_BUDGET", 1 << 14),
                     timing = env_int("MANDELBROT_CUDA_TIMING", 0),
                     park = env_int("MANDELBROT_CUDA_PARK", 8),  // Room to park n / park orbits
                     chunk = env_int("MANDELBROT_CUDA_CHUNK", 16);  // Items per claim
    const int64_t cap = std::min<int64_t>(n, std::max<int64_t>(1 << 16, n / park));
    Mem<O> parked(cap, true), next(cap, true);
    Mem<int64_t> items(cap, true), next_items(cap, true);
    Mem<unsigned long long> counters(13, true);
    counters.zero();
    cudaEvent_t e0, e1, e2;
    cuda_check(cudaEventCreate(&e0)); cuda_check(cudaEventCreate(&e1)); cuda_check(cudaEventCreate(&e2));
    cuda_check(cudaEventRecord(e0, stream()));
    // Register budget: 65536 / (256 · min_blocks) per thread
    const int min_blocks = min_blocks_env ? min_blocks_env : task.min_blocks;
    const int grid = blocks_per_sm * num_sms(), threads = grid * block;
    const int32_t budget = int32_t(std::max<int64_t>(1, budget_steps / task.burst));
    if (timing)
      engine_detail::launch_orbit_kernel<Task, true>(min_blocks, grid, task, n, stride, budget, parked.p, items.p,
                                                     cap, counters.p, chunk);
    else
      engine_detail::launch_orbit_kernel<Task, false>(min_blocks, grid, task, n, stride, budget, parked.p, items.p,
                                                      cap, counters.p, chunk);
    cuda_check(cudaGetLastError());
    cuda_check(cudaEventRecord(e1, stream()));
    // Rounds: settle pending orbits together, then resume the rest until they are done or pending again
    int64_t count = std::min<int64_t>(cap, int64_t(counters.get(2))), rounds = 0;
    stats.overflow = count;
    while (count) {
      rounds++;
      const int g = int(std::min<int64_t>(grid, (count + block - 1) / block));
      engine_detail::settle_kernel<Task><<<g, block, 0, stream()>>>(task, parked.p, items.p, count, counters.p);
      cuda_check(cudaMemsetAsync(counters.p + 8, 0, sizeof(unsigned long long), stream()));
      cuda_check(cudaMemsetAsync(counters.p + 12, 0, sizeof(unsigned long long), stream()));
      engine_detail::resume_kernel<Task><<<g, block, 0, stream()>>>(task, parked.p, items.p, count, next.p,
                                                                   next_items.p, cap, counters.p);
      cuda_check(cudaGetLastError());
      count = std::min<int64_t>(cap, int64_t(counters.get(12)));
      std::swap(parked.p, next.p);
      std::swap(items.p, next_items.p);
    }
    cuda_check(cudaEventRecord(e2, stream()));
    unsigned long long h[13];
    counters.to_host(h, 13);
    stats.iters = int64_t(h[1]);
    if (timing) {
      float main_ms, over_ms;
      cuda_check(cudaEventElapsedTime(&main_ms, e0, e1));
      cuda_check(cudaEventElapsedTime(&over_ms, e1, e2));
      print("    cuda run: %d items, %d threads, main %.1f ms, %d parked, %d rounds %.1f ms, %.3g it/s; "
            "iterations per thread mean %.3g, max %.3g", n, threads, main_ms, stats.overflow, rounds, over_ms,
            double(h[1]) / ((main_ms + over_ms) * 1e-3), double(h[1]) / threads, double(h[3]));
      print("      main pass: %.3g it/s, %.1f%% of thread cycles in run, SIMT efficiency at run %.1f%%, "
            "lane steps / warp steps %.1f%%; rounds %.3g it/s", double(h[1] - h[9]) / (main_ms * 1e-3),
            100.0 * double(h[4]) / double(h[5]), 100.0 * double(h[6]) / double(h[7]),
            100.0 * double(h[10]) / double(h[11]), over_ms > 0 ? double(h[9]) / (over_ms * 1e-3) : 0.0);
    }
    cuda_check(cudaEventDestroy(e0)); cuda_check(cudaEventDestroy(e1)); cuda_check(cudaEventDestroy(e2));
#else
    die("run_orbits: built without CUDA");
#endif
  }
  stats.secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return stats;
}

// Call f(i) for i < n in parallel; f must be __host__ __device__ and trivially copyable
template<class F> void for_each(const int64_t n, const F& f, const bool cuda) {
  if (!n) return;
  if (!cuda) {
    const int threads = cpu_threads();
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; t++)
      pool.emplace_back([&, t]() {
        const int64_t lo = n * t / threads, hi = n * (t + 1) / threads;
        for (int64_t i = lo; i < hi; i++) f(i);
      });
    for (auto& t : pool) t.join();
  } else {
#ifdef __CUDACC__
    engine_detail::for_each_kernel<F><<<8 * num_sms(), 256, 0, stream()>>>(n, f);
    cuda_check(cudaGetLastError());
#else
    die("for_each: built without CUDA");
#endif
  }
}

}  // namespace mandelbrot
