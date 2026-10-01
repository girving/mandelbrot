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
#include <typeinfo>
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
#ifdef __CUDACC__
// Atomic max, on the device.  (CUDA's atomics take unsigned long long, which is not uint64_t on Linux, so the
// casts live only in these helpers.)
__device__ static inline void atomic_max(uint64_t* p, const uint64_t v) {
  atomicMax(reinterpret_cast<unsigned long long*>(p), static_cast<unsigned long long>(v));
}
#endif

// A buffer of trivially copyable T on the host or the device
template<class T> struct Mem : public Noncopyable {
  T* p = nullptr;
  int64_t n = 0;
  bool cuda = false;

  Mem(const int64_t n, const bool cuda) : n(n), cuda(cuda) {
    static_assert(std::is_trivially_copyable_v<T>);
    if (!n) return;
    if (cuda) {
      IF_CUDA(p = static_cast<T*>(cuda_malloc(n * sizeof(T))));
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
    if (cuda) {
      IF_CUDA(cuda_check(cudaMemcpyAsync(dst, p, count * sizeof(T), cudaMemcpyDeviceToHost, stream()));
              cuda_check(cudaStreamSynchronize(stream())));
    } else {
      memcpy(dst, p, count * sizeof(T));
    }
  }
  void from_host(const T* src, const int64_t count) {
    slow_assert(count <= n);
    if (!count) return;
    if (cuda) {
      IF_CUDA(cuda_check(cudaMemcpyAsync(p, src, count * sizeof(T), cudaMemcpyHostToDevice, stream()));
              cuda_check(cudaStreamSynchronize(stream())));
    } else {
      memcpy(p, src, count * sizeof(T));
    }
  }
  T get(const int64_t i) const {
    T x;
    slow_assert(i < n);
    if (cuda) {
      IF_CUDA(cuda_check(cudaMemcpyAsync(&x, p + i, sizeof(T), cudaMemcpyDeviceToHost, stream()));
              cuda_check(cudaStreamSynchronize(stream())));
    } else {
      x = p[i];
    }
    return x;
  }
};

// Statistics of one run_orbits call
struct RunStats {
  int64_t iters = 0;     // Total iterations
  int64_t overflow = 0;  // Orbits deferred to the second pass (GPU only)
  int64_t cpu_tail = 0;  // Parked orbits finished on CPU threads (GPU only)
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
// 9: iterations in rounds, 10-11 with timing: lane steps and warp steps, 12: parked in this round, 13-16 with
// timing: resume_kernel's cycles inside run, total cycles, lane steps and warp steps].

// Per-warp queue of ready orbits in shared memory, restocked 32 at a time (one per lane), so that lanes going idle
// in different bursts copy a ready orbit instead of each preparing one (starting a sample: a few hundred
// instructions of placement, hashing, and the cardioid test) with the rest of the warp waiting.  (Not for
// resume_kernel: restocking 32 at a time hoards parked orbits in busy warps at the end of each round.)  Item
// indices fit in 32 bits (n < 2^31), which saves registers.
template<class State> struct WarpQueue {
  State* slots;            // This warp's 32 slots
  int32_t* items;
  unsigned ready = 0;      // Slots holding ready orbits (warp-uniform)
  bool exhausted = false;  // Nothing left to restock from (warp-uniform)

  // Give idle (done) lanes ready orbits while there are any.  restock(slot, item, filled) fills this lane's slot,
  // setting filled if it holds a ready orbit, and returns whether the source is now exhausted (warp-uniform).
  template<class Restock> __device__ void refill(bool& done, State& o, int32_t& i, int32_t& bursts,
                                                 const Restock& restock) {
    const int lane = int(threadIdx.x & 31);
    for (unsigned need = __ballot_sync(0xffffffff, done); need && (ready || !exhausted);) {
      if (!ready) {
        bool filled = false;
        exhausted = restock(slots[lane], items[lane], filled);
        ready = __ballot_sync(0xffffffff, filled);
        __syncwarp();
        continue;
      }
      // The r-th idle lane takes the r-th ready slot
      const int r = __popc(need & ((1u << lane) - 1)), avail = __popc(ready);
      if (((need >> lane) & 1) && r < avail) {
        const int slot = int(__fns(ready, 0, r + 1));
        o = slots[slot];
        i = items[slot];
        done = false;
        bursts = 0;
      }
      const int used = min(__popc(need), avail);
      ready = used < avail ? ready & ~((1u << __fns(ready, 0, used + 1)) - 1) : 0;
      need = __ballot_sync(0xffffffff, done);
      __syncwarp();
    }
  }
};

// A task's compact finish record (Task::Record, with record and finish_record), or char if it has none
template<class Task> struct RecordOf { typedef char type; };
template<class Task> requires requires { typename Task::Record; } struct RecordOf<Task> {
  typedef typename Task::Record type;
};

// Persistent threads.  Orbits that are pending or have run `budget` bursts are parked (if room) and finished in
// rounds of settle_kernel and resume_kernel.
template<class Task, bool timing, int min_blocks> __global__ void __launch_bounds__(256, min_blocks)
orbit_kernel(const Task task, const int64_t n, const int64_t stride, const int32_t budget,
             typename Task::State* overflow, int64_t* overflow_items, const int64_t overflow_cap,
             uint64_t* counters) {
  // Lanes stay in the loop until their whole warp is out of work, and reconverge before each burst, so that
  // lanes refilling at different times do not split the warp into groups that each step half empty.  Refills
  // come from a WarpQueue, restocked by claiming 32 items and starting them in place (so that no second state
  // occupies registers); items decided at once (the cardioid) finish there.
  //
  // Tasks with a Record buffer their finishes the same way: done lanes store compact records, and the warp
  // finishes 32 at once (when the buffer would overflow, and at the end), instead of finishing a few lanes per
  // burst with the rest idle.  Only if the buffer fits in static shared memory alongside the queue.
  typedef typename Task::State State;
  typedef typename RecordOf<Task>::type Record;
  constexpr bool buffered = requires { typename Task::Record; } &&
                            256 * (sizeof(State) + sizeof(Record) + 8) <= 48 * 1024;
  __shared__ alignas(16) unsigned char queue_bytes[256 * sizeof(State)];
  __shared__ int32_t queue_items[256];
  __shared__ alignas(16) unsigned char record_bytes[buffered ? 256 * sizeof(Record) : 16];
  __shared__ int32_t record_items[buffered ? 256 : 1];
  const int lane = int(threadIdx.x & 31);
  WarpQueue<State> queue{reinterpret_cast<State*>(queue_bytes) + (threadIdx.x & ~31u), queue_items + (threadIdx.x & ~31u)};
  [[maybe_unused]] Record* const records = reinterpret_cast<Record*>(record_bytes) + (buffered ? threadIdx.x & ~31u : 0);
  [[maybe_unused]] int32_t* const ritems = record_items + (buffered ? threadIdx.x & ~31u : 0);
  [[maybe_unused]] int records_n = 0;  // Buffered records (warp-uniform)
  const auto flush = [&]() {
    if constexpr (buffered) {
      if (lane < records_n) task.finish_record(records[lane], ritems[lane]);
      records_n = 0;
      __syncwarp();
    }
  };
  State o;
  int32_t i = -1, bursts = 0;  // i = current item or -1
  uint64_t it = 0, run_cycles = 0, active = 0, slots = 0, lane_steps = 0, warp_steps = 0;
  const long long t0 = timing ? clock64() : 0;
  bool done = true;
  const auto restock = [&](State& q, int32_t& qi, bool& filled) {
    uint64_t j0 = 0;
    if (!lane) j0 = atomic_fetch_add(counters, uint64_t(32));
    j0 = __shfl_sync(0xffffffff, j0, 0);
    const int64_t j = int64_t(j0) + lane;
    if (j < n) {
      const int32_t k = int32_t(scramble(j, stride, n));
      if (task.start(q, k)) {  // Decided at once
        task.finish(q, k);
        it += task.iters(q);
      } else {
        qi = k;
        filled = true;
      }
    }
    return int64_t(j0) + 32 >= n;
  };
  for (;;) {
    if constexpr (buffered) {
      const bool has = done && i >= 0;
      const unsigned m = __ballot_sync(0xffffffff, has);
      if (m) {
        if (records_n + __popc(m) > 32) flush();
        if (has) {
          const int r = records_n + __popc(m & ((1u << lane) - 1));
          records[r] = task.record(o);
          ritems[r] = i;
          it += task.iters(o);
          i = -1;
        }
        records_n += __popc(m);
        __syncwarp();
      }
    } else if (done && i >= 0) {
      task.finish(o, i);
      it += task.iters(o);
      i = -1;
    }
    queue.refill(done, o, i, bursts, restock);
    if (__all_sync(0xffffffff, done)) break;  // Out of work
    if (!done) {
      if constexpr (timing) {
        const unsigned mask = __activemask();
        const bool leader = (threadIdx.x & 31) == __ffs(mask) - 1;
        if (leader) { active += __popc(mask); slots += 32; }
        const int64_t p0 = task.progress(o);
        const long long r0 = clock64();
        done = task.run(o);
        if constexpr (requires { task.immediate(o); }) if (done && task.immediate(o)) done = task.settle(o);
        run_cycles += clock64() - r0;
        // Steps this lane took, against the warp's longest: idle lanes within bursts
        const unsigned steps = unsigned(task.progress(o) - p0), longest = __reduce_max_sync(mask, steps);
        lane_steps += steps;
        if (leader) warp_steps += uint64_t(32) * longest;
      } else {
        done = task.run(o);
        // Settles due at once (an overflowed block), for all such lanes of the warp together
        if constexpr (requires { task.immediate(o); }) if (done && task.immediate(o)) done = task.settle(o);
      }
      // Park orbits that are pending (so that lanes settle together later) or have used their step budget.
      // Parked orbits' iterations are counted when they finish.
      const bool pending = done && task.pending(o);
      if (pending || (!done && ++bursts == budget)) {
        const int64_t k = int64_t(atomic_fetch_add(counters + 2, 1));
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
    __syncwarp();
  }
  flush();
  atomic_add(counters + 1, it);
  atomic_max(counters + 3, it);
  if constexpr (timing) {
    atomic_add(counters + 4, run_cycles);
    atomic_add(counters + 5, uint64_t(clock64() - t0));
    atomic_add(counters + 6, active);
    atomic_add(counters + 7, slots);
    atomic_add(counters + 10, lane_steps);
    atomic_add(counters + 11, warp_steps);
  }
}

// Launch orbit_kernel with register pressure chosen at run time: min_blocks resident 256-thread blocks per SM
template<class Task, bool timing> void launch_orbit_kernel(const int min_blocks, const int blocks, const Task& task,
                                                           const int64_t n, const int64_t stride, const int32_t budget,
                                                           typename Task::State* overflow, int64_t* items,
                                                           const int64_t cap, uint64_t* counters) {
#define LAUNCH(b) orbit_kernel<Task, timing, b><<<blocks, 256, 0, stream()>>>(task, n, stride, budget, overflow, \
                                                                             items, cap, counters)
  switch (min_blocks) {
    case 1: LAUNCH(1); break;
    case 2: LAUNCH(2); break;
    case 3: LAUNCH(3); break;
    case 4: LAUNCH(4); break;
    default: die("MANDELBROT_CUDA_MIN_BLOCKS must be 1 to 4, got %d", min_blocks);
  }
#undef LAUNCH
}

// Counting sort of parked orbits by task.settle_key (the work settle will do, in [0, kSettleKeys)) into a
// permutation, so that each warp of settle_kernel settles similar work: settles are long serial loops (Newton over
// the candidate period) that would otherwise run at the pace of each warp's slowest lane.  Each block ranks a tile
// of keys in shared memory and reserves each bucket's range with one global atomic per (tile, bucket).
constexpr int kSettleKeys = 512, kSortTile = 256 * 8;
template<class Task> __global__ void settle_keys_kernel(const Task task, const typename Task::State* parked,
                                                        const int64_t* items, const int64_t count, uint16_t* keys,
                                                        uint64_t* sizes) {
  __shared__ uint32_t hist[kSettleKeys];
  for (int b = threadIdx.x; b < kSettleKeys; b += blockDim.x) hist[b] = 0;
  __syncthreads();
  for (int64_t k = int64_t(blockIdx.x) * blockDim.x + threadIdx.x; k < count; k += int64_t(blockDim.x) * gridDim.x) {
    const int key = items[k] < 0 ? 0 : min(max(task.settle_key(parked[k]), 0), kSettleKeys - 1);
    keys[k] = uint16_t(key);
    atomicAdd(hist + key, 1u);
  }
  __syncthreads();
  for (int b = threadIdx.x; b < kSettleKeys; b += blockDim.x)
    if (hist[b]) atomic_add(sizes + b, uint64_t(hist[b]));
}
template<class Task> __global__ void settle_sort_kernel(const uint16_t* keys, const int64_t count, uint64_t* offsets,
                                                        int32_t* perm) {
  __shared__ uint32_t hist[kSettleKeys];
  __shared__ uint64_t base[kSettleKeys];
  constexpr int per = kSortTile / 256;
  for (int64_t t0 = int64_t(blockIdx.x) * kSortTile; t0 < count; t0 += int64_t(gridDim.x) * kSortTile) {
    for (int b = threadIdx.x; b < kSettleKeys; b += blockDim.x) hist[b] = 0;
    __syncthreads();
    uint32_t rank[per];
    int key[per];
    for (int j = 0; j < per; j++) {
      const int64_t k = t0 + j * 256 + threadIdx.x;
      key[j] = k < count ? keys[k] : -1;
      if (key[j] >= 0) rank[j] = atomicAdd(hist + key[j], 1u);
    }
    __syncthreads();
    for (int b = threadIdx.x; b < kSettleKeys; b += blockDim.x)
      if (hist[b]) base[b] = atomic_fetch_add(offsets + b, uint64_t(hist[b]));
    __syncthreads();
    for (int j = 0; j < per; j++)
      if (key[j] >= 0) perm[base[key[j]] + rank[j]] = int32_t(t0 + j * 256 + threadIdx.x);
    __syncthreads();
  }
}

// Settle all pending parked orbits at once (in the order perm, if given), finishing those that are done (marked by
// item -1)
template<class Task> __global__ void settle_kernel(const Task task, typename Task::State* parked, int64_t* items,
                                                   const int64_t count, uint64_t* counters,
                                                   const int32_t* perm = nullptr) {
  uint64_t it = 0;
  for (int64_t t = int64_t(blockIdx.x) * blockDim.x + threadIdx.x; t < count; t += int64_t(blockDim.x) * gridDim.x) {
    const int64_t k = perm ? perm[t] : t;
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
  atomic_add(counters + 1, it);
  atomic_add(counters + 9, it);
}

// Resume parked orbits (skipping finished ones) with persistent warp-synchronous threads like orbit_kernel's,
// until done or pending again, when they park into next
template<class Task, bool timing> __global__ void resume_kernel(const Task task, const typename Task::State* parked,
                                                   const int64_t* items, const int64_t count,
                                                   typename Task::State* next, int64_t* next_items, const int64_t cap,
                                                   uint64_t* counters) {
  typename Task::State o;
  int64_t i = -1;
  uint64_t it = 0, run_cycles = 0, lane_steps = 0, warp_steps = 0;
  const long long t0 = timing ? clock64() : 0;
  bool done = true, out = false;
  for (;;) {
    while (done && !out) {
      if (i >= 0) { task.finish(o, i); it += task.iters(o); }
      i = -1;
      const int64_t k = int64_t(atomic_fetch_add(counters + 8, 1));
      if (k >= count) { out = true; break; }
      if (items[k] < 0) continue;
      o = parked[k];
      i = items[k];
      done = false;
    }
    if (__all_sync(0xffffffff, out)) break;
    __syncwarp();
    if (!out) {
      const unsigned mask = timing ? __activemask() : 0;
      const int64_t p0 = timing ? task.progress(o) : 0;
      const long long r0 = timing ? clock64() : 0;
      done = task.run(o);
      if constexpr (timing) {
        run_cycles += clock64() - r0;
        const unsigned steps = unsigned(task.progress(o) - p0), longest = __reduce_max_sync(mask, steps);
        lane_steps += steps;
        if ((threadIdx.x & 31) == __ffs(mask) - 1) warp_steps += uint64_t(32) * longest;
      }
      if constexpr (requires { task.immediate(o); }) if (done && task.immediate(o)) done = task.settle(o);
      if (done && task.pending(o)) {
        const int64_t k = int64_t(atomic_fetch_add(counters + 12, 1));
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
  atomic_add(counters + 1, it);
  atomic_add(counters + 9, it);
  if constexpr (timing) {
    atomic_add(counters + 13, run_cycles);
    atomic_add(counters + 14, uint64_t(clock64() - t0));
    atomic_add(counters + 15, lane_steps);
    atomic_add(counters + 16, warp_steps);
  }
}

// Finish orbits whose states were completed on the host (CPU tail), skipping finished ones (item -1)
template<class Task> __global__ void finish_kernel(const Task task, const typename Task::State* states,
                                                   const int64_t* items, const int64_t count) {
  for (int64_t k = int64_t(blockIdx.x) * blockDim.x + threadIdx.x; k < count; k += int64_t(blockDim.x) * gridDim.x)
    if (items[k] >= 0) task.finish(states[k], items[k]);
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
                     budget_env = env_int("MANDELBROT_CUDA_BUDGET", -1),
                     timing = env_int("MANDELBROT_CUDA_TIMING", 0),
                     park = env_int("MANDELBROT_CUDA_PARK", 8),  // Room to park n / park orbits
                     // Finish on CPU threads once this few orbits remain: a lone GPU lane steps a sequential
                     // orbit ~40× slower than a CPU core, so late rounds of a few long orbits idle the GPU.
                     cpu_tail = env_int("MANDELBROT_CUDA_CPU_TAIL", 1024);
    const int64_t cap = std::min<int64_t>(n, std::max<int64_t>(1 << 16, n / park));
    Mem<O> parked(cap, true), next(cap, true);
    Mem<int64_t> items(cap, true), next_items(cap, true);
    Mem<uint64_t> counters(17, true);
    counters.zero();
    cudaEvent_t e0, e1, e2;
    cuda_check(cudaEventCreate(&e0)); cuda_check(cudaEventCreate(&e1)); cuda_check(cudaEventCreate(&e2));
    cuda_check(cudaEventRecord(e0, stream()));
    // Register budget: 65536 / (256 · min_blocks) per thread
    const int min_blocks = min_blocks_env ? min_blocks_env : task.min_blocks;
    const int grid = blocks_per_sm * num_sms(), threads = grid * block;
    // Steps before parking: the task's park_steps (its first Newton step, so that Newton runs in the rounds with
    // all lanes of a warp at once), overridden by MANDELBROT_CUDA_BUDGET, else 2^14
    int64_t budget_steps = int64_t(1) << 14;
    if constexpr (requires { task.park_steps; }) if (task.park_steps > 0) budget_steps = task.park_steps;
    if (budget_env > 0) budget_steps = budget_env;
    const int32_t budget = int32_t(std::max<int64_t>(1, budget_steps / task.burst));
    if (timing)
      engine_detail::launch_orbit_kernel<Task, true>(min_blocks, grid, task, n, stride, budget, parked.p, items.p,
                                                     cap, counters.p);
    else
      engine_detail::launch_orbit_kernel<Task, false>(min_blocks, grid, task, n, stride, budget, parked.p, items.p,
                                                      cap, counters.p);
    cuda_check(cudaGetLastError());
    cuda_check(cudaEventRecord(e1, stream()));
    // Rounds: settle pending orbits together, then resume the rest until they are done or pending again
    int64_t count = std::min<int64_t>(cap, int64_t(counters.get(2))), rounds = 0;
    stats.overflow = count;
    int64_t cpu_iters = 0;
    while (count) {
      if (count <= cpu_tail) {
        // Few orbits left: copy them to the host, finish them on CPU threads, and write results on the device
        std::vector<O> h(count);
        std::vector<int64_t> hi(count);
        parked.to_host(h.data(), count);
        items.to_host(hi.data(), count);
        std::atomic<int64_t> next_k(0), iters(0);
        std::vector<std::thread> pool;
        for (int t = 0; t < cpu_threads(); t++)
          pool.emplace_back([&]() {
            int64_t it = 0;
            for (int64_t k; (k = next_k.fetch_add(1)) < count;) {
              if (hi[k] < 0) continue;
              O o = h[k];
              for (bool done = false; !done;) {
                if (task.pending(o)) done = task.settle(o);
                else done = task.run(o) && !task.pending(o);
              }
              h[k] = o;
              it += task.iters(o);
            }
            iters += it;
          });
        for (auto& t : pool) t.join();
        parked.from_host(h.data(), count);
        const int g = int(std::min<int64_t>(grid, (count + block - 1) / block));
        engine_detail::finish_kernel<Task><<<g, block, 0, stream()>>>(task, parked.p, items.p, count);
        cuda_check(cudaGetLastError());
        cpu_iters = iters;
        stats.cpu_tail = count;
        break;
      }
      rounds++;
      const int g = int(std::min<int64_t>(grid, (count + block - 1) / block));
      const auto r0 = std::chrono::steady_clock::now();
      if constexpr (requires { task.settle_key(*parked.p); }) {
        // Sort by settle work (see settle_keys_kernel), then settle in that order
        Mem<uint16_t> keys(count, true);
        Mem<int32_t> perm(count, true);
        Mem<uint64_t> sizes(engine_detail::kSettleKeys, true);
        sizes.zero();
        engine_detail::settle_keys_kernel<Task><<<g, block, 0, stream()>>>(task, parked.p, items.p, count, keys.p,
                                                                           sizes.p);
        uint64_t h[engine_detail::kSettleKeys];
        sizes.to_host(h, engine_detail::kSettleKeys);
        for (uint64_t b = 0, offset = 0; b < uint64_t(engine_detail::kSettleKeys); b++) {
          const uint64_t size = h[b];
          h[b] = offset;
          offset += size;
        }
        sizes.from_host(h, engine_detail::kSettleKeys);
        const int sg = int(std::min<int64_t>(grid, (count + engine_detail::kSortTile - 1) / engine_detail::kSortTile));
        engine_detail::settle_sort_kernel<Task><<<sg, block, 0, stream()>>>(keys.p, count, sizes.p, perm.p);
        engine_detail::settle_kernel<Task><<<g, block, 0, stream()>>>(task, parked.p, items.p, count, counters.p,
                                                                      perm.p);
      } else {
        engine_detail::settle_kernel<Task><<<g, block, 0, stream()>>>(task, parked.p, items.p, count, counters.p);
      }
      if (timing) cuda_check(cudaStreamSynchronize(stream()));
      const auto r1 = std::chrono::steady_clock::now();
      cuda_check(cudaMemsetAsync(counters.p + 8, 0, sizeof(uint64_t), stream()));
      cuda_check(cudaMemsetAsync(counters.p + 12, 0, sizeof(uint64_t), stream()));
      if (timing)
        engine_detail::resume_kernel<Task, true><<<g, block, 0, stream()>>>(task, parked.p, items.p, count, next.p,
                                                                            next_items.p, cap, counters.p);
      else
        engine_detail::resume_kernel<Task, false><<<g, block, 0, stream()>>>(task, parked.p, items.p, count, next.p,
                                                                             next_items.p, cap, counters.p);
      cuda_check(cudaGetLastError());
      const int64_t was = count;
      count = std::min<int64_t>(cap, int64_t(counters.get(12)));
      if (timing) {
        const auto r2 = std::chrono::steady_clock::now();
        print("      round %d: %d orbits, settle %.1f ms, resume %.1f ms, %d still pending", rounds, was,
              1e3 * std::chrono::duration<double>(r1 - r0).count(), 1e3 * std::chrono::duration<double>(r2 - r1).count(),
              count);
      }
      std::swap(parked.p, next.p);
      std::swap(items.p, next_items.p);
    }
    cuda_check(cudaEventRecord(e2, stream()));
    uint64_t h[17];
    counters.to_host(h, 17);
    stats.iters = int64_t(h[1]) + cpu_iters;
    if (timing) {
      float main_ms, over_ms;
      cuda_check(cudaEventElapsedTime(&main_ms, e0, e1));
      cuda_check(cudaEventElapsedTime(&over_ms, e1, e2));
      print("    cuda run (%s): %d items, %d threads, main %.1f ms, %d parked, %d rounds %.1f ms, %.3g it/s; "
            "iterations per thread mean %.3g, max %.3g", typeid(Task).name(), n, threads, main_ms, stats.overflow, rounds, over_ms,
            double(h[1]) / ((main_ms + over_ms) * 1e-3), double(h[1]) / threads, double(h[3]));
      print("      main pass: %.3g it/s, %.1f%% of thread cycles in run, SIMT efficiency at run %.1f%%, "
            "lane steps / warp steps %.1f%%; rounds %.3g it/s", double(h[1] - h[9]) / (main_ms * 1e-3),
            100.0 * double(h[4]) / double(h[5]), 100.0 * double(h[6]) / double(h[7]),
            100.0 * double(h[10]) / double(h[11]), over_ms > 0 ? double(h[9]) / (over_ms * 1e-3) : 0.0);
      if (h[14])
        print("      resume: %.1f%% of thread cycles in run, lane steps / warp steps %.1f%%, %.3g thread cycles",
              100.0 * double(h[13]) / double(h[14]), 100.0 * double(h[15]) / double(h[16]), double(h[14]));
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
