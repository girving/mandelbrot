// Parallel execution of resumable orbits and simple loops, on CPU threads or the GPU
//
// The same task code runs in both places, so CPU runs are exact references for GPU runs.  A task describes
// n independent items, each an orbit (Orbit, OrbitDE, ...) that is started, advanced in short bursts, and
// finished.  Workers claim items in scrambled order (slow orbits cluster spatially, so consecutive claims
// should land far apart) and refill each lane or GPU thread as soon as its orbit finishes.
//
// On the GPU, each lane's orbit is a small state machine inside one persistent kernel (pool_kernel): started
// from a per-warp queue of ready orbits, run in bursts, then finished (through a per-warp buffer of compact
// records if the task has them), or, when it stops for work that should run together with other lanes (a
// Newton step), pushed onto a device ring keyed by that work (settle_key: Newton's period).  Warps pop 32
// orbits with the same key, settle them together, and run on the ones that continue.  No host rounds or
// barriers: the kernel ends when no items, rings, or running orbits remain, or, when few remain, hands them
// to CPU threads, which step a sequential orbit far faster than a lone GPU lane.
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
//   int min_blocks;                          // GPU: resident 256-thread blocks per SM (1 to 4), budgeting registers
// Optional: int settle_key(State& o) (in [0, kSettleKeys), grouping similar settles; it may cache work in o),
// bool immediate(const State& o) (a pending orbit the lane should settle at once), and a compact finish record
// (typedef Record, Record record(const State&), void finish_record(const Record&, int64_t i)).
#pragma once

#include "cutil.h"
#include "debug.h"
#include "noncopyable.h"
#include "print.h"
#include "rings.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <memory>
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
  void to_host(T* dst, const int64_t count, const int64_t offset = 0) const {
    slow_assert(offset + count <= n);
    if (!count) return;
    if (cuda) {
      IF_CUDA(cuda_check(cudaMemcpyAsync(dst, p + offset, count * sizeof(T), cudaMemcpyDeviceToHost, stream()));
              cuda_check(cudaStreamSynchronize(stream())));
    } else {
      memcpy(dst, p + offset, count * sizeof(T));
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
  int64_t overflow = 0;  // Orbits pushed onto settle rings (GPU only)
  int64_t cpu_tail = 0;  // Orbits finished on CPU threads (GPU only)
  double secs = 0;
};

// Settle keys (Task::settle_key), and the bound on orbits waiting in settle queues at once
constexpr int kSettleKeys = 512;
constexpr int64_t kSettleBound = 1 << 15;

static inline uint64_t next_pow2(const uint64_t x) {
  uint64_t p = 1;
  while (p < x) p *= 2;
  return p;
}

#ifdef __CUDACC__
namespace engine_detail {
// An IntRing holding 0, ..., count - 1 (kernels in this header are templates, so that each program has one copy)
template<class Ring> __global__ void int_ring_reset(const Ring r, const int64_t count) {
  for (int64_t i = int64_t(blockIdx.x) * blockDim.x + threadIdx.x; i < int64_t(r.cap);
       i += int64_t(blockDim.x) * gridDim.x) {
    r.slots[i] = int32_t(i);
    r.seq[i] = uint64_t(i) + (i < count);
    if (!i) { r.ends[0] = 0; r.ends[1] = uint64_t(count); *r.avail = count; }
  }
}
template<class I> __global__ void set_int(I* p, const I v) { *p = v; }
// out[i] = the value at position pos[i] of queue key[i]
template<class T> __global__ void queues_gather(const Queues<T> q, const int32_t* key, const uint64_t* pos,
                                                const int64_t n, T* out) {
  for (int64_t i = int64_t(blockIdx.x) * blockDim.x + threadIdx.x; i < n; i += int64_t(blockDim.x) * gridDim.x) {
    const uint64_t p = pos[i], x = q.table[uint64_t(key[i]) * q.entries + ((p / Queues<T>::S) & (q.entries - 1))];
    out[i] = q.values[(x & 0xffffffff) * Queues<T>::S + p % Queues<T>::S];
  }
}
}  // namespace engine_detail

// Device storage for an IntRing of capacity cap; reset(count) makes it hold 0, ..., count - 1
struct IntRingMem : public Noncopyable {
  const uint64_t cap;
  Mem<int32_t> slots;
  Mem<uint64_t> seq, ends;
  Mem<int64_t> avail;
  explicit IntRingMem(const uint64_t cap)
    : cap(cap), slots(int64_t(cap), true), seq(int64_t(cap), true), ends(2, true), avail(1, true) {
    slow_assert(cap && !(cap & (cap - 1)));
  }
  void reset(const int64_t count = 0) {
    slow_assert(uint64_t(count) <= cap);
    engine_detail::int_ring_reset<<<int(std::min<uint64_t>(1024, (cap + 255) / 256)), 256, 0, stream()>>>(ring(),
                                                                                                         count);
    cuda_check(cudaGetLastError());
  }
  IntRing ring() const { return {slots.p, seq.p, ends.p, avail.p, cap}; }
};

// Device storage for Queues: keys queues with at most bound values outstanding, with segments and table entries
// to spare.  Allocated once and reset (on the device) for each use.
template<class T> struct QueuesMem : public Noncopyable {
  static constexpr int S = Queues<T>::S;
  const int keys;
  const int64_t bound, segments;
  const uint64_t entries;
  Mem<T> values;
  Mem<uint64_t> vseq, table, ends;  // ends: [head[keys], tail[keys]]
  Mem<uint32_t> consumed, nonempty;
  Mem<int64_t> avail, space;
  IntRingMem free;

  QueuesMem(const int keys, const int64_t bound)
    : keys(keys), bound(bound), segments(bound + 256), entries(next_pow2(4 * bound / S)),
      values(segments * S, true), vseq(segments * S, true), table(keys * int64_t(entries), true),
      ends(2 * keys, true), consumed(segments, true), nonempty((keys + 31) / 32, true), avail(keys, true),
      space(1, true), free(next_pow2(uint64_t(segments))) {
    reset();
  }

  // Empty every queue, with every segment free
  void reset() {
    vseq.zero();
    consumed.zero();
    nonempty.zero();
    ends.zero();
    avail.zero();
    cuda_check(cudaMemsetAsync(table.p, 0xff, keys * entries * sizeof(uint64_t), stream()));
    engine_detail::set_int<int64_t><<<1, 1, 0, stream()>>>(space.p, bound);
    free.reset(segments);
  }

  Queues<T> queues() const {
    return {values.p, vseq.p, table.p, consumed.p, ends.p, ends.p + keys, avail.p, space.p, nonempty.p,
            free.ring(), keys, entries};
  }

  // Every value still queued, once the kernels using the queues have finished (gathered on the device)
  std::vector<T> contents() const {
    std::vector<uint64_t> e(2 * keys);
    ends.to_host(e.data(), 2 * keys);
    std::vector<int32_t> key;
    std::vector<uint64_t> pos;
    for (int k = 0; k < keys; k++)
      for (uint64_t p = e[k]; p < e[keys + k]; p++) { key.push_back(k); pos.push_back(p); }
    const int64_t n = int64_t(key.size());
    std::vector<T> all(n);
    if (!n) return all;
    Mem<int32_t> dk(n, true);
    Mem<uint64_t> dp(n, true);
    Mem<T> out(n, true);
    dk.from_host(key.data(), n);
    dp.from_host(pos.data(), n);
    engine_detail::queues_gather<<<int(std::min<int64_t>(1024, (n + 255) / 256)), 256, 0, stream()>>>(
        queues(), dk.p, dp.p, n, out.p);
    cuda_check(cudaGetLastError());
    out.to_host(all.data(), n);
    return all;
  }
};
#endif  // __CUDACC__

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

// Counters: [0: next item claim, 1: iterations, 2: active warps, 3: drain flag, 4: ring pushes, 5: settles alone
// (a full ring), 6: orbits dumped to the CPU tail, 7-12 with timing: cycles inside run, total cycles, lane steps,
// warp steps, cycles settling popped orbits, cycles idle, 13-15 with timing: cycles refilling, pushing, finishing].
constexpr int kCounters = 18;

// An orbit and its item, as the settle rings hold them
template<class State> struct Parked {
  State o;
  int32_t item;
};

// Per-warp queue of ready orbits in shared memory, restocked 32 at a time (one per lane), so that lanes going idle
// in different bursts copy a ready orbit instead of each preparing one (starting a sample: a few hundred
// instructions of placement, hashing, and the cardioid test; or settling a popped orbit) with the rest of the warp
// waiting.  Item indices fit in 32 bits (n < 2^31), which saves registers.
template<class State> struct WarpQueue {
  State* slots;        // This warp's 32 slots
  int32_t* items;
  unsigned ready = 0;  // Slots holding ready orbits (warp-uniform)

  // Give idle (done) lanes ready orbits while there are any.  restock(need), warp-collective, fills slots (each
  // lane its own) and ready, and returns false if it found no work.
  template<class Restock> __device__ void refill(bool& done, State& o, int32_t& i, const Restock& restock) {
    const int lane = int(threadIdx.x & 31);
    for (unsigned need = __ballot_sync(0xffffffff, done); need;) {
      if (!ready) {
        if (!restock(need)) return;
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

// The engine: persistent threads, each lane's orbit a state machine (see the top of this file).  Lanes stay in
// the loop until their warp is out of work, and reconverge before each burst, so that lanes refilling at
// different times do not split the warp into groups that each step half empty.
//
// Refills: a WarpQueue restocked, in order of preference, by popping a full batch (32 orbits with the same
// settle key, announced on the ripe queue by the push that completed it) and settling it, by claiming 32 items
// and starting them, or, once items run out, by popping whatever a nonempty ring holds (at most as many as
// there are idle lanes, so that busy warps do not hoard the last orbits), found through the rings' nonempty
// mask.
//
// Finishes: tasks with a Record store compact records, and the warp finishes 32 at once, instead of a few lanes
// per burst with the rest idle (only if the buffer fits in static shared memory alongside the queue).
//
// Termination: counters[2] counts active warps.  An idle warp leaves the count, and rejoins it before popping a
// ring it saw nonempty; only active warps push.  So once the count is zero with every ring empty, nothing more
// can arrive, and the warp exits.  Once items are out and at most drain_warps warps are active, the drain flag
// sends every warp's orbits (running and queued) to dump and the warp home, for the host to finish on CPU
// threads with whatever the rings still hold.
template<class Task, bool timing, int min_blocks> __global__ void __launch_bounds__(256, min_blocks)
pool_kernel(const Task task, const int64_t n, const int64_t stride, const Queues<Parked<typename Task::State>> rings,
            const IntRing ripe, Parked<typename Task::State>* dump, uint64_t* counters, const int drain_warps) {
  typedef typename Task::State State;
  typedef typename RecordOf<Task>::type Record;
  constexpr bool keyed = requires(State& s) { task.settle_key(s); };
  constexpr bool buffered = requires { typename Task::Record; } &&
                            256 * (sizeof(State) + sizeof(Record) + 8) <= 48 * 1024;
  __shared__ alignas(16) unsigned char queue_bytes[256 * sizeof(State)];
  __shared__ int32_t queue_items[256];
  __shared__ alignas(16) unsigned char record_bytes[buffered ? 256 * sizeof(Record) : 16];
  __shared__ int32_t record_items[buffered ? 256 : 1];
  const int lane = int(threadIdx.x & 31);
  const unsigned lanes_below = (1u << lane) - 1;
  WarpQueue<State> queue{reinterpret_cast<State*>(queue_bytes) + (threadIdx.x & ~31u),
                         queue_items + (threadIdx.x & ~31u)};
  [[maybe_unused]] Record* const records = reinterpret_cast<Record*>(record_bytes) + (buffered ? threadIdx.x & ~31u : 0);
  [[maybe_unused]] int32_t* const ritems = record_items + (buffered ? threadIdx.x & ~31u : 0);
  [[maybe_unused]] int records_n = 0;  // Buffered records (warp-uniform)
  State o;
  int32_t i = -1;  // Current item, or -1
  bool done = true, items_out = false, active = true;  // items_out and active are warp-uniform
  uint64_t it = 0, pushes = 0, alone = 0;
  uint64_t run_cycles = 0, lane_steps = 0, warp_steps = 0, settle_cycles = 0, idle_cycles = 0, refill_cycles = 0,
           push_cycles = 0, finish_cycles = 0, key_cycles = 0, ripe_cycles = 0;
  const long long t0 = timing ? clock64() : 0;
  if (!lane) atomic_add(counters + 2, 1);

  // Finish one orbit (each lane with has), through the record buffer if any
  const auto finish = [&](const bool has, const State& s, const int32_t item) {
    if constexpr (buffered) {
      const unsigned m = __ballot_sync(0xffffffff, has);
      if (!m) return;
      if (records_n + __popc(m) > 32) {
        if (lane < records_n) task.finish_record(records[lane], ritems[lane]);
        records_n = 0;
        __syncwarp();
      }
      if (has) {
        const int r = records_n + __popc(m & lanes_below);
        records[r] = task.record(s);
        ritems[r] = item;
      }
      records_n += __popc(m);
      __syncwarp();
    } else {
      if (has) task.finish(s, item);
    }
    if (has) it += task.iters(s);
  };
  const auto flush = [&]() {
    if constexpr (buffered) {
      if (lane < records_n) task.finish_record(records[lane], ritems[lane]);
      records_n = 0;
      __syncwarp();
    }
  };

  // Restock the queue: settle popped orbits, else start new items (see above).  Full batches come from the ripe
  // queue (keys whose rings gained 32 more orbits); only once items run out does a warp scan every ring.
  const auto restock = [&](const unsigned need) {
    State& q = queue.slots[lane];
    int32_t key = -1;
    if (!lane && !ripe.pop(key)) key = -1;
    key = __shfl_sync(0xffffffff, key, 0);
    if (key < 0 && items_out) key = rings.any();
    if (key >= 0) {
      Parked<State> p;
      const int got = rings.pop(key, items_out ? min(__popc(need), 32) : 32, p);
      if (got) {
        const long long c0 = timing ? clock64() : 0;
        bool fin = false, ok = false;
        if (lane < got) {
          q = p.o;
          fin = task.settle(q);
          if (!fin) { queue.items[lane] = p.item; ok = true; }
        }
        finish(fin, q, p.item);
        queue.ready = __ballot_sync(0xffffffff, ok);
        if constexpr (timing) settle_cycles += clock64() - c0;
        return true;
      }
    }
    if (items_out) return false;
    uint64_t j0 = 0;
    if (!lane) j0 = atomic_fetch_add(counters, uint64_t(32));
    j0 = __shfl_sync(0xffffffff, j0, 0);
    const int64_t j = int64_t(j0) + lane;
    bool started = false, decided = false;
    int32_t k = -1;
    if (j < n) {
      // Started in place, so that no second state occupies registers
      k = int32_t(scramble(j, stride, n));
      decided = task.start(q, k);
      if (!decided) { queue.items[lane] = k; started = true; }
    }
    finish(decided, q, k);  // Decided at once (the cardioid)
    items_out = int64_t(j0) + 32 >= n;
    queue.ready = __ballot_sync(0xffffffff, started);
    return int64_t(j0) < n;
  };
  // A counter as the whole warp sees it (lane 0's read)
  const auto uniform = [&](const uint64_t* c) {
    uint64_t v = 0;
    if (!lane) v = ring_load(c);
    return __shfl_sync(0xffffffff, v, 0);
  };

  for (;;) {
    // Drain: hand this warp's orbits to the host
    if (items_out && uniform(counters + 3)) {
      const bool run_mine = !done && i >= 0;
      const unsigned queued = queue.ready;
      const int mine = int(run_mine) + int((queued >> lane) & 1);
      const unsigned any = __ballot_sync(0xffffffff, mine > 0);
      if (any) {
        // Warp-aggregated reservation of this warp's dump slots
        int total = 0, before = 0;
        for (int l = 0; l < 32; l++) {
          const int c = __shfl_sync(0xffffffff, mine, l);
          if (l < lane) before += c;
          total += c;
        }
        uint64_t base = 0;
        if (!lane) base = atomic_fetch_add(counters + 6, uint64_t(total));
        base = __shfl_sync(0xffffffff, base, 0);
        int at = int(base) + before;
        if (run_mine) { dump[at].o = o; dump[at].item = i; at++; }
        if ((queued >> lane) & 1) { dump[at].o = queue.slots[lane]; dump[at].item = queue.items[lane]; }
      }
      break;
    }
    const long long f0 = timing ? clock64() : 0;
    queue.refill(done, o, i, restock);
    if constexpr (timing) refill_cycles += clock64() - f0;
    if (__all_sync(0xffffffff, done)) {
      // Idle: out of items with the rings empty for now.  Leave the active count, and wait until a ring has
      // work (rejoining before popping it), everything is done, or the drain begins.
      const long long c0 = timing ? clock64() : 0;
      if (active && !lane) atomic_add(counters + 2, uint64_t(-1));
      active = false;
      bool leave = false;
      for (unsigned sleep = 1000;; sleep = min(2 * sleep, 64000u)) {
        if (uniform(counters + 3)) break;
        if (rings.any() >= 0) {
          if (!lane) atomic_add(counters + 2, 1);
          active = true;
          break;
        }
        const uint64_t a = uniform(counters + 2);
        if (!a && rings.any() < 0) { leave = true; break; }
        if (a <= uint64_t(drain_warps) && !lane) atomicExch(reinterpret_cast<unsigned long long*>(counters + 3), 1ull);
        __nanosleep(sleep);  // Backing off, so that idle warps' scans do not load the rings the rest are using
      }
      if constexpr (timing) idle_cycles += clock64() - c0;
      if (leave) break;
      continue;
    }
    if (!done) {
      const long long r0 = timing ? clock64() : 0;
      const int64_t p0 = timing ? task.progress(o) : 0;
      done = task.run(o);
      // Settles due at once (an overflowed block), for all such lanes of the warp together
      if constexpr (requires { task.immediate(o); }) if (done && task.immediate(o)) done = task.settle(o);
      if constexpr (timing) {
        run_cycles += clock64() - r0;
        const unsigned mask = __activemask(), steps = unsigned(task.progress(o) - p0),
                       longest = __reduce_max_sync(mask, steps);
        lane_steps += steps;
        if ((threadIdx.x & 31) == __ffs(mask) - 1) warp_steps += uint64_t(32) * longest;
      }
    }
    __syncwarp();
    const long long u0 = timing ? clock64() : 0;
    // Pending orbits go onto the ring for their settle key; if it is full, the lane settles alone
    const bool pend = !done ? false : i >= 0 && task.pending(o);
    int key = 0;
    if constexpr (keyed) if (pend) key = task.settle_key(o);
    __syncwarp();
    const long long k1 = timing ? clock64() : 0;
    Parked<State> p;
    if (pend) { p.o = o; p.item = i; }
    int batches = 0;
    const bool pushed = rings.push(pend, key, p, batches);
    const long long k2 = timing ? clock64() : 0;
    // Announce each full batch on the ripe queue (dropping announcements if it is half full)
    for (; batches > 0; batches--) ripe.push_if_room(key);
    if constexpr (timing) { key_cycles += k1 - u0; ripe_cycles += clock64() - k2; }
    if (pend) {
      if (pushed) { i = -1; pushes++; }
      else { done = task.settle(o); alone++; }
    }
    const long long u1 = timing ? clock64() : 0;
    // Finished orbits
    const bool fin = done && i >= 0;
    finish(fin, o, i);
    if (fin) i = -1;
    if constexpr (timing) { push_cycles += u1 - u0; finish_cycles += clock64() - u1; }
  }
  flush();
  atomic_add(counters + 1, it);
  atomic_add(counters + 4, pushes);
  atomic_add(counters + 5, alone);
  if constexpr (timing) {
    atomic_add(counters + 7, run_cycles);
    atomic_add(counters + 8, uint64_t(clock64() - t0));
    atomic_add(counters + 9, lane_steps);
    atomic_add(counters + 10, warp_steps);
    atomic_add(counters + 11, settle_cycles);
    atomic_add(counters + 12, idle_cycles);
    atomic_add(counters + 13, refill_cycles);
    atomic_add(counters + 14, push_cycles);
    atomic_add(counters + 15, finish_cycles);
    atomic_add(counters + 16, key_cycles);
    atomic_add(counters + 17, ripe_cycles);
  }
}

// Launch pool_kernel with register pressure chosen at run time (min_blocks resident 256-thread blocks per SM),
// on as many blocks as fit at once: the kernel's termination counts every started warp as active until it idles,
// so blocks that could only start once others exit would do nothing
template<class Task, bool timing> int launch_pool_kernel(const int min_blocks, const Task& task, const int64_t n,
                                                         const int64_t stride,
                                                         const Queues<Parked<typename Task::State>>& rings,
                                                         const IntRing& ripe,
                                                         Parked<typename Task::State>* dump, uint64_t* counters,
                                                         const int drain_warps, const bool launch) {
  int per_sm = 0;
#define LAUNCH(b) { \
    const auto kernel = pool_kernel<Task, timing, b>; \
    cuda_check(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&per_sm, kernel, 256, 0)); \
    if (launch) \
      kernel<<<per_sm * num_sms(), 256, 0, stream()>>>(task, n, stride, rings, ripe, dump, counters, drain_warps); }
  switch (min_blocks) {
    case 1: LAUNCH(1); break;
    case 2: LAUNCH(2); break;
    case 3: LAUNCH(3); break;
    case 4: LAUNCH(4); break;
    default: die("MANDELBROT_CUDA_MIN_BLOCKS must be 1 to 4, got %d", min_blocks);
  }
#undef LAUNCH
  return per_sm * num_sms();
}

// Every ring slot free for its first push
struct InitSeq {
  uint64_t* seq;
  uint64_t cap;
  __device__ void operator()(const int64_t i) const { seq[i] = uint64_t(i) & (cap - 1); }
};

// Finish orbits whose states were completed on the host (CPU tail)
template<class Task> __global__ void finish_kernel(const Task task, const Parked<typename Task::State>* orbits,
                                                   const int64_t count) {
  for (int64_t k = int64_t(blockIdx.x) * blockDim.x + threadIdx.x; k < count; k += int64_t(blockDim.x) * gridDim.x)
    task.finish(orbits[k].o, orbits[k].item);
}

template<class F> __global__ void for_each_kernel(const int64_t n, const F f);

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
    typedef engine_detail::Parked<O> P;
    static const int min_blocks_env = env_int("MANDELBROT_CUDA_MIN_BLOCKS", 0),  // Override task.min_blocks
                     timing = env_int("MANDELBROT_CUDA_TIMING", 0),
                     // Finish on CPU threads once this few orbits remain (in about 1/32 as many warps): a lone GPU
                     // lane steps a sequential orbit ~40× slower than a CPU core
                     cpu_tail = env_int("MANDELBROT_CUDA_CPU_TAIL", 1024);
    constexpr bool keyed = requires(O& o) { task.settle_key(o); };
    const int keys = keyed ? kSettleKeys : 1, min_blocks = min_blocks_env ? min_blocks_env : task.min_blocks;
    // Settle queues and the ripe queue (about one key per 32 queued), allocated once per thread and reset per run
    static thread_local std::unique_ptr<QueuesMem<P>> queues_cache[2];
    static thread_local std::unique_ptr<IntRingMem> ripe_cache;
    auto& queues_ptr = queues_cache[keyed];
    if (!queues_ptr) queues_ptr.reset(new QueuesMem<P>(keys, kSettleBound));
    else queues_ptr->reset();
    if (!ripe_cache) ripe_cache.reset(new IntRingMem(next_pow2(2 * uint64_t(kSettleBound) / 32 + 4096)));
    ripe_cache->reset();
    QueuesMem<P>& queues = *queues_ptr;
    IntRingMem& ripe_mem = *ripe_cache;
    Mem<uint64_t> counters(engine_detail::kCounters, true);
    counters.zero();
    auto rings = queues.queues();
    if (env_int("MANDELBROT_NO_SPACE", 0)) rings.space = nullptr;  // Experiment: no bound on outstanding orbits
    const auto ripe = ripe_mem.ring();
    const int drain_warps = cpu_tail > 0 ? std::max(1, cpu_tail / 32) : -1;
    const int grid = timing ? engine_detail::launch_pool_kernel<Task, true>(min_blocks, task, n, stride, rings,
                                                                            ripe, nullptr, counters.p, drain_warps, false)
                            : engine_detail::launch_pool_kernel<Task, false>(min_blocks, task, n, stride, rings,
                                                                             ripe, nullptr, counters.p, drain_warps, false);
    Mem<P> dump(int64_t(grid) * 256 * 2, true);  // Each lane's running and queued orbits
    cudaEvent_t e0, e1;
    cuda_check(cudaEventCreate(&e0)); cuda_check(cudaEventCreate(&e1));
    cuda_check(cudaEventRecord(e0, stream()));
    if (timing)
      engine_detail::launch_pool_kernel<Task, true>(min_blocks, task, n, stride, rings, ripe, dump.p, counters.p,
                                                    drain_warps, true);
    else
      engine_detail::launch_pool_kernel<Task, false>(min_blocks, task, n, stride, rings, ripe, dump.p, counters.p,
                                                     drain_warps, true);
    cuda_check(cudaGetLastError());
    cuda_check(cudaEventRecord(e1, stream()));
    uint64_t h[engine_detail::kCounters];
    counters.to_host(h, engine_detail::kCounters);
    stats.overflow = int64_t(h[4]);
    int64_t cpu_iters = 0;
    if (h[3]) {
      // Drained: finish the dumped orbits and those left in the rings on CPU threads, then write their results
      std::vector<P> left(h[6]);
      dump.to_host(left.data(), int64_t(h[6]));
      for (const auto& p : queues.contents()) left.push_back(p);
      const int64_t count = int64_t(left.size());
      std::atomic<int64_t> next_k(0), iters(0);
      std::vector<std::thread> pool;
      for (int t = 0; t < cpu_threads(); t++)
        pool.emplace_back([&]() {
          int64_t it = 0;
          for (int64_t k; (k = next_k.fetch_add(1)) < count;) {
            O& o = left[k].o;
            for (bool done = false; !done;) {
              if (task.pending(o)) done = task.settle(o);
              else done = task.run(o) && !task.pending(o);
            }
            it += task.iters(o);
          }
          iters += it;
        });
      for (auto& t : pool) t.join();
      if (count) {
        Mem<P> back(count, true);
        back.from_host(left.data(), count);
        engine_detail::finish_kernel<Task><<<int(std::min<int64_t>(grid, (count + 255) / 256)), 256, 0,
                                             stream()>>>(task, back.p, count);
        cuda_check(cudaGetLastError());
        cuda_sync();
      }
      cpu_iters = iters;
      stats.cpu_tail = count;
    }
    stats.iters = int64_t(h[1]) + cpu_iters;
    if (timing) {
      float ms;
      cuda_check(cudaEventElapsedTime(&ms, e0, e1));
      const double cycles = double(h[8]);
      print("    cuda run (%s): %d items, %d threads, %.1f ms, %.3g it/s; %d ring pushes, %d settled alone, "
            "%d to the CPU tail", typeid(Task).name(), n, grid * 256, ms, double(h[1]) / (ms * 1e-3), h[4], h[5],
            stats.cpu_tail);
      print("      thread cycles: %.1f%% in run, %.1f%% refilling (%.1f%% settling popped orbits), %.1f%% pushing, "
            "%.1f%% finishing, %.1f%% idle; lane steps / warp steps %.1f%%", 100 * double(h[7]) / cycles,
            100 * double(h[13]) / cycles, 100 * double(h[11]) / cycles, 100 * double(h[14]) / cycles,
            100 * double(h[15]) / cycles, 100 * double(h[12]) / cycles, 100.0 * double(h[9]) / double(h[10]));
      print("      pushing: %.1f%% settle keys, %.1f%% ripe announcements", 100 * double(h[16]) / cycles,
            100 * double(h[17]) / cycles);
    }
    cuda_check(cudaEventDestroy(e0)); cuda_check(cudaEventDestroy(e1));
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
