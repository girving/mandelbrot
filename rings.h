// Multi-producer multi-consumer queues in device memory, one per key, with warp-cooperative batch pushes and pops
//
// The GPU engine keeps orbits waiting for a settle (Newton) in one queue per settle key, so that a warp can pop 32
// orbits with the same key and settle them together.  Thousands of warps push onto a few popular keys at once,
// so every contended step is a fetch-and-add (a contended CAS succeeds for one warp per round trip).
//
// Queues: positions per key come from fetch-and-add on tail (pushes) and head (pops); position p of key k lives
// in slot p mod S of segment table[k][p / S mod entries], a block of S slots taken from a shared pool by the first
// push to reach it and returned once all S of its values are popped.  Slots are written once per segment lifetime, so
// a push never waits for a pop, and a pop waits only for pushes already under way.  A pop takes its count from
// avail (published values, returning any overshoot), then claims that many positions from head.  A semaphore
// (space) bounds the values outstanding, which bounds the segments and table entries in use, so neither runs
// out; a push finding no space reports failure instead.  (Pops claim the oldest positions, whose pushes may still
// be allocating, so the values outstanding can each hold a segment: the pool needs one segment per value it may
// hold, not one per S.)
//
// IntRing: a bounded MPMC ring of ints for at most cap outstanding values (the pool's free segments, and the
// queue of ripe keys), so that a push's slot is always free or being read.
#pragma once

#include "cutil.h"
#include <cstdint>
namespace mandelbrot {

#ifdef __CUDACC__
// Shared state is read and written with relaxed GPU-scope operations, as volatile asm: SMs' L1 caches are not
// coherent (a plain load can stay stale), and the compiler must not hoist loads out of spin loops
__device__ static inline uint64_t ring_load(const uint64_t* p) {
  uint64_t v;
  asm volatile("ld.relaxed.gpu.global.u64 %0, [%1];" : "=l"(v) : "l"(p) : "memory");
  return v;
}
__device__ static inline int64_t ring_load(const int64_t* p) {
  return int64_t(ring_load(reinterpret_cast<const uint64_t*>(p)));
}
__device__ static inline uint32_t ring_load(const uint32_t* p) {
  uint32_t v;
  asm volatile("ld.relaxed.gpu.global.u32 %0, [%1];" : "=r"(v) : "l"(p) : "memory");
  return v;
}
__device__ static inline int32_t ring_load(const int32_t* p) {
  return int32_t(ring_load(reinterpret_cast<const uint32_t*>(p)));
}
__device__ static inline void ring_store(uint64_t* p, const uint64_t v) {
  asm volatile("st.relaxed.gpu.global.u64 [%0], %1;" :: "l"(p), "l"(v) : "memory");
}
__device__ static inline void ring_store(int32_t* p, const int32_t v) {
  asm volatile("st.relaxed.gpu.global.u32 [%0], %1;" :: "l"(p), "r"(v) : "memory");
}
__device__ static inline uint64_t ring_add(uint64_t* p, const uint64_t v) {
  return atomicAdd(reinterpret_cast<unsigned long long*>(p), static_cast<unsigned long long>(v));
}
__device__ static inline int64_t ring_add(int64_t* p, const int64_t v) {
  return int64_t(atomicAdd(reinterpret_cast<unsigned long long*>(p), static_cast<unsigned long long>(v)));
}
// Copy a value written by another SM, in 8- or 4-byte words
template<class T> __device__ static inline void ring_copy(T& dst, const T* src) {
  if constexpr (sizeof(T) % 8 == 0) {
    const uint64_t* s = reinterpret_cast<const uint64_t*>(src);
    uint64_t* d = reinterpret_cast<uint64_t*>(&dst);
    for (int j = 0; j < int(sizeof(T) / 8); j++) d[j] = ring_load(s + j);
  } else {
    static_assert(sizeof(T) % 4 == 0);
    const uint32_t* s = reinterpret_cast<const uint32_t*>(src);
    uint32_t* d = reinterpret_cast<uint32_t*>(&dst);
    for (int j = 0; j < int(sizeof(T) / 4); j++) d[j] = ring_load(s + j);
  }
}

// Spin loops: with -DMANDELBROT_RING_DEBUG, give up after 2^20 iterations, counting where in ring_stuck
#ifdef MANDELBROT_RING_DEBUG
__device__ unsigned long long ring_stuck[8];
#define RING_SPIN(where, cond) \
  for (uint32_t spins_ = 0; (cond); spins_++) \
    if (spins_ == (1u << 20)) { atomicAdd(ring_stuck + (where), 1ull); break; }
#else
#define RING_SPIN(where, cond) while (cond) {}
#endif
#endif  // __CUDACC__

// Bounded ring of ints, lane by lane, for at most cap values outstanding.  Slot p mod cap holds sequence number p
// when free for the push at p, p + 1 once written, and p + cap once popped (Vyukov).
struct IntRing {
  int32_t* slots;
  uint64_t* seq;    // [cap], initially i
  uint64_t* ends;   // [head, tail]
  int64_t* avail;   // Published values not yet taken
  uint64_t cap;     // A power of 2

#ifdef __CUDACC__
  // Push v, unless the ring is half full (so that pushes in flight never exceed cap): returns whether it did
  __device__ bool push_if_room(const int32_t v) const {
    if (ring_load(ends + 1) - ring_load(ends) >= cap / 2) return false;
    push(v);
    return true;
  }
  __device__ void push(const int32_t v) const {
    const uint64_t t = ring_add(ends + 1, uint64_t(1)), s = t & (cap - 1);
    RING_SPIN(0, ring_load(seq + s) != t);  // At most cap outstanding: free, or its pop is reading it
    ring_store(slots + s, v);
    __threadfence();
    ring_store(seq + s, t + 1);
    __threadfence();
    ring_add(avail, int64_t(1));
  }
  __device__ bool pop(int32_t& v) const {
    if (ring_load(avail) <= 0) return false;
    if (ring_add(avail, int64_t(-1)) <= 0) {
      ring_add(avail, int64_t(1));
      return false;
    }
    const uint64_t h = ring_add(ends, uint64_t(1)), s = h & (cap - 1);
    RING_SPIN(1, ring_load(seq + s) != h + 1);  // Its push is writing it
    __threadfence();
    v = ring_load(slots + s);
    __threadfence();
    ring_store(seq + s, h + cap);
    return true;
  }
#endif  // __CUDACC__
};

template<class T> struct Queues {
  static constexpr int S = 16;  // Slots per segment
  T* values;            // [segments * S]
  uint64_t* vseq;       // [segments * S]: position + 1 once written, 0 once popped (and initially)
  uint64_t* table;      // [keys * T]: (block << 32) | segment, or ~0 if none
  uint32_t* consumed;   // [segments]: values popped from each segment in use
  uint64_t* head;       // [keys]
  uint64_t* tail;       // [keys]
  int64_t* avail;       // [keys]: published values not yet taken by a pop
  int64_t* space;       // Values that may still be pushed (initially the bound on outstanding values)
  uint32_t* nonempty;   // [keys / 32]: bit k set while queue k may hold values (see any)
  IntRing free;         // Free segments
  int keys;
  uint64_t entries;     // Table entries per key, a power of 2

#ifdef __CUDACC__
  // The segment holding position p of key k, allocated by the first push to reach it
  __device__ int32_t segment(const int k, const uint64_t p, const bool push) const {
    const uint64_t block = p / S, tag = block << 32;
    uint64_t* e = table + uint64_t(k) * entries + (block & (entries - 1));
#ifdef MANDELBROT_RING_DEBUG
    for (uint32_t spins = 0;; spins++) {
      if (spins == (1u << 20)) { atomicAdd(ring_stuck + (push ? 2 : 3), 1ull); return 0; }
#else
    for (;;) {
#endif
      const uint64_t v = ring_load(e);
      if (v != ~uint64_t(0) && (v & ~uint64_t(0xffffffff)) == tag) return int32_t(v & 0xffffffff);
      if (!push || v != ~uint64_t(0)) continue;  // Not yet allocated (pops wait for its push)
      int32_t id;
      RING_SPIN(4, !free.pop(id));  // The space bound leaves segments to spare
      const uint64_t old = atomicCAS(reinterpret_cast<unsigned long long*>(e), ~0ull,
                                     static_cast<unsigned long long>(tag | uint32_t(id)));
      if (old == ~uint64_t(0)) return id;
      free.push(id);  // Another push allocated it first
    }
  }

  // Lanes with has push v onto queue key.  Returns whether this lane's push succeeded (false if the bound on
  // outstanding values is reached).  ripe counts, in one lane per key pushed, how many multiples of 32 that key's
  // avail crossed: each marks 32 more values, a full batch for a warp to pop.
  __device__ bool push(const bool has, const int key, const T& v, int& ripe) const {
    const int lane = int(threadIdx.x & 31);
    bool ok = false;
    ripe = 0;
    for (unsigned todo = __ballot_sync(0xffffffff, has); todo;) {
      // Lanes pushing the same key as the lowest remaining one, as a group
      const int lead = __ffs(todo) - 1, k = __shfl_sync(0xffffffff, key, lead);
      const unsigned group = __ballot_sync(0xffffffff, has && key == k) & todo;
      const bool mine = (group >> lane) & 1;
      const int count = __popc(group), rank = __popc(group & ((1u << lane) - 1));
      uint64_t t = ~uint64_t(0);
      if (lane == lead) {
        if (ring_add(space, -int64_t(count)) >= count) t = ring_add(tail + k, uint64_t(count));
        else ring_add(space, int64_t(count));  // No space: give it back
      }
      t = __shfl_sync(0xffffffff, t, lead);
      if (t != ~uint64_t(0)) {
        if (mine) {
          const uint64_t p = t + rank, s = uint64_t(segment(k, p, true)) * S + (p & (S - 1));
          values[s] = v;
          __threadfence();
          ring_store(vseq + s, p + 1);
          ok = true;
        }
        __syncwarp();
        if (lane == lead) {
          __threadfence();
          const int64_t a = ring_add(avail + k, int64_t(count));
          ripe = int((a + count) / 32 - a / 32);
          uint32_t* w = nonempty + k / 32;
          const uint32_t bit = 1u << (k & 31);
          if (!(ring_load(w) & bit)) atomicOr(w, bit);
        }
      }
      todo &= ~group;
    }
    return ok;
  }

  // Pop up to want ≤ 32 values from queue key (warp-uniform) into out, one per lane: lanes below the returned
  // count receive one
  __device__ int pop(const int key, const int want, T& out) const {
    const int lane = int(threadIdx.x & 31);
    int got = 0;
    uint64_t h = 0;
    if (!lane) {
      const int64_t a = ring_load(avail + key);
      int64_t left = 0;
      if (a > 0) {
        const int64_t take = a < want ? a : want, before = ring_add(avail + key, -take);
        got = int(before >= take ? take : before > 0 ? before : 0);
        if (got < take) ring_add(avail + key, take - got);  // Return the overshoot
        if (got) h = ring_add(head + key, uint64_t(got));
        left = before - take;
      }
      if (left <= 0) {
        // Possibly empty: clear its bit, then look again, since a push may have landed in between.  (A push sets
        // the bit after adding to avail, so a pop taking that value first leaves it set: the next pop to find
        // the queue empty clears it here.)
        uint32_t* w = nonempty + key / 32;
        const uint32_t bit = 1u << (key & 31);
        if (ring_load(w) & bit) {
          atomicAnd(w, ~bit);
          __threadfence();
          if (ring_load(avail + key) > 0) atomicOr(w, bit);
        }
      }
    }
    got = __shfl_sync(0xffffffff, got, 0);
    h = __shfl_sync(0xffffffff, h, 0);
    if (lane < got) {
      const uint64_t p = h + lane;
      const int32_t seg = segment(key, p, false);
      const uint64_t s = uint64_t(seg) * S + (p & (S - 1));
      RING_SPIN(5, ring_load(vseq + s) != p + 1);  // Published, or its push is writing it (avail counted it)
      __threadfence();
      ring_copy(out, values + s);
      ring_store(vseq + s, 0);
      __threadfence();
      if (atomicAdd(consumed + seg, 1u) == S - 1) {
        // The segment's last value: return it to the pool
        consumed[seg] = 0;
        uint64_t* e = table + uint64_t(key) * entries + ((p / S) & (entries - 1));
        atomicExch(reinterpret_cast<unsigned long long*>(e), ~0ull);
        __threadfence();
        free.push(seg);
      }
    }
    __syncwarp();
    if (!lane && got) ring_add(space, int64_t(got));
    return got;
  }

  // Some key whose queue may hold values, else -1: one load per lane for up to 1024 keys
  __device__ int any() const {
    const int lane = int(threadIdx.x & 31);
    const uint32_t w = lane * 32 < keys ? ring_load(nonempty + lane) : 0;
    const unsigned has = __ballot_sync(0xffffffff, w != 0);
    if (!has) return -1;
    const int who = __ffs(has) - 1;
    return who * 32 + __ffs(__shfl_sync(0xffffffff, w, who)) - 1;
  }
#endif  // __CUDACC__
};

}  // namespace mandelbrot
