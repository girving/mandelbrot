// Bounded multi-producer multi-consumer rings in device memory, one per key, with warp-cooperative batch pushes
// and pops
//
// The GPU engine keeps orbits waiting for a settle (Newton) in one ring per settle key, so that a warp can pop 32
// orbits with the same key and settle them together.  Every step is a fetch-and-add, never a CAS loop, since
// thousands of warps push onto a few popular keys at once, and a contended CAS succeeds for one of them per round
// trip.  A push reserves positions from tail, waits for its slots to be free, writes them, publishes them
// (Vyukov's per-slot sequence numbers: slot p mod cap holds p when free for the push at position p, p + 1 once
// written, and p + cap once popped), and adds them to avail.  A pop takes up to its want from avail (returning any
// overshoot), claims that many positions from head, and waits for their writes (already under way: avail never
// exceeds the reservations).  Space counts slots neither reserved nor still being popped: a push takes its
// slots from it first (reporting a full ring rather than waiting), and a pop returns them once read, so that a
// push waits only for pops already reading its slots, and a pop only for pushes already writing them.
//
// All device functions are warp-collective: every lane of a converged warp calls them.
#pragma once

#include "cutil.h"
#include <cstdint>
namespace mandelbrot {

template<class T> struct Rings {
  T* slots;         // [key * cap + position mod cap]
  uint64_t* seq;    // Sequence numbers, as slots
  uint64_t* head;   // [key]: next position to pop
  uint64_t* tail;   // [key]: next position to push
  int64_t* avail;   // [key]: published values not yet taken by a pop
  int64_t* space;   // [key]: free slots not yet reserved (initially cap)
  int keys;
  uint64_t cap;     // A power of 2
  uint32_t* nonempty = nullptr;  // Optional: bit k of word k / 32 set while ring k may hold values (see any)

#ifdef __CUDACC__
  // Shared state is read and written at L2 (ld/st.cg): SMs' L1 caches are not coherent, so a cached head, tail,
  // sequence number, or slot from an earlier lap could be stale indefinitely
  __device__ static uint64_t load(const uint64_t* p) {
    return __ldcg(reinterpret_cast<const unsigned long long*>(p));
  }
  __device__ static int64_t load(const int64_t* p) { return __ldcg(reinterpret_cast<const long long*>(p)); }
  __device__ static uint32_t load32(const uint32_t* p) { return __ldcg(reinterpret_cast<const unsigned*>(p)); }
  __device__ static void store(uint64_t* p, const uint64_t v) {
    asm volatile("st.global.cg.u64 [%0], %1;" :: "l"(p), "l"(v) : "memory");
  }
  __device__ static uint64_t fetch_add(uint64_t* p, const uint64_t v) {
    return atomicAdd(reinterpret_cast<unsigned long long*>(p), static_cast<unsigned long long>(v));
  }
  __device__ static int64_t fetch_add(int64_t* p, const int64_t v) {
    return int64_t(atomicAdd(reinterpret_cast<unsigned long long*>(p), static_cast<unsigned long long>(v)));
  }
  // Copy a slot's value at L2, in 8- or 4-byte words
  __device__ static void load_value(T& dst, const T* src) {
    if constexpr (sizeof(T) % 8 == 0) {
      const unsigned long long* s = reinterpret_cast<const unsigned long long*>(src);
      unsigned long long* d = reinterpret_cast<unsigned long long*>(&dst);
      for (int j = 0; j < int(sizeof(T) / 8); j++) d[j] = __ldcg(s + j);
    } else {
      static_assert(sizeof(T) % 4 == 0);
      const unsigned* s = reinterpret_cast<const unsigned*>(src);
      unsigned* d = reinterpret_cast<unsigned*>(&dst);
      for (int j = 0; j < int(sizeof(T) / 4); j++) d[j] = __ldcg(s + j);
    }
  }
  __device__ uint64_t slot(const int key, const uint64_t pos) const { return uint64_t(key) * cap + (pos & (cap - 1)); }

  // Lanes with has push v onto ring key.  Returns whether this lane's push succeeded (false if its ring was
  // nearly full).  ripe counts, in one lane per key pushed, how many multiples of 32 that key's avail crossed:
  // each marks 32 more values, a full batch for a warp to pop.
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
      // Reserve space, then positions
      uint64_t t = ~uint64_t(0);
      if (lane == lead) {
        if (fetch_add(space + k, -int64_t(count)) >= count) t = fetch_add(tail + k, uint64_t(count));
        else fetch_add(space + k, int64_t(count));  // Full
      }
      t = __shfl_sync(0xffffffff, t, lead);
      if (t != ~uint64_t(0)) {
        if (mine) {
          const uint64_t pos = t + rank, s = slot(k, pos);
          while (load(seq + s) != pos) {}  // Free once the previous lap's pop (under way) has read it
          slots[s] = v;
          __threadfence();
          store(seq + s, pos + 1);
          ok = true;
        }
        __syncwarp();
        if (lane == lead) {
          __threadfence();
          const int64_t a = fetch_add(avail + k, int64_t(count));
          ripe = int((a + count) / 32 - a / 32);
          if (nonempty) {
            uint32_t* w = nonempty + k / 32;
            const uint32_t bit = 1u << (k & 31);
            if (!(load32(w) & bit)) atomicOr(w, bit);
          }
        }
      }
      todo &= ~group;
    }
    return ok;
  }

  // Pop up to want ≤ 32 values from ring key (warp-uniform) into out, one per lane: lanes below the returned
  // count receive one
  __device__ int pop(const int key, const int want, T& out) const {
    const int lane = int(threadIdx.x & 31);
    int got = 0;
    uint64_t h = 0;
    if (!lane) {
      const int64_t a = load(avail + key);
      if (a > 0) {
        const int64_t take = a < want ? a : want, before = fetch_add(avail + key, -take);
        got = int(before >= take ? take : before > 0 ? before : 0);
        if (got < take) fetch_add(avail + key, take - got);  // Return the overshoot
        if (got) h = fetch_add(head + key, uint64_t(got));
        if (nonempty && before - take <= 0) {
          // Possibly emptied it: clear its bit, then look again, since a push may have landed in between
          uint32_t* w = nonempty + key / 32;
          const uint32_t bit = 1u << (key & 31);
          atomicAnd(w, ~bit);
          __threadfence();
          if (load(avail + key) > 0) atomicOr(w, bit);
        }
      }
    }
    got = __shfl_sync(0xffffffff, got, 0);
    h = __shfl_sync(0xffffffff, h, 0);
    if (lane < got) {
      const uint64_t pos = h + lane, s = slot(key, pos);
      while (load(seq + s) != pos + 1) {}  // Published, or about to be (avail counted it)
      __threadfence();
      load_value(out, slots + s);
      __threadfence();
      store(seq + s, pos + cap);
    }
    __syncwarp();
    if (!lane && got) fetch_add(space + key, int64_t(got));
    return got;
  }

  // Some key whose ring may hold values (with nonempty), else -1: one load per lane for up to 1024 keys
  __device__ int any() const {
    const int lane = int(threadIdx.x & 31);
    const uint32_t w = lane * 32 < keys ? load32(nonempty + lane) : 0;
    const unsigned has = __ballot_sync(0xffffffff, w != 0);
    if (!has) return -1;
    const int who = __ffs(has) - 1;
    return who * 32 + __ffs(__shfl_sync(0xffffffff, w, who)) - 1;
  }
#endif  // __CUDACC__
};

}  // namespace mandelbrot
