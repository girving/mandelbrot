// Bounded multi-producer multi-consumer rings in device memory, one per key, with warp-cooperative batch pushes
// and pops
//
// The GPU engine keeps orbits waiting for a settle (Newton) in one ring per settle key, so that a warp can pop 32
// orbits with the same key and settle them together.  Each slot carries a sequence number (Vyukov's bounded MPMC
// queue): slot p mod cap holds p when free for the push at position p, p + 1 once that push has written it, and
// p + cap once the pop at position p has read it (freeing it for the push at p + cap).  Pushes and pops claim
// consecutive positions with one CAS per warp and key, then write or read their slots and publish them.
//
// All device functions are warp-collective: every lane of a converged warp calls them.
#pragma once

#include "cutil.h"
#include <cstdint>
namespace mandelbrot {

template<class T> struct Rings {
  T* slots;        // [key * cap + position mod cap]
  uint64_t* seq;   // Sequence numbers, as slots
  uint64_t* head;  // [key]: next position to pop
  uint64_t* tail;  // [key]: next position to push
  int keys;
  uint64_t cap;    // A power of 2
  uint32_t* nonempty = nullptr;  // Optional: bit k of word k / 32 set while ring k may hold values (see any)

#ifdef __CUDACC__
  __device__ static uint64_t load(const uint64_t* p) { return *reinterpret_cast<const volatile uint64_t*>(p); }
  __device__ static void store(uint64_t* p, const uint64_t v) { *reinterpret_cast<volatile uint64_t*>(p) = v; }
  __device__ uint64_t slot(const int key, const uint64_t pos) const { return uint64_t(key) * cap + (pos & (cap - 1)); }
  __device__ static bool cas(uint64_t* p, const uint64_t old, const uint64_t v) {
    return atomicCAS(reinterpret_cast<unsigned long long*>(p), static_cast<unsigned long long>(old),
                     static_cast<unsigned long long>(v)) == old;
  }

  // Lanes with has push v onto ring key.  Returns whether this lane's push succeeded (false if its ring was full).
  // ripe counts, in one lane per key pushed, how many multiples of 32 that key's tail crossed: each marks 32 more
  // values, a full batch for a warp to pop.
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
      for (;;) {
        uint64_t t = 0;
        if (lane == lead) t = load(tail + k);
        t = __shfl_sync(0xffffffff, t, lead);
        const uint64_t pos = t + rank;
        // Slot states: free for pos, still occupied from the previous lap (full), or already claimed (stale t)
        int state = 0;
        if (mine) {
          const uint64_t s = load(seq + slot(k, pos));
          state = s == pos ? 0 : s < pos ? 1 : 2;
        }
        if (__any_sync(0xffffffff, state == 1)) break;  // Full
        if (__any_sync(0xffffffff, state == 2)) continue;
        bool won = false;
        if (lane == lead) won = cas(tail + k, t, t + count);
        if (__shfl_sync(0xffffffff, won, lead)) {
          if (lane == lead) ripe = int((t + count) / 32 - t / 32);
          if (mine) {
            slots[slot(k, pos)] = v;
            __threadfence();
            store(seq + slot(k, pos), pos + 1);
            ok = true;
          }
          if (nonempty && lane == lead) {
            uint32_t* w = nonempty + k / 32;
            const uint32_t bit = 1u << (k & 31);
            if (!(*reinterpret_cast<volatile uint32_t*>(w) & bit)) atomicOr(w, bit);
          }
          break;
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
    for (;;) {
      uint64_t h = 0;
      if (!lane) h = load(head + key);
      h = __shfl_sync(0xffffffff, h, 0);
      const uint64_t pos = h + lane;
      // Slot states: written for pos (ready), not yet written (empty), or already popped (stale h)
      int state = 1;
      if (lane < want) {
        const uint64_t s = load(seq + slot(key, pos));
        state = s == pos + 1 ? 0 : s <= pos ? 1 : 2;
      }
      const unsigned ready = __ballot_sync(0xffffffff, state == 0);
      const int got = ready == 0xffffffff ? 32 : __ffs(~ready) - 1;  // Ready slots from the head on
      if (!got) {
        if (__shfl_sync(0xffffffff, state, 0) == 2) continue;
        return 0;
      }
      bool won = false;
      if (!lane) won = cas(head + key, h, h + got);
      if (!__shfl_sync(0xffffffff, won, 0)) continue;
      if (lane < got) {
        __threadfence();
        out = slots[slot(key, pos)];
        __threadfence();
        store(seq + slot(key, pos), pos + cap);
      }
      if (nonempty && !lane && load(tail + key) == h + got) {
        // Emptied it: clear its bit, then look again, since a push may have landed in between
        uint32_t* w = nonempty + key / 32;
        const uint32_t bit = 1u << (key & 31);
        atomicAnd(w, ~bit);
        __threadfence();
        if (load(tail + key) != load(head + key)) atomicOr(w, bit);
      }
      return got;
    }
  }

  // Some key whose ring may hold values (with nonempty), else -1: one load per lane for up to 1024 keys
  __device__ int any() const {
    const int lane = int(threadIdx.x & 31);
    const uint32_t w = lane * 32 < keys ? *reinterpret_cast<const volatile uint32_t*>(nonempty + lane) : 0;
    const unsigned has = __ballot_sync(0xffffffff, w != 0);
    if (!has) return -1;
    const int who = __ffs(has) - 1;
    return who * 32 + __ffs(__shfl_sync(0xffffffff, w, who)) - 1;
  }

#endif  // __CUDACC__
};

}  // namespace mandelbrot
