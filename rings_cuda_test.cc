// Stress test of the device settle rings: every warp pushes values onto rings chosen by a hash and pops from
// rings in turn, concurrently, with rings small enough to fill and wrap many times.  Every value must be popped
// exactly once.

#include "engine.h"
#include "rings.h"
#include "tests.h"
namespace mandelbrot {
namespace {

__global__ void stress(const Rings<uint32_t> rings, const int per_lane, uint32_t* seen, uint64_t* popped,
                       const uint64_t total, uint64_t* full) {
  const int lane = int(threadIdx.x & 31);
  const uint32_t gid = blockIdx.x * blockDim.x + threadIdx.x, warp = gid / 32;
  int next = 0;
  for (uint32_t step = 0;; step++) {
    const bool has = next < per_lane;
    const uint32_t v = gid * uint32_t(per_lane) + uint32_t(next);
    const int key = int(mix64(v) % uint64_t(rings.keys));
    int ripe;
    const bool ok = rings.push(has, key, v, ripe);
    if (ok) next++;
    if (has && !ok) atomicAdd(reinterpret_cast<unsigned long long*>(full), 1ull);
    uint32_t w = 0;
    const int got = rings.pop(int((warp * 7 + step) % uint32_t(rings.keys)), 1 + int((step + warp) % 32), w);
    if (lane < got) atomicAdd(seen + w, 1u);
    if (!lane && got) atomicAdd(reinterpret_cast<unsigned long long*>(popped), static_cast<unsigned long long>(got));
    uint64_t p = 0;
    if (!lane) p = Rings<uint32_t>::load(popped);
    if (__shfl_sync(0xffffffff, p, 0) >= total && __all_sync(0xffffffff, next == per_lane)) break;
  }
}

struct Init {
  uint64_t* seq;
  uint64_t cap;
  __device__ void operator()(const int64_t i) const { seq[i] = uint64_t(i) & (cap - 1); }
};

TEST(rings) {
  const int keys = 7, blocks = 264, per_lane = 64;
  const uint64_t cap = 64;
  const uint64_t total = uint64_t(blocks) * 256 * per_lane;
  Mem<uint32_t> slots(keys * cap, true), seen(total, true);
  Mem<uint64_t> seq(keys * cap, true), ends(2 * keys, true), counts(2, true);
  ends.zero();
  seen.zero();
  counts.zero();
  for_each(keys * int64_t(cap), Init{seq.p, cap}, true);
  const Rings<uint32_t> rings{slots.p, seq.p, ends.p, ends.p + keys, keys, cap};
  stress<<<blocks, 256, 0, stream()>>>(rings, per_lane, seen.p, counts.p, total, counts.p + 1);
  cuda_check(cudaGetLastError());
  std::vector<uint32_t> h(total);
  seen.to_host(h.data(), int64_t(total));
  uint64_t c[2];
  counts.to_host(c, 2);
  ASSERT_EQ(c[0], total);
  ASSERT_LT(uint64_t(0), c[1]) << "rings never filled";
  for (uint64_t v = 0; v < total; v++) ASSERT_EQ(h[v], 1u) << tfm::format("value %d", v);
  // Every ring empty, with head and tail having wrapped many laps
  std::vector<uint64_t> e(2 * keys);
  ends.to_host(e.data(), 2 * keys);
  for (int k = 0; k < keys; k++) {
    ASSERT_EQ(e[k], e[keys + k]);
    ASSERT_LT(10 * cap, e[k]);
  }
}

}  // namespace
}  // namespace mandelbrot
