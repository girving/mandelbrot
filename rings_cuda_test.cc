// Stress test of the device settle queues: every warp pushes values onto queues chosen by a hash and pops from
// queues in turn, concurrently, with a bound on outstanding values small enough that pushes often fail and
// segments recycle many times.  Every value must be popped exactly once.
#define MANDELBROT_RING_DEBUG

#include "engine.h"
#include "rings.h"
#include "tests.h"
namespace mandelbrot {
namespace {

__global__ void stress(const Queues<uint32_t> rings, const int per_lane, uint32_t* seen, uint64_t* popped,
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
    if (!lane) p = ring_load(popped);
    if (__shfl_sync(0xffffffff, p, 0) >= total && __all_sync(0xffffffff, next == per_lane)) break;
    if (step > (1u << 20)) {  // Watchdog: report instead of hanging
      if (!lane) atomicAdd(reinterpret_cast<unsigned long long*>(full + 1), 1ull);
      break;
    }
  }
}

TEST(rings) {
  // Every block must be resident at once, since warps spin until every value is popped
  int per_sm = 0;
  cuda_check(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&per_sm, stress, 256, 0));
  ASSERT_LE(1, per_sm);
  const int keys = 7, blocks = num_sms(), per_lane = 32;
  print("rings: %d blocks", blocks);
  const int64_t bound = 4096;
  const uint64_t total = uint64_t(blocks) * 256 * per_lane;
  QueuesMem<uint32_t> queues(keys, bound);
  Mem<uint32_t> seen(total, true);
  Mem<uint64_t> counts(3, true);
  seen.zero();
  counts.zero();
  stress<<<blocks, 256, 0, stream()>>>(queues.queues(), per_lane, seen.p, counts.p, total, counts.p + 1);
  cuda_check(cudaGetLastError());
  std::vector<uint32_t> h(total);
  seen.to_host(h.data(), int64_t(total));
  uint64_t c[3];
  counts.to_host(c, 3);
  std::vector<uint64_t> e(2 * keys);
  queues.ends.to_host(e.data(), 2 * keys);
  if (c[2]) {
    std::vector<int64_t> av(keys);
    queues.avail.to_host(av.data(), keys);
    print("watchdog: %d warps stopped, %d of %d popped", c[2], c[0], total);
    for (int k = 0; k < keys; k++) print("  key %d: head %d, tail %d, avail %d", k, e[k], e[keys + k], av[k]);
  }
  unsigned long long stuck[8];
  cuda_check(cudaMemcpyFromSymbol(stuck, ring_stuck, sizeof(stuck)));
  for (int j = 0; j < 8; j++) if (stuck[j]) print("stuck in spin %d: %d lanes", j, stuck[j]);
  ASSERT_EQ(c[0], total);
  ASSERT_LT(uint64_t(0), c[1]) << "pushes never failed";
  for (uint64_t v = 0; v < total; v++) ASSERT_EQ(h[v], 1u) << tfm::format("value %d", v);
  // Every queue empty, with segments recycled many times
  for (int k = 0; k < keys; k++) {
    ASSERT_EQ(e[k], e[keys + k]);
    ASSERT_LT(uint64_t(10 * bound), e[k]);
  }
  ASSERT_TRUE(queues.contents().empty());
}

}  // namespace
}  // namespace mandelbrot
