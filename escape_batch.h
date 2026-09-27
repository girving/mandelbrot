// Batched escape-time sampling of leaf cells, on CPU threads or the GPU
//
// Sample i lies in leaf i / m, in stratum (i % m) % strata^2 of the leaf's strata × strata grid, jittered by
// counter-based randomness of (seed, i).  Points are therefore a pure function of i, so CPU and GPU runs,
// and float and double runs, classify identical points.
#pragma once

#include "orbit.h"
#include <cstdint>
#include <span>
namespace mandelbrot {

// An axis-aligned cell [x, x + w] × [y, y + h]
struct Leaf { double x, y, w, h; };

struct SampleParams {
  int m;             // Samples per leaf
  int strata;        // Each group of strata^2 samples is jittered on a strata × strata grid
  uint64_t seed;
  int64_t max_iter;
  int K;             // Thresholds 2^-ks[k] for k < K ≤ 32
  int ks[32];
};

// splitmix64 finalizer
__host__ __device__ static inline uint64_t mix64(uint64_t z) {
  z += 0x9e3779b97f4a7c15;
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
  z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
  return z ^ (z >> 31);
}

// Uniform double in [0, 1) from (seed, i, j)
__host__ __device__ static inline double uniform(const uint64_t seed, const uint64_t i, const int j) {
  return double(mix64(seed ^ mix64(2 * i + j)) >> 11) * 0x1p-53;
}

// Location of sample i
__host__ __device__ static inline void sample_point(const Leaf& l, const SampleParams& p, const int64_t i,
                                                    double& x, double& y) {
  const int j = int(i % p.m) % (p.strata * p.strata), jx = j % p.strata, jy = j / p.strata;
  x = l.x + (jx + uniform(p.seed, i, 0)) / p.strata * l.w;
  y = l.y + (jy + uniform(p.seed, i, 1)) / p.strata * l.h;
}

// Bit k is below(e, ks[k])
__host__ __device__ static inline uint32_t below_bits(const Escape& e, const SampleParams& p) {
  uint32_t b = 0;
  for (int k = 0; k < p.K; k++)
    b |= uint32_t(e.steps < 0 || e.log2g < -p.ks[k]) << k;
  return b;
}

// CPU threads to use: $MANDELBROT_THREADS if set (say, a container's CPU request), else all hardware threads
int cpu_threads();

// Classify all leaves.size() * m samples into bits (one word per sample), returning total iterations.
// T is the orbit precision (float or double).
template<class T> int64_t sample_leaves_cpu(std::span<const Leaf> leaves, const SampleParams& p,
                                            std::span<uint32_t> bits);
template<class T> int64_t sample_leaves_cuda(std::span<const Leaf> leaves, const SampleParams& p,
                                             std::span<uint32_t> bits);  // Dies without CUDA

}  // namespace mandelbrot
