// Cuda utilities
#ifdef __CUDACC__

#include "cutil.h"
namespace mandelbrot {

void cuda_check_fail(cudaError_t code, const char* function, const char* file, unsigned int line,
                     const char* expression, const string& message) {
  die("%s:%d:%s: %s = %d (%s)%s",
      file, line, function, expression, int(code), cudaGetErrorString(code),
      message.size() ? ", " + message : "");
}

// Move to header if we want to go to multiple streams later
struct Stream {
private:
  shared_ptr<CUstream_st> s;
public:

  Stream() {
    // Keep freed stream-ordered allocations (Mem) in the device's pool: by default the pool returns them to the OS
    // at each synchronize, so every batch's buffers (gigabytes of leaf bits) would be mapped afresh
    static const bool pooled = [] {
      int device;
      cuda_check(cudaGetDevice(&device));
      cudaMemPool_t pool;
      cuda_check(cudaDeviceGetDefaultMemPool(&pool, device));
      uint64_t threshold = UINT64_MAX;
      cuda_check(cudaMemPoolSetAttribute(pool, cudaMemPoolAttrReleaseThreshold, &threshold));
      return true;
    }();
    (void)pooled;
    CUstream p;
    cuda_check(cudaStreamCreateWithFlags(&p, cudaStreamNonBlocking));
    s.reset(p, [](CUstream p) { cuda_check(cudaStreamDestroy(p)); });
  }

  CUstream stream() const { return s.get(); }
  operator CUstream() const { return stream(); }
};

CUstream stream() {
  const bool synchronous = false;
  if (synchronous)
    return 0;
  else {
    // One non-blocking stream per host thread, so that concurrent batches (TreeParams::overlap) overlap
    static thread_local Stream s;
    return s;
  }
}

void cuda_sync() {
  cuda_check(cudaStreamSynchronize(stream()));
}

static cudaMemPool_t default_pool() {
  int device;
  cuda_check(cudaGetDevice(&device));
  cudaMemPool_t pool;
  cuda_check(cudaDeviceGetDefaultMemPool(&pool, device));
  return pool;
}

string gpu_memory() {
  uint64_t used = 0, high = 0, reserved = 0;
  const auto pool = default_pool();
  cuda_check(cudaMemPoolGetAttribute(pool, cudaMemPoolAttrUsedMemCurrent, &used));
  cuda_check(cudaMemPoolGetAttribute(pool, cudaMemPoolAttrUsedMemHigh, &high));
  cuda_check(cudaMemPoolGetAttribute(pool, cudaMemPoolAttrReservedMemCurrent, &reserved));
  size_t free = 0, total = 0;
  cuda_check(cudaMemGetInfo(&free, &total));
  return tfm::format("GPU memory in use %.1f GB (peak %.1f), pool %.1f GB, device %.1f of %.1f GB", used * 1e-9,
                     high * 1e-9, reserved * 1e-9, (total - free) * 1e-9, total * 1e-9);
}

void* cuda_malloc(const size_t bytes) {
  void* p = nullptr;
  cudaError_t e = cudaMallocAsync(&p, bytes, stream());
  if (e == cudaErrorMemoryAllocation) {
    // The pool keeps freed memory (see Stream), possibly in pieces cached for other streams: return it all to the
    // device and try once more
    (void)cudaGetLastError();
    cuda_check(cudaDeviceSynchronize());
    cuda_check(cudaMemPoolTrimTo(default_pool(), 0));
    e = cudaMallocAsync(&p, bytes, stream());
  }
  if (e != cudaSuccess) die("cuda_malloc: %.3g GB: %s (%s)", bytes * 1e-9, cudaGetErrorString(e), gpu_memory());
  return p;
}

int num_sms() {
  static const int num_sms = []() {
    const int device = 0;  // This will need to change if we go multi-device
    int num_sms;
    cuda_check(cudaDeviceGetAttribute(&num_sms, cudaDevAttrMultiProcessorCount, device));
    return num_sms;
  }();
  return num_sms;
}

}  // namespace mandelbrot
#endif  // __CUDACC_
