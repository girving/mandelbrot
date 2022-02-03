// Cuda utilities
#pragma once

#include "arith.h"
#include "debug.h"
#include "device.h"
#include "noncopyable.h"
#include <cuda.h>
#include <driver_types.h>
#include <memory>
#include <type_traits>
namespace mandelbrot {

using std::add_const_t;
using std::conditional_t;
using std::is_const_v;
using std::is_signed_v;
using std::min;
using std::remove_const_t;
using std::shared_ptr;
using std::type_identity_t;
using std::unique_ptr;

void __attribute__((noreturn, cold))
cuda_check_fail(cudaError_t code, const char* function, const char* file, unsigned int line,
                const char* expression, const string& message);
#define cuda_check(code, ...) ({ \
  auto _code = (code); \
  if (_code != cudaSuccess) \
    cuda_check_fail(_code, __PRETTY_FUNCTION__, __FILE__, __LINE__, #code, format(__VA_ARGS__)); })

// For now, we share one stream for simplicity
CUstream stream();
void cuda_sync();

struct CudaDeleter {
  template<class T> void operator()(T* p) const { cuda_check(cudaFreeAsync(p, stream())); }
};
template<class T> using DevicePtr = unique_ptr<Device<T>,CudaDeleter>;

// span<Device<T>> but with ownership
template<class T> struct DeviceArray : public Noncopyable {
  typedef remove_const_t<T> S;
  typedef conditional_t<is_const_v<T>,const Device<S>,Device<S>> DT;
  typedef Device<T>* Data;

protected:
  DevicePtr<T> x;
  int64_t size_;
public:

  DeviceArray(const int64_t size) : size_(relu(size)) {
    if (size_) {
      Device<T>* p;
      cuda_check(cudaMallocAsync(&p, size_ * sizeof(T), stream()));
      x.reset(p);
    }
  }

  void swap(DeviceArray& y) { x.swap(y.x); std::swap(size_, y.size_); }
  void clear() { x.reset(); size_ = 0; }

  friend auto device_get(const DeviceArray& p) { return device_get(p.x.get()); }
  auto span() const { return std::span<DT>(x.get(), size_t(size_)); }
  operator std::span<DT>() const { return span(); }
  operator std::span<const DT>() const { return span(); }
  Device<T>* data() const { return x.get(); }
};

// Unpack a Device<T>* into a T* for use in kernel invocations
template<class T> static inline T* device_get(Device<T>* p) { return reinterpret_cast<T*>(p); }
template<class T> static inline const T* device_get(const Device<T>* p) { return reinterpret_cast<const T*>(p); }
template<class T> static inline auto device_get(const DevicePtr<T>& p) { return device_get(p.get()); }
template<class T> static inline auto device_get(span<Device<T>> p) { return device_get(p.data()); }
template<class T> static inline auto device_get(span<const Device<T>> p) { return device_get(p.data()); }

// Host to device and back.
// We use synchronous copies since our high performance code will be entirely GPU resident.
template<class T> static inline void host_to_device(span<Device<T>> dst, type_identity_t<span<const T>> src) {
  slow_assert(dst.size() == src.size());
  cuda_check(cudaMemcpy(device_get(dst.data()), src.data(), src.size()*sizeof(T), cudaMemcpyHostToDevice));
}
template<class T> static inline void device_to_host(span<T> dst, type_identity_t<span<const Device<T>>> src) {
  slow_assert(dst.size() == src.size());
  cuda_check(cudaMemcpy(dst.data(), device_get(src.data()), src.size()*sizeof(T), cudaMemcpyDeviceToHost));
}

// For device to device, we copy asynchronously
template<class T> static inline void device_to_device(span<Device<T>> dst, type_identity_t<span<const Device<T>>> src) {
  slow_assert(dst.size() == src.size());
  cuda_check(cudaMemcpyAsync(device_get(dst.data()), device_get(src.data()),
                             src.size()*sizeof(T), cudaMemcpyDeviceToDevice, stream()));
}

// Number of SMs
int num_sms();

// One-dimensional grid-stride loops
#define INVOKE_GRID_STRIDE_LOOP(name, n, ...) ({ \
  const int _n = (n);  /* For now, assume we fit in int32_t */ \
  const int threads = (_n + 31) >> 5 << 5;  /* Round up to a multiple of 32 */ \
  const int grids = min(32*num_sms(), (_n + threads - 1) / threads); \
  name<<<grids, threads>>>(_n, __VA_ARGS__); })

// For now, assume we fit in int32_t
template<class I> __device__ static inline int grid_stride_loop_size(const I n) {
  static_assert(is_signed_v<I>);
  return int(n);
}

#define GRID_STRIDE_LOOP(n, i) \
  for (int _n = grid_stride_loop_size(n), \
           _stride = blockDim.x * gridDim.x, \
           i = blockIdx.x * blockDim.x + threadIdx.x; \
       i < _n; i += _stride)

}  // namespace mandelbrot
