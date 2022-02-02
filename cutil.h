// Cuda utilities
#pragma once

#include "debug.h"
#include "device.h"
#include <cuda.h>
#include <driver_types.h>
#include <memory>
#include <type_traits>
namespace mandelbrot {

using std::add_const_t;
using std::conditional_t;
using std::is_const_v;
using std::remove_const_t;
using std::shared_ptr;
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
  template<class T> void operator()(T* p) const { cuda_check(cudaFree(p)); }
};
template<class T> using DevicePtr = unique_ptr<Device<T>,CudaDeleter>;

// span<Device<T>> but with ownership
template<class T> struct DeviceArray {
  typedef remove_const_t<T> S;
  typedef conditional_t<is_const_v<T>,const Device<S>,Device<S>> DT;

  DevicePtr<T> x;
  int64_t n;

  void clear() { x.reset(); n = 0; }
  operator span<DT>() const { return span<DT>(x.get(), size_t(n)); }
  operator span<const DT>() const { return span<const DT>(x.get(), size_t(n)); }
};

template<class T> DeviceArray<T> alloc(const int64_t n) {
  slow_assert(n >= 0);
  Device<T>* t;
  cuda_check(cudaMalloc(&t, n * sizeof(T)));
  return DeviceArray<T>{DevicePtr<T>(t), n};
}

// Unpack a Device<T>* into a T* for use in kernel invocations
template<class T> static inline T* device_get(Device<T>* p) { return reinterpret_cast<T*>(p); }
template<class T> static inline const T* device_get(const Device<T>* p) { return reinterpret_cast<const T*>(p); }
template<class T> static inline auto device_get(const DevicePtr<T>& p) { return device_get(p.get()); }
template<class T> static inline auto device_get(const DeviceArray<T>& p) { return device_get(p.x.get()); }
template<class T> static inline auto device_get(span<Device<T>> p) { return device_get(p.data()); }
template<class T> static inline auto device_get(span<const Device<T>> p) { return device_get(p.data()); }

// Host to device and back.
// We use synchronous copies since our high performance code will be entirely GPU resident.
template<class T> static inline void host_to_device(span<Device<T>> dst, span<const T> src) {
  slow_assert(dst.size() == src.size());
  cuda_check(cudaMemcpy(device_get(dst.data()), src.data(), src.size()*sizeof(T), cudaMemcpyHostToDevice));
}
template<class T> static inline void device_to_host(span<T> dst, span<const Device<T>> src) {
  slow_assert(dst.size() == src.size());
  cuda_check(cudaMemcpy(dst.data(), device_get(src.data()), src.size()*sizeof(T), cudaMemcpyDeviceToHost));
}

}  // namespace mandelbrot
