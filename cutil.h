// Cuda utilities
#pragma once

#ifndef __CUDACC__

#define __host__
#define __device__
#define cuda_check(...) ((void)0)
#define IF_CUDA(...)
#define CUDA_OR_DIE(...) die("No CUDA")
#include "debug.h"

#else  // __CUDACC__

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

#define IF_CUDA(...) __VA_ARGS__
#define CUDA_OR_DIE(...) __VA_ARGS__

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

// Unpack a Device<T>* into a T* for use in kernel invocations
template<class T> static inline T* device_get(Device<T>* p) { return reinterpret_cast<T*>(p); }
template<class T> static inline const T* device_get(const Device<T>* p) { return reinterpret_cast<const T*>(p); }
template<class C> static inline auto device_get(C&& c) { return device_get(c.data()); }

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

// One slow write
template<class T> static inline void single_host_to_device(Device<T>* dst, const T src) {
  cuda_check(cudaMemcpy(device_get(dst), &src, sizeof(T), cudaMemcpyHostToDevice));
}

// For device to device, we copy asynchronously
template<class T> static inline void device_to_device(span<Device<T>> dst, type_identity_t<span<const Device<T>>> src) {
  slow_assert(dst.size() == src.size());
  cuda_check(cudaMemcpyAsync(device_get(dst.data()), device_get(src.data()),
                             src.size()*sizeof(T), cudaMemcpyDeviceToDevice, stream()));
}

// Number of SMs
int num_sms();

}  // namespace mandelbrot
#endif  // __CUDACC__

#include "device.h"
#include "debug.h"
#include "span.h"
namespace mandelbrot {

using std::type_identity_t;

// Host to host or device to device
template<class T> static inline void same_to_same(span<T> dst, type_identity_t<span<const T>> src) {
  if constexpr (is_device<T>)
    CUDA_OR_DIE(device_to_device(dst, src));
  else {
    slow_assert(dst.size() == src.size());
    std::copy(src.data(), src.data() + src.size(), dst.data());
  }
}

}  // namespace mandelbrot
