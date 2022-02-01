// Cuda utilities

#include "debug.h"
#include <cuda.h>
#include <memory>
namespace mandelbrot {

void __attribute__((noreturn, cold))
cuda_check_fail(cudaError_t code, const char* function, const char* file, unsigned int line,
                const char* expression, const string& message);
#define cuda_check(code, ...) ({ \
  auto _code = (code); \
  if (_code != cudaSuccess) \
    cuda_check_fail(_code, __PRETTY_FUNCTION__, __FILE__, __LINE__, #code, format(__VA_ARGS__)); })

using std::shared_ptr;
using std::unique_ptr;

struct Stream {
private:
  shared_ptr<CUstream_st> s;
public:

  Stream() {
    cudaStream_t p;
    cuda_check(cudaStreamCreate(&p));
    s.reset(p, [](cudaStream_t p) { cuda_check(cudaStreamDestroy(p)); });
  }

  cudaStream_t stream() const { return s.get(); }
  operator cudaStream_t() const { return stream(); }
  void sync() const { cuda_check(cudaStreamSynchronize(stream())); }
};

struct AsyncDeleter {
  cudaStream_t s;
  template<class T> void operator()(T* p) const { cuda_check(cudaFreeAsync(p, s)); }
};
template<class T> using DevicePtr = unique_ptr<T[],AsyncDeleter>;

template<class T> DevicePtr<T> alloc(const int64_t n, cudaStream_t s) {
  slow_assert(n >= 0);
  T* t;
  cuda_check(cudaMallocAsync(&t, n * sizeof(T), s));
  return DevicePtr<T>(t, AsyncDeleter{s});
}

}  // namespace mandelbrot
