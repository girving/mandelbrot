// CUDA utilities

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
    CUstream p;
    cuda_check(cudaStreamCreate(&p));
    s.reset(p, [](CUstream p) { cuda_check(cudaStreamDestroy(p)); });
  }

  CUstream stream() const { return s.get(); }
  operator CUstream() const { return stream(); }
};

CUstream stream() {
  static Stream s;
  return s;
}

void cuda_sync() {
  cuda_check(cudaStreamSynchronize(stream()));
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
