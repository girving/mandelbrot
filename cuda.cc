// CUDA utilities

#include "cuda.h"
namespace mandelbrot {

void cuda_check_fail(cudaError_t code, const char* function, const char* file, unsigned int line,
                     const char* expression, const string& message) {
  die("%s:%d:%s: %s = %d (%s)%s",
      file, line, function, expression, code, cudaGetErrorString(code),
      message.size() ? ", " + message : "");
}

}  // namespace mandelbrot
