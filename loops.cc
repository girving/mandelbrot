// Cuda/host loops

#include "loops.h"
namespace mandelbrot {

using std::make_tuple;
using std::min;

// Chop a loop into [start,end) chunks
tuple<int64_t,int64_t> partition_loop(const int64_t steps, const int threads, const int thread) {
  slow_assert(threads > 0 && unsigned(thread) < unsigned(threads) && steps >= 0);
  const int64_t steps_per_thread = steps / threads;  // Round down, so some threads will get one more step
  const int extra_steps = steps % threads;  // The first extra_steps threads will get one extra step
  const int64_t start = steps_per_thread * thread + min(extra_steps, thread);
  const int64_t end = start + steps_per_thread + (unsigned(thread) < unsigned(extra_steps));
  return make_tuple(start, end);
}

}  // namespace mandelbrot
