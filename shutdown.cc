// Graceful shutdown on exit

#include "shutdown.h"
#include <vector>
namespace mandelbrot {

using std::vector;

static vector<function<void()>>& list() {
  static vector<function<void()>> list;
  return list;
}

void on_shutdown(const function<void()>& f) {
  list().push_back(f);
}

void shutdown() {
  for (const auto& f : list())
    f();
}

}  // namespace mandelbrot
