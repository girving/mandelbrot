// Assertions

#include "debug.h"
#include <signal.h>
#include <iostream>
namespace mandelbrot {

using std::runtime_error;

void assertion_failed(const char* function, const char* file, unsigned int line,
                      const char* condition, const string& message) {
  const string error = format("%s:%d:%s: %s, condition = %s", file, line, function,
                              message.size() ? message : "Assertion failed", condition);
  static const bool break_on_assert = getenv("BREAK_ON_ASSERT") != 0;
  if (break_on_assert) {
    std::cout << std::flush;
    std::cerr << "\n\n*** Error: " << error << '\n' << std::endl;
    raise(SIGINT);
  }
  throw runtime_error(error);
}

void die(const string& msg) {
  if (msg.size())
    std::cerr << "\nerror: " << msg << std::endl;
  if (getenv("BREAK_ON_ASSERT"))
    raise(SIGINT);
  exit(1);
}

}  // namespace mandelbrot
