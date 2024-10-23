// Assertions

#include "debug.h"
#include "print.h"
#include <signal.h>
#include <iostream>
namespace mandelbrot {

using std::runtime_error;

void assertion_failed(const char* function, const char* file, unsigned int line,
                      const char* condition, const string& message) {
  const string error = tfm::format("%s:%d:%s: %s, condition = %s", file, line, function,
                                   message.size() ? message : "Assertion failed", condition);
  static const bool break_on_assert = getenv("BREAK_ON_ASSERT") != 0;
  if (break_on_assert) {
    print_error("\n\n*** Error: %s\n", error);
    raise(SIGINT);
  }
  throw runtime_error(error);
}

void die(const string& msg) {
  if (msg.size())
    print_error("\nerror: %s", msg);
  if (getenv("BREAK_ON_ASSERT"))
    raise(SIGINT);
  exit(1);
}

}  // namespace mandelbrot
