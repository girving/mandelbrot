// Very simple testing infrastructure

#include "tests.h"
#include "print.h"
#include "shutdown.h"
#include "wall_time.h"
#include <functional>
namespace mandelbrot {

using std::exception;
using std::function;
using std::tuple;

typedef vector<tuple<string,function<void()>>> Tests;
static Tests& tests() { static Tests tests; return tests; }

int register_test(const char* name, void (*test)()) {
  tests().emplace_back(name, test);
  return 0;
}

void test_throw_fail(const char* sx, const char* se, const char* function, const int line) {
  print("%s %s didn't throw %s", red(format("%s:%d:", function, line)), sx, se);
  throw test_error();
}

// See https://stackoverflow.com/questions/2616906/how-do-i-output-coloured-text-to-a-linux-terminal
static string color(const string& s, const int color) { return format("\033[%dm%s\033[0m", color, s); }
string green(const string& s) { return color(s, 32); }
string red(const string& s) { return color(s, 31); }

static int run_tests(const vector<string>& args) {
  const auto count = [](const int n) { return format("%d %s", n, n == 1 ? "test" : "tests"); };
  wall_time_t total;
  print("%s %s: %s", green("[==========]"), args[0], count(tests().size()));
  int passed = 0;
  vector<string> failed;
  for (const auto& [name, test] : tests()) {
    print("%s %s", green("[ RUN      ]"), name);
    const auto start = wall_time();
    bool good;
    try {
      test();
      good = true;
    } catch (const test_error& e) {
      good = false;
    } catch (const exception& e) {
      print("  %s %s, %s", red("exception:"), typeid(e).name(), e.what());
      good = false;
    }
    const auto elapsed = wall_time() - start;
    total += elapsed;
    print("%s %s (%d ms)", good ? green("[       OK ]") : red("[  FAILED  ]"),
          name, int(rint(elapsed.milliseconds())));
    if (good) passed++;
    else failed.push_back(name);
  }

  // Finish up
  print("%s %s ran (%d ms total)", green("[==========]"), count(tests().size()), int(rint(total.milliseconds())));
  print("%s %s", green("[  PASSED  ]"), count(passed));
  if (failed.size()) {
    print("%s %s", red("[  FAILED  ]"), count(failed.size()));
    for (const auto& name : failed)
      print("%s %s", red("[  FAILED  ]"), name);
  }
  shutdown();
  return failed.size() ? 1 : 0;
}

}  // namespace mandelbrot
using namespace mandelbrot;

int main(int argc, char** argv) {
  try {
    return run_tests({argv, argv + argc});
  } catch (const exception& e) {
    die(e.what());
  }
}
