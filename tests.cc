// Very simple testing infrastructure

#include "tests.h"
#include "print.h"
#include "shutdown.h"
#include "wall_time.h"
#include <functional>
#include <unistd.h>
#include <unordered_set>
namespace mandelbrot {

using std::exception;
using std::function;
using std::tuple;
using std::unordered_set;

typedef vector<tuple<string,function<void()>>> Tests;
static Tests& tests() { static Tests tests; return tests; }

int register_test(const char* name, void (*test)()) {
  tests().emplace_back(name, test);
  return 0;
}

void test_throw(const function<void()>& f, const char* sx, const char* se, const char* function, const int line) {
  try {
    f();
  } catch (const exception& e) {
    const string s = typeid(e).name();
    if (s.find(string(se)) != string::npos)
      return;
    print("%s %s threw %s, not %s", red(format("%s:%d:", function, line)), sx, s, se);
    throw test_error();
  }
  print("%s %s didn't throw %s", red(format("%s:%d:", function, line)), sx, se);
  throw test_error();
}

Tmpfile::Tmpfile(string_view prefix) {
  string p = "/tmp/" + string(prefix) + "XXXXXX";
  const int fd = mkstemp(p.data());
  slow_assert(fd >= 0, strerror(errno));
  const_cast<string&>(path) = p;
}

Tmpfile::~Tmpfile() { if (path.size()) unlink(path.c_str()); }

// See https://stackoverflow.com/questions/2616906/how-do-i-output-coloured-text-to-a-linux-terminal
static string color(const string& s, const int color) { return format("\033[%dm%s\033[0m", color, s); }
string red(const string& s) { return color(s, 31); }
string green(const string& s) { return color(s, 32); }
string blue(const string& s) { return color(s, 34); }

static int run_tests(const vector<string>& args) {
  const auto& tests = mandelbrot::tests();
  const unordered_set<string> chosen(args.begin() + 1, args.end());
  const auto skip = [&chosen](const string& name) {
    return chosen.size() && chosen.find(name) == chosen.end();
  };
  const auto count = [](const int n) { return format("%d %s", n, n == 1 ? "test" : "tests"); };
  wall_time_t total;
  print("%s %s: %s", green("[==========]"), args[0], count(tests.size()));

  // Log skipped tests
  for (const auto& [name, test] : tests)
    if (skip(name))
      print("%s %s", blue("[   SKIP   ]"), name);

  // Run non-skipped tests
  int passed = 0;
  vector<string> failed;
  for (const auto& [name, test] : tests) {
    if (skip(name)) continue;
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
  print("%s %s ran (%d ms total)", green("[==========]"), count(tests.size()), int(rint(total.milliseconds())));
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
