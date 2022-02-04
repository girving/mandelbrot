// Mandelbrot area via Arb

#include "arb_area.h"
#include "area.h"
#include "debug.h"
#include "device.h"
#include "expansion.h"
#include <exception>
#include <functional>
#include <map>
#include <vector>

using std::function;
using std::map;
using std::vector;
using namespace mandelbrot;

int main(int argc, char** argv) {
  try {
    const vector<string> args(argv+1, argv+argc);
    string cmd = "exp2";
    slow_assert(args.size() <= 1, "Expected zero or one arguments, not %s", args);
    if (args.size())
      cmd = args[0];

    const int max_k = 22;
    const double tol = 1;
    const int prec = 2000;
    const map<string,function<void()>> cmds = {
      {"arb", [&]() { arb_areas(max_k, prec); }},
      {"double", [&]() { areas<double>(max_k, tol); }},
      {"cuda-double", [&]() { areas<Device<double>>(max_k, tol); }},
      {"exp2", [&]() { areas<Expansion<2>>(max_k, tol); }},
      {"cuda-exp2", [&]() { areas<Device<Expansion<2>>>(max_k, tol); }},
    };

    const auto it = cmds.find(cmd);
    if (it != cmds.end()) {
      it->second();
      return 0;
    } else {
      vector<string> options;
      for (const auto& p : cmds)
        options.push_back(p.first);
      die("Unknown command '%s'.  Expected one of %s.", cmd, options);
    }
  } catch (const std::exception& e) {
    die(e.what());
  }
}
