// Mandelbrot area via Arb

#include "arb_area.h"
#include "area.h"
#include "argparse.hpp"
#include "debug.h"
#include "device.h"
#include "expansion.h"
#include "print.h"
#include <exception>
#include <functional>
#include <map>
#include <vector>

using std::function;
using std::map;
using std::min;
using std::numeric_limits;
using std::vector;
using namespace mandelbrot;
static const double inf = numeric_limits<double>::infinity();

int main(int argc, char** argv) {
  try {
    argparse::ArgumentParser program("mandelbrot", "0.1");

    // Options and computation modes
    const auto mode = [&]() { return program.get("mode"); };
    const auto k = [&]() { return int(min(1000.0, program.get<double>("k"))); };
    const auto tol = [&]() { return program.get<double>("tol"); };
    const auto prec = [&]() { return program.get<int>("prec"); };
    const map<string,function<void()>> modes = {
      {"arb", [&]() { arb_areas(k(), prec()); }},
      {"double", [&]() { areas<double>(k(), tol()); }},
      {"exp2", [&]() { areas<Expansion<2>>(k(), tol()); }},
      IF_CUDA({"cuda-double", [&]() { areas<Device<double>>(k(), tol()); }},)
      IF_CUDA({"cuda-exp2", [&]() { areas<Device<Expansion<2>>>(k(), tol()); }},)
    };
    vector<string> mode_names;
    for (const auto& m : modes) mode_names.push_back(m.first);

    // Parse arguments
    program.add_argument("mode").help(format("computation method (one of %s)", mode_names))
        .default_value(string("exp2"));
    program.add_argument("-p", "--prec").help("precision, if we're using arb").scan<'i',int>().default_value(2000);
    program.add_argument("-t", "--tol").help("bail if error exceeds this").scan<'g',double>().default_value(inf);
    program.add_argument("-k").help("stop after 2^k terms").scan<'g',double>().default_value(inf);
    program.parse_args(argc, argv);
    if (program.is_used("--prec") && mode() != "arb")
      die("--prec only makes sense for mode == 'arb'");

    // Compute!
    print("mode = %s\n", mode());
    const auto it = modes.find(mode());
    slow_assert(it != modes.end());
    it->second();
    return 0;
  } catch (const std::exception& e) {
    die(e.what());
  }
}
