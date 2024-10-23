// Mandelbrot area via Arb

#include "arb_area.h"
#include "area.h"
#include "argparse.hpp"
#include "debug.h"
#include "device.h"
#include "expansion.h"
#include "join.h"
#include "print.h"
#include <algorithm>
#include <exception>
#include <functional>
#include <map>
#include <sys/stat.h>
#include <vector>

using std::endl;
using std::function;
using std::map;
using std::min;
using std::numeric_limits;
using std::vector;
using namespace mandelbrot;
static const double inf = numeric_limits<double>::infinity();

template<class T> void areas(const optional<string>& input, const optional<string>& output,
                             const string& mode, const int max_k, const double tol) {
  // Initialize, either from scratch or from file
  Series<T> g(input ? read_bottcher<T>(*input) : bottcher_base<T>());

  // Learn more terms
  const int k0 = known_to_k(g.known());
  for (int k = k0+1; k <= max_k; k++) {
    const auto [f, mu] = bottcher_step(g, tol);
    if (output)
      write_bottcher<T>(*output, mode, mu, f, g);
  }
}

int main(const int argc, const char** argv) {
  try {
    argparse::ArgumentParser program("mandelbrot", "0.1");

    // Options and computation modes
    const auto mode = [&]() { return program.get("mode"); };
    const auto k = [&]() { return int(std::min(1000.0, program.get<double>("k"))); };
    const auto tol = [&]() { return program.get<double>("tol"); };
    const auto prec = [&]() { return program.get<int>("prec"); };
    const auto input = [&]() { return program.present("input"); };
    const auto output = [&]() { return program.present("output"); };
    #define AREAS(T) [&]() { areas<T>(input(), output(), mode(), k(), tol()); }
    const map<string,function<void()>> modes = {
      {"arb", [&]() { arb_areas(k(), prec()); }},
      {"double", AREAS(double)},
      {"exp2", AREAS(Expansion<2>)},
      IF_CUDA({"cuda-double", AREAS(Device<double>)},)
      IF_CUDA({"cuda-exp2", AREAS(Device<Expansion<2>>)},)
    };
    vector<string> mode_names;
    for (const auto& m : modes) mode_names.push_back(m.first);

    // Parse arguments
    program.add_argument("mode").help(tfm::format("computation method (one of %s)", mode_names))
        .default_value(string("exp2"));
    program.add_argument("-k").help("stop after 2^k terms").scan<'g',double>().default_value(inf);
    program.add_argument("-p", "--prec").help("precision, if we're using arb").scan<'i',int>().default_value(2000);
    program.add_argument("-t", "--tol").help("bail if error exceeds this").scan<'g',double>().default_value(inf);
    program.add_argument("-i", "--input").help("g = log f series to resume from");
    program.add_argument("-o", "--output").help("directory to write results to");
    program.parse_args(argc, argv);
    if (program.is_used("--prec") && mode() != "arb")
      die("--prec only makes sense for mode == 'arb'");

    // Make output directory
    if (output()) {
      const int r = mkdir(output()->c_str(), 0777);
      if (r) die("failed to create directory '%s': %s", *output(), strerror(errno));
      tee(*output() + "/log");
    }

    // Log command line options
    print("cmd = %s\n", join(span<const char*>(argv, argv+argc), " "));
    print("mode = %s", mode());
    print("k = %g", program.get<double>("k"));
    print("tol = %g", tol());
    if (mode() == "arb") print("prec = %d", prec());
    if (input()) print("input = '%s'", *input());
    if (output()) print("output = '%s'", *output());
    print();

    // Compute!
    const auto it = modes.find(mode());
    slow_assert(it != modes.end());
    it->second();
    return 0;
  } catch (const std::exception& e) {
    die(e.what());
  }
}
