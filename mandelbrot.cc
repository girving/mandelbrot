// Mandelbrot area via Arb

#include "arb_area.h"
#include "area.h"
#include "argparse.hpp"
#include "debug.h"
#include "device.h"
#include "expansion.h"
#include "print.h"
#include <algorithm>
#include <exception>
#include <fstream>
#include <functional>
#include <map>
#include <sys/stat.h>
#include <vector>

using std::endl;
using std::function;
using std::map;
using std::min;
using std::ofstream;
using std::numeric_limits;
using std::vector;
using namespace mandelbrot;
static const double inf = numeric_limits<double>::infinity();

// Write a series to a file in a simple text format
template<class T,bool v> void write_series(const string& path, vector<string> comments, const Series<T,v>& x) {
  ofstream out(path);

  // From https://stackoverflow.com/questions/16357999/current-date-and-time-as-string/16358264
  const auto t = std::time(nullptr);
  const auto tm = *std::localtime(&t);
  comments.push_back(format("time = %s", std::put_time(&tm, "%F %T")));

  for (const auto& c : comments)
    out << "# " << c << endl;

  // Write series terms in plain text
  for (const auto& t : host_copy(x))
    out << safe(t) << endl;
}

template<class T> void areas(const string& output, const string& mode, const int max_k, const double tol) {
  auto g = bottcher_base<T>();
  for (int k = 1; k <= max_k; k++) {
    const auto [f, mu] = bottcher_step(g, tol);
    if (output.size()) {
      const auto write = [&output,&mode,k,mu=mu](const string& n, const string& name, const auto& x) {
        write_series(
            format("%s/%c-k%d", output, n, k),
            {name,
             format("mode = %s", mode),
             format("k = %d", k),
             format("terms = %d", x.known()),
             format("mu = %s", safe(mu))},
            x);
      };
      write("g", "g = log(f)", g);
      write("f", "f = f(z) = 1/phi(1/z)", f);
    }
  }
}

int main(int argc, char** argv) {
  try {
    argparse::ArgumentParser program("mandelbrot", "0.1");

    // Options and computation modes
    const auto mode = [&]() { return program.get("mode"); };
    const auto k = [&]() { return int(std::min(1000.0, program.get<double>("k"))); };
    const auto tol = [&]() { return program.get<double>("tol"); };
    const auto prec = [&]() { return program.get<int>("prec"); };
    const auto output = [&]() { return program.is_used("--output") ? program.get("output") : ""; };
    #define AREAS(T) [&]() { areas<T>(output(), mode(), k(), tol()); }
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
    program.add_argument("mode").help(format("computation method (one of %s)", mode_names))
        .default_value(string("exp2"));
    program.add_argument("-k").help("stop after 2^k terms").scan<'g',double>().default_value(inf);
    program.add_argument("-p", "--prec").help("precision, if we're using arb").scan<'i',int>().default_value(2000);
    program.add_argument("-t", "--tol").help("bail if error exceeds this").scan<'g',double>().default_value(inf);
    program.add_argument("-o", "--output").help("directory to write results to");
    program.parse_args(argc, argv);
    if (program.is_used("--prec") && mode() != "arb")
      die("--prec only makes sense for mode == 'arb'");

    // Make output directory
    if (output().size()) {
      const int r = mkdir(output().c_str(), 0777);
      if (r) die("failed to create directory '%s': %s", output(), strerror(errno));
      tee(output() + "/log");
    }

    // Log command line options
    print("mode = %s", mode());
    print("k = %g", program.get<double>("k"));
    print("tol = %g", tol());
    if (mode() == "arb") print("prec = %d", prec());
    if (output().size()) print("output = '%s'", output());
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
