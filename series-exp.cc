// Exponentiate a series

#include "expansion_arith.h"
#include "series.h"
#include <unistd.h>

using namespace mandelbrot;

void exp(const string& input, const string& output) {
  // Verify that output file doesn't exist
  print("%s = exp(%s)", output, input);
  slow_assert(access(output.c_str(), F_OK) == -1, "%s exists", output);

  // Read series, assuming exp2 for now
#ifndef __CUDACC__
  typedef Expansion<2> S;
#else
  typedef Device<Expansion<2>> S;
#endif
  print("reading...");
  const auto [_, x] = read_series<S>(input);

  // Exponentiate
  print("exponentiating...");
  Series<S> y(x.nonzero());
  y = exp(x);

  // Write output
  print("writing...");
  write_series_npy<S>(output, y);
  print("done");
}

int main(const int argc, const char** argv) {
  try {
    slow_assert(argc == 3, "usage %s <input> <output.npy>", argv[0]);
    exp(argv[1], argv[2]);
    return 0;
  } catch (const std::exception& e) {
    die(e.what());
  }
}
