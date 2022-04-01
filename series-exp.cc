// Exponentiate a series

#include "expansion_arith.h"
#include "series.h"
#include <unistd.h>

using namespace mandelbrot;

void exp(const string& input, const string& output) {
  // Verify that output file doesn't exist
  print("%s = exp(%s)", input, output);
  slow_assert(access(output.c_str(), F_OK) == -1, "%s exists", output);

  // Read series, assuming exp2 for now
  typedef Expansion<2> S;
  const auto [comments, x] = read_series<S>(input); 

  // Exponentiate
  Series<S> y(x.nonzero());
  y = exp(x);

  // Write output
  write_series<S>(output, comments, y);
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
