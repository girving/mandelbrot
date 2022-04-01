// Convert series text files to .npy

#include "expansion.h"
#include "series.h"
#include <unistd.h>

using namespace mandelbrot;

void convert(const string& input, const string& output) {
  // Verify that output file doesn't exist
  print("converting %s to %s", input, output);
  slow_assert(output.ends_with(".npy"));
  slow_assert(access(output.c_str(), F_OK) == -1, "%s exists", output);

  // Read series, assuming exp2 for now
  typedef Expansion<2> S;
  const auto [_, x] = read_series<S>(input); 

  // Write .npy file
  write_series_npy<S>(output, x);
}

int main(const int argc, const char** argv) {
  try {
    slow_assert(argc == 3, "usage %s <input> <output.npy>", argv[0]);
    convert(argv[1], argv[2]);
    return 0;
  } catch (const std::exception& e) {
    die(e.what());
  }
}
