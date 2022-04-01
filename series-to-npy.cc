// Convert series text files to .npy

#include "expansion.h"
#include "series.h"
#include <unistd.h>

using namespace mandelbrot;

string numpy_header(const vector<int>& shape) {
  // Assume little endian
  const char endian = '<';
  // Assume double dtype
  const char letter = 'f';
  const int bytes = 8;
  string header("\x93NUMPY\x01\x00??", 10);
  header += format("{'descr': '%c%c%d', 'fortran_order': False, 'shape': (", endian, letter, bytes);
  for (const int n : shape)
    header += format("%d,", n);
  header += "), }";
  while ((header.size()+1) & 15)
    header.push_back(' ');
  header.push_back('\n');
  const auto header_size = uint16_t(header.size()-10);
  memcpy(header.data() + 8, &header_size, 2);
  return header;
} 

void convert(const string& input, const string& output) {
  // Verify that output file doesn't exist
  print("converting %s to %s", input, output);
  slow_assert(output.ends_with(".npy"));
  slow_assert(access(output.c_str(), F_OK) == -1, "%s exists", output);

  // Read series, assuming exp2 for now
  typedef Expansion<2> S;
  const auto [_, x] = read_series<S>(input); 

  // Write .npy file
  FILE* f = fopen(output.c_str(), "wb");
  slow_assert(f);
  const auto header = numpy_header({x.nonzero(), 2});
  fwrite(header.c_str(), header.size(), 1, f);
  fwrite(x.data(), sizeof(S)*x.nonzero(), 1, f);
  fclose(f);
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
