// .npy format

#include "numpy.h"
#include "debug.h"
#include "format.h"
#include <cstring>
#include <fstream>
namespace mandelbrot {

string numpy_header(const vector<int>& shape) {
  // Assume little endian
  const char endian = '<';
  // Assume double dtype
  const char letter = 'f';
  const int bytes = 8;
  string header("\x93NUMPY\x01\x00??", 10);
  header += tfm::format("{'descr': '%c%c%d', 'fortran_order': False, 'shape': (", endian, letter, bytes);
  for (const int n : shape)
    header += tfm::format("%d,", n);
  header += "), }";
  while ((header.size()+1) & 15)
    header.push_back(' ');
  header.push_back('\n');
  const auto header_size = uint16_t(header.size()-10);
  memcpy(header.data() + 8, &header_size, 2);
  return header;
} 

Numpy read_numpy(const string& path) {
  std::ifstream in(path, std::ios::binary);
  slow_assert(in, "failed to open %s", path);
  char magic[8];
  in.read(magic, 8);
  slow_assert(in && !memcmp(magic, "\x93NUMPY", 6), "%s: not a .npy file", path);
  const int version = magic[6];
  slow_assert(version == 1 || version == 2, "%s: unsupported .npy version %d", path, version);
  uint32_t header_size = 0;
  in.read(reinterpret_cast<char*>(&header_size), version == 1 ? 2 : 4);
  string header(header_size, ' ');
  in.read(header.data(), header_size);
  slow_assert(in, "%s: truncated header", path);
  slow_assert(header.find("'descr': '<f8'") != string::npos, "%s: expected dtype <f8, got %s", path, header);
  slow_assert(header.find("'fortran_order': False") != string::npos, "%s: expected C order", path);

  // Parse shape tuple, e.g. (3, 4) or (5,) or ()
  const auto s = header.find("'shape': (");
  slow_assert(s != string::npos, "%s: no shape in header", path);
  Numpy r;
  int64_t size = 1;
  for (size_t i = s + 10; header[i] != ')';) {
    if (header[i] == ',' || header[i] == ' ') { i++; continue; }
    slow_assert(isdigit(header[i]), "%s: bad shape in header %s", path, header);
    int64_t n = 0;
    for (; isdigit(header[i]); i++)
      n = 10 * n + (header[i] - '0');
    r.shape.push_back(n);
    size *= n;
  }

  r.data.resize(size);
  in.read(reinterpret_cast<char*>(r.data.data()), size * sizeof(double));
  slow_assert(in, "%s: expected %d doubles", path, size);
  return r;
}

}  // namespace mandelbrot
