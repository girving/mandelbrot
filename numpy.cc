// .npy format

#include "numpy.h"
#include "format.h"
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

}  // namespace mandelbrot
