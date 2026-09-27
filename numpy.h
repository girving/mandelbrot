// .npy format

#include <cstdint>
#include <string>
#include <vector>
namespace mandelbrot {

using std::string;
using std::vector;

string numpy_header(const vector<int>& shape);

// Read a little endian float64 C-order .npy file (format version 1 or 2)
struct Numpy {
  vector<int64_t> shape;
  vector<double> data;
};
Numpy read_numpy(const string& path);

}  // namespace mandelbrot
