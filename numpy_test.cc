// .npy tests

#include "numpy.h"
#include "tests.h"
#include <fstream>
namespace mandelbrot {
namespace {

using std::runtime_error;

void write(const string& path, const vector<int>& shape, const vector<double>& data) {
  std::ofstream out(path, std::ios::binary);
  const auto h = numpy_header(shape);
  out.write(h.data(), h.size());
  out.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(double));
}

TEST(roundtrip) {
  const Tmpfile tmp("numpy_test");
  for (const auto& shape : {vector<int>{2, 3}, vector<int>{5}, vector<int>{1, 1, 4}, vector<int>{0, 7}}) {
    int64_t size = 1;
    for (const int n : shape) size *= n;
    vector<double> data(size);
    for (int64_t i = 0; i < size; i++)
      data[i] = 0.5 * i - 1.25;
    write(tmp.path, shape, data);
    const auto r = read_numpy(tmp.path);
    ASSERT_EQ(r.shape.size(), shape.size());
    for (size_t i = 0; i < shape.size(); i++)
      ASSERT_EQ(r.shape[i], shape[i]);
    ASSERT_EQ(r.data, data);
  }
}

TEST(truncated) {
  const Tmpfile tmp("numpy_test");
  write(tmp.path, {4}, {1, 2, 3});  // One double short
  ASSERT_THROW(read_numpy(tmp.path), runtime_error);
}

}  // namespace
}  // namespace mandelbrot
