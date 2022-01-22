// Series area tests

#include "area.h"
#include "gtest/gtest.h"
namespace mandelbrot {
namespace {

TEST(area, double) {
  const int max_k = 7;
  areas<double>(max_k);
}

}  // namespace
}  // namespace mandelbrot
