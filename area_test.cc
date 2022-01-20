// Arb area tests

#include "area.h"
#include "gtest/gtest.h"
namespace mandelbrot {
namespace {

TEST(area, area) {
  const int max_k = 7;
  const int prec = 200;
  areas(max_k, prec);
}

}  // namespace
}  // namespace mandelbrot
