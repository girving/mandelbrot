// Arb area tests

#include "arb_area.h"
#include "gtest/gtest.h"
namespace mandelbrot {
namespace {

TEST(arb_area, areas) {
  const int max_k = 7;
  const int prec = 200;
  arb_areas(max_k, prec);
}

}  // namespace
}  // namespace mandelbrot
