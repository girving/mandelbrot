// Arb area tests

#include "arb_area.h"
#include "tests.h"
namespace mandelbrot {
namespace {

TEST(areas) {
  const int max_k = 7;
  const int prec = 200;
  arb_areas(max_k, prec);
}

}  // namespace
}  // namespace mandelbrot
