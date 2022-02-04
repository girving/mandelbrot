// Series area tests

#include "area.h"
#include "expansion.h"
#include "tests.h"
namespace mandelbrot {
namespace {

TEST(double) {
  const int max_k = 7;
  const double tol = 1e-14;
  areas<double>(max_k, tol);
}

TEST(expansion2) {
  const int max_k = 7;
  const double tol = 2e-31;
  areas<Expansion<2>>(max_k, tol);
}

}  // namespace
}  // namespace mandelbrot
