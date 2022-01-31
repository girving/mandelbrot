// Series area tests

#include "area.h"
#include "expansion.h"
#include "gtest/gtest.h"
namespace mandelbrot {
namespace {

TEST(area, double) {
  const int max_k = 7;
  const double tol = 1e-6;
  areas<double>(max_k, tol);
}

TEST(area, expansion2) {
  const int max_k = 7;
  const double tol = 6e-32;
  areas<Expansion<2>>(max_k, tol);
}

}  // namespace
}  // namespace mandelbrot
