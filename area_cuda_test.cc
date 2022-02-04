// Series area CUDA tests

#include "area.h"
#include "device.h"
#include "expansion.h"
#include "tests.h"
namespace mandelbrot {
namespace {

TEST(cuda_double) {
  const int max_k = 7;
  const double tol = 1e-14;
  areas<Device<double>>(max_k, tol);
}

/*
// TODO: Re-enable
TEST(cuda_expansion2) {
  const int max_k = 7;
  const double tol = 6e-32;
  areas<Device<Expansion<2>>>(max_k, tol);
}
*/

}  // namespace
}  // namespace mandelbrot
