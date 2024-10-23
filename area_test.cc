// Series area tests

#include "area.h"
#include "expansion.h"
#include "tests.h"
namespace mandelbrot {
namespace {

template<class S> void area_test(const double tol) {
  const int max_k = 7;
  Tmpdir tmp("area");

  // Compute from scratch
  S last_mu;
  Series<S> last_g;
  {
    auto g = bottcher_base<S>();
    for (int k = 1; k <= max_k; k++) {
      const auto [f, mu] = bottcher_step(g, tol);
      if (k == max_k - 1)
        write_bottcher<S>(tmp.path, "test", mu, f, g);
      if (k == max_k) {
        last_mu = mu;
        g.swap(last_g);
      }
    }
    print();
  }

  // Recompute last step, and check equality
  {
    auto g = read_bottcher<S>(tfm::format("%s/g-k%d", tmp.path, max_k-1));
    const auto mu = get<1>(bottcher_step(g, tol));
    ASSERT_EQ(g.known(), last_g.known());
    ASSERT_EQ(mu, last_mu);
    for (int64_t i = 0; i < g.known(); i++)
      ASSERT_EQ(g[i], last_g[i]);
  }
}

TEST(double) { area_test<double>(1e-14); }
TEST(expansion2) { area_test<Expansion<2>>(2e-31); }

}  // namespace
}  // namespace mandelbrot
