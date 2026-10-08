// Bulb area tests

#include "bulb.h"
#include "cutil.h"
#include "expansion_arith.h"
#include "nearest.h"
#include "tests.h"
#include <numeric>
namespace mandelbrot {
namespace {

// The 1/2 bulb is the disk |c + 1| < 1/4: area π/16, F = 1
TEST(half) {
  const auto r = bulb_areas({bulb_job(1, Complex<double>(0), 1, 2)}, BulbParams())[0];
  ASSERT_EQ(r.status, bulb_ok);
  const E2 pi = nearest_pi<E2>();
  ASSERT_LE(abs(double((r.area - pi / E2(int64_t(16))) / r.area)), 1e-29);
  ASSERT_LE(abs(double(r.F - E2(1.0))), 1e-29);
}

// Areas don't depend on the boundary resolution, polish beyond 2 steps, or (for the cardioid's) a parent given
// at period 1 vs cusp coordinates; and a few reference values
TEST(converged) {
  vector<BulbJob> jobs = {bulb_job(1, Complex<double>(0), 1, 3), bulb_job(1, Complex<double>(0), 13, 32),
                          bulb_job(1, Complex<double>(0), 2, 1449),
                          bulb_job(3, Complex<double>(-0.12256116687665358, 0.74486176661974424), 4, 13)};
  BulbParams p;
  const auto a = bulb_areas(jobs, p);
  p.N = 128; p.polish = 3;
  const auto b = bulb_areas(jobs, p);
  for (size_t i = 0; i < jobs.size(); i++) {
    ASSERT_EQ(a[i].status, bulb_ok);
    ASSERT_LE(abs(double((a[i].area - b[i].area) / a[i].area)), 1e-20);  // Near-cusp bulbs (2/1449): ~1e-23 at N = 64
  }
  ASSERT_LE(abs(double(a[0].F) - 0.9649068234500382), 1e-15);
  ASSERT_LE(abs(double(a[2].F) - 1.2393681904127829), 1e-15);
}

TEST(cuda_matches_cpu) {
  IF_CUDA({
    vector<BulbJob> jobs;
    for (int q = 2; q <= 40; q++)
      for (int p = 1; p < q; p++)
        if (std::gcd(p, q) == 1) {
          if (2 * p <= q) jobs.push_back(bulb_job(1, Complex<double>(0), p, q));
          if (q <= 20) jobs.push_back(bulb_job(3, Complex<double>(-0.12256116687665358, 0.74486176661974424), p, q));
        }
    BulbParams p;
    const auto cpu = bulb_areas(jobs, p);
    p.cuda = true;
    const auto gpu = bulb_areas(jobs, p);
    for (size_t i = 0; i < jobs.size(); i++) {
      ASSERT_EQ(cpu[i].status, gpu[i].status);
      if (cpu[i].status != bulb_ok) continue;
      ASSERT_EQ(cpu[i].area.x[0], gpu[i].area.x[0]);
      ASSERT_EQ(cpu[i].area.x[1], gpu[i].area.x[1]);
      ASSERT_EQ(cpu[i].F.x[0], gpu[i].F.x[0]);
    }
    print("cuda_matches_cpu: %zu bulbs identical", jobs.size());
  })
}

}  // namespace
}  // namespace mandelbrot
