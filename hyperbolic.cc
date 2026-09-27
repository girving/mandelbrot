// Hyperbolic component centers and areas via multiplier-map continuation
//
// For each period p, we find the centers (roots of f_c^p(0)) by Newton from Hubbard-Schleicher-Sutherland
// starting rings, then for each component continue the solution (z,c) of
//   f_c^p(z) = z,  (f_c^p)'(z) = λ
// from λ = 0 to |λ| = 1 in double, polish the boundary points in Expansion<2>, and integrate
//   area = 1/2 ∮ Im(c̄ dc) = 1/2 ∫ Re(c̄ λ c'(λ)) dθ
// with the trapezoid rule, which is spectrally accurate since c(λ) is analytic.

#include "complex.h"
#include "debug.h"
#include "expansion_arith.h"
#include "nearest.h"
#include "print.h"
#include "wall_time.h"
#include <algorithm>
#include <vector>
namespace mandelbrot {
namespace {

using std::max;
using std::min;
using std::vector;
typedef Expansion<2> E;

// OEIS A000740: number of hyperbolic components of exact period p
const int a000740[] = {0, 1, 1, 3, 6, 15, 27, 63, 120, 252, 495, 1023, 2010};

template<class S> Complex<S> cdiv(const Complex<S> a, const Complex<S> b) {
  const S d = sqr(b.r) + sqr(b.i);
  const auto n = a * conj(b);
  return Complex<S>(n.r / d, n.i / d);
}
Complex<double> to_double(const Complex<E> z) { return Complex<double>(double(z.r), double(z.i)); }
Complex<E> to_e(const Complex<double> z) { return Complex<E>(E(z.r), E(z.i)); }
Complex<double> to_double(const Complex<double> z) { return z; }
double cabs(const Complex<double> z) { return abs(z); }

// Newton step for g(c) = f_c^p(0)
template<class S> Complex<S> center_step(const Complex<S> c, const int p) {
  Complex<S> z = c, dz(1);
  for (int k = 1; k < p; k++) {
    // Far from the roots, the remaining iterations are z → z², dz → 2z dz up to O(|c|/|z|²), so
    // z_p/dz_p = z_k/(2^(p-k) dz_k).  Returning early avoids overflow for large p.
    if (cabs(to_double(z)) > 1e30) return cdiv(z, ldexp(dz, p - k));
    dz = twice(z * dz) + Complex<S>(1);
    z = sqr(z) + c;
  }
  return cdiv(z, dz);
}

struct Centers {
  vector<Complex<E>> c;  // Exact period p only
  int distinct;          // Distinct roots of f_c^p(0) over all periods dividing p
  int starts;
  double max_step;       // Size of the final Expansion<2> Newton step
};

Centers centers(const int p) {
  const int d = 1 << (p - 1);  // Degree of f_c^p(0)
  const double pi = M_PI;

  // Hubbard-Schleicher-Sutherland starting points, scaled since the roots lie in |c| <= 2
  const double logd = std::log(max(d, 2));
  const int s = max(1, int(std::ceil(0.26632 * logd)));
  const int n = max(16, int(std::ceil(8.32547 * d * logd)));
  vector<Complex<double>> roots;
  int starts = 0;
  for (int v = 1; v <= s; v++) {
    const double r = 2 * (1 + std::sqrt(2.)) * std::pow((d - 1.) / d, (2*v - 1) / (4.*s));
    for (int j = 0; j < n; j++) {
      starts++;
      const double t = 2 * pi * (j + 0.5*v) / n;
      Complex<double> c(r * std::cos(t), r * std::sin(t));
      bool ok = false;
      for (int it = 0; it < 40*d + 200; it++) {
        const auto dc = center_step(c, p);
        c -= dc;
        if (!(cabs(dc) < 1e3)) break;
        if (cabs(dc) < 1e-14) { ok = true; break; }
      }
      if (!ok) continue;
      bool dup = false;
      for (const auto& o : roots)
        if (cabs(c - o) < 1e-9) { dup = true; break; }
      if (!dup) roots.push_back(c);
    }
  }

  // Keep exact period p: reject if f_c^k(0) ≈ 0 for a proper divisor k
  Centers C;
  C.distinct = int(roots.size());
  C.starts = starts;
  C.max_step = 0;
  for (const auto& c : roots) {
    Complex<double> z = c;
    bool lower = false;
    for (int k = 1; k < p; k++) {
      if (p % k == 0 && cabs(z) < 1e-8) { lower = true; break; }
      z = sqr(z) + c;
    }
    if (lower) continue;
    auto ce = to_e(c);
    for (int it = 0; it < 3; it++) {
      const auto dc = center_step(ce, p);
      ce -= dc;
      if (it == 2) C.max_step = max(C.max_step, cabs(to_double(dc)));
    }
    C.c.push_back(ce);
  }
  return C;
}

// One Newton step for (z,c) at multiplier lam.  Returns the step size and residual; sets dc_dlam = c'(λ).
template<class S> struct NewtonResult { double step, residual; Complex<S> dc_dlam; };

template<class S> NewtonResult<S> boundary_step(Complex<S>& z, Complex<S>& c, const Complex<S> lam, const int p) {
  typedef Complex<S> C;
  C zk = z, dzz(1), dzc(0), P(1), dPz(0), dPc(0);
  for (int k = 0; k < p; k++) {
    const C P2 = twice(zk * P);
    dPz = twice(dzz * P + zk * dPz);
    dPc = twice(dzc * P + zk * dPc);
    P = P2;
    dzz = twice(zk * dzz);
    dzc = twice(zk * dzc) + C(1);
    zk = sqr(zk) + c;
  }
  const C F1 = zk - z, F2 = P - lam;
  const C a = dzz - C(1), b = dzc, cc = dPz, d = dPc;
  const C det = a * d - b * cc;
  const C dz = cdiv(d * F1 - b * F2, det), dc = cdiv(a * F2 - cc * F1, det);
  z -= dz;
  c -= dc;
  const double step = max(cabs(to_double(dz)), cabs(to_double(dc)));
  const double res = max(cabs(to_double(F1)), cabs(to_double(F2)));
  return {step, res, cdiv(a, det)};  // J (z',c') = (0,1)  =>  c' = a / det
}

struct Area {
  bool ok;
  double area_d, res_d;  // Double: area, max final residual
  E area_e;              // Expansion<2> area
  double res_e;          // Max final Expansion<2> residual
};

// Area of the period p component with center c0, using N boundary points
Area area(const Complex<E> c0, const int p, const int N, const vector<Complex<E>>& lams,
          const E pi, const int steps = 60) {
  Area A{true, 0, 0, E(0), 0};
  E sum_e(0);
  double sum_d = 0;
  const Complex<double> c0d = to_double(c0);
  for (int j = 0; j < N; j++) {
    const auto lam_e = lams[j];
    const auto lam_d = to_double(lam_e);
    // Continue radially in double
    Complex<double> z(0), c = c0d;
    NewtonResult<double> r{0, 0, Complex<double>(0)};
    for (int t = 1; t <= steps; t++) {
      const auto lam = (double(t) / steps) * lam_d;
      for (int it = 0; it < 8; it++) {
        r = boundary_step(z, c, lam, p);
        if (!(r.step < 1e3)) { A.ok = false; return A; }
        if (r.step < 1e-15) break;
      }
    }
    // One more step to get residual and c'(λ) at the final point
    // Double only needs to land in Newton's basin; convergence is judged after the Expansion<2> polish.
    r = boundary_step(z, c, lam_d, p);
    if (!(r.residual < 1e-6)) { A.ok = false; return A; }
    A.res_d = max(A.res_d, r.residual);
    sum_d += (conj(c) * r.dc_dlam * lam_d).r;

    // Polish in Expansion<2>
    auto ze = to_e(z), ce = to_e(c);
    NewtonResult<E> re{0, 0, Complex<E>(0)};
    for (int it = 0; it < 3; it++)
      re = boundary_step(ze, ce, lam_e, p);
    re = boundary_step(ze, ce, lam_e, p);  // Residual and c'(λ) at the polished point
    if (!(re.residual < 1e-20)) { A.ok = false; return A; }
    A.res_e = max(A.res_e, re.residual);
    sum_e += (conj(ce) * re.dc_dlam * lam_e).r;
  }
  A.area_d = M_PI * sum_d / N;
  A.area_e = pi * sum_e / E(int64_t(N));
  return A;
}

vector<Complex<E>> lambdas(const int N) {
  // λ_j = exp(2πi (j + 1/2) / N), offset by a half step so we never hit the root λ = 1
  vector<Complex<E>> lams(N);
  for (int j = 0; j < N; j++)
    lams[j] = nearest_twiddle<E>(2*j + 1, 2*N);
  return lams;
}

void run(const int min_p, const int max_p, const int N) {
  const auto pi = nearest_pi<E>();
  auto t0 = wall_time();
  const auto lams = lambdas(N), lams2 = lambdas(2*N);
  print("twiddles for N = %d, %d: %.3f s", N, 2*N, (wall_time() - t0).seconds());

  E cum(0);
  for (int p = min_p; p <= max_p; p++) {
    t0 = wall_time();
    const auto C = centers(p);
    const double tc = (wall_time() - t0).seconds();
    const int want = a000740[p];
    print("p = %d: %d starts, %d distinct roots (want %d), %d exact-period centers (want %d), "
          "max final center step %.1e, %.3f s",
          p, C.starts, C.distinct, 1 << (p - 1), int(C.c.size()), want, C.max_step, tc);
    slow_assert(C.distinct == (1 << (p - 1)) && int(C.c.size()) == want, "center count mismatch at p = %d", p);

    t0 = wall_time();
    double sum_d = 0, res_d = 0, res_e = 0;
    E sum_e(0);
    int failed = 0;
    for (const auto& c0 : C.c) {
      const auto A = area(c0, p, N, lams, pi);
      if (!A.ok) {
        failed++;
        print("  FAILED to converge: center %s", to_double(c0));
        continue;
      }
      sum_d += A.area_d;
      sum_e += A.area_e;
      res_d = max(res_d, A.res_d);
      res_e = max(res_e, A.res_e);
    }
    const double ta = (wall_time() - t0).seconds();

    // Trapezoid convergence check at 2N, only for the top period since it doubles the cost
    t0 = wall_time();
    string check = "skipped";
    if (p == max_p) {
      E sum2(0);
      for (const auto& c0 : C.c) {
        const auto A = area(c0, p, 2*N, lams2, pi);
        if (A.ok) sum2 += A.area_e;
      }
      check = tfm::format("%.2e", double(sum2 - sum_e));
    }
    const double t2 = (wall_time() - t0).seconds();

    cum += sum_e;
    print("  area exp2 %s", safe(sum_e));
    print("  area dbl  %.17g,  exp2 - dbl %.2e,  2N - N (exp2) %s", sum_d, double(sum_e - E(sum_d)), check);
    print("  cum from p = %d %.15g,  max resid dbl %.1e exp2 %.1e,  failed %d,  area %.3f s, 2N check %.3f s",
          min_p, double(cum), res_d, res_e, failed, ta, t2);
    if (p == 1) print("  cardioid err vs 3π/8: %.2e", double(sum_e - E(3) * pi / E(int64_t(8))));
    if (p == 2) print("  disk err vs π/16: %.2e", double(sum_e - pi / E(int64_t(16))));
  }
}

}  // namespace
}  // namespace mandelbrot

int main(const int argc, const char** argv) {
  using namespace mandelbrot;
  try {
    const int max_p = argc > 1 ? atoi(argv[1]) : 8;
    const int N = argc > 2 ? atoi(argv[2]) : 1024;
    const int min_p = argc > 3 ? atoi(argv[3]) : 1;
    slow_assert(1 <= min_p && min_p <= max_p && max_p <= 12 && N >= 16,
                "usage: %s [max_p <= 12] [N] [min_p]", argv[0]);
    const auto t0 = wall_time();
    run(min_p, max_p, N);
    print("total %.3f s", (wall_time() - t0).seconds());
    return 0;
  } catch (const std::exception& e) {
    die(e.what());
  }
}
