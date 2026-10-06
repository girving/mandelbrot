// Tests for filled Julia set areas via the area transfer operator

#include "julia.h"
#include "expansion_arith.h"
#include "tests.h"
#include <cmath>
#include <random>
namespace mandelbrot {
namespace {

using std::abs;

TEST(trig_row) {
  // Barycentric interpolation on nt odd equispaced angles reproduces trig polynomials of degree ≤ (nt-1)/2
  const int nt = 9, K = 4;
  std::mt19937 rng(3);
  std::uniform_real_distribution<double> u(-1, 1);
  vector<double> ac(K + 1), as(K + 1), rx(nt), ry(nt), f(nt), row(nt);
  for (int k = 0; k <= K; k++) { ac[k] = u(rng); as[k] = u(rng); }
  const auto F = [&](const double t) { double s = 0; for (int k = 0; k <= K; k++) s += ac[k] * cos(k * t) + as[k] * sin(k * t); return s; };
  for (int j = 0; j < nt; j++) { rx[j] = cos(M_PI * j / nt); ry[j] = sin(M_PI * j / nt); f[j] = F(2 * M_PI * j / nt); }
  for (int i = 0; i < 100; i++) {
    const double t = 4 * M_PI * u(rng);
    trig_row(rx, ry, cos(t / 2), sin(t / 2), row.data());
    double s = 0;
    for (int j = 0; j < nt; j++) s += row[j] * f[j];
    ASSERT_LE(abs(s - F(t)), 1e-13);
  }
}

TEST(cheb_row) {
  // Barycentric interpolation at first kind Chebyshev points reproduces polynomials of degree < n
  const int n = 7;
  vector<double> x(n), w(n), f(n), row(n);
  const auto F = [](const double s) { return ((((0.3 * s - 1) * s + 2) * s - 0.5) * s + 1) * s * s + 0.25; };
  for (int k = 0; k < n; k++) {
    x[k] = 1 + 2 * cos((2 * k + 1) * M_PI / (2 * n));
    w[k] = (k & 1 ? -1 : 1) * sin((2 * k + 1) * M_PI / (2 * n));
    f[k] = F(x[k]);
  }
  for (const double s : {-1.0, 0.3, 1.0, 2.7, 3.0}) {
    cheb_row(x, w, s, row.data());
    double v = 0;
    for (int k = 0; k < n; k++) v += row[k] * f[k];
    ASSERT_LE(abs(v - F(s)), 1e-12);
  }
}

TEST(disk) {
  // c = 0: K is the closed unit disk, area π
  JuliaParams p;
  p.nr = 24; p.nt = 5;
  const double e = julia_area<double>(p).area - M_PI;
  p.nr = 40;
  const auto r2 = julia_area<Expansion<2>>(p);
  const auto r3 = julia_area<Expansion<3>>(p);
  const Expansion<3> pi(string("3.14159265358979323846264338327950288419716939937510582097494459"));
  const double e2 = double(Expansion<3>(r2.area.x[0], r2.area.x[1], 0.0, nonoverlap) - pi), e3 = double(r3.area - pi);
  print("  π errors: double %.3g, Expansion<2> %.3g, Expansion<3> %.3g", e, e2, e3);
  ASSERT_LE(abs(e), 1e-14);
  ASSERT_LE(abs(e2), 1e-30);
  ASSERT_LE(abs(e3), 1e-44);
}

TEST(small_c) {
  // c = -0.2 against the numpy prototype (scratch/julia/transfer.py: 3.0673885651, converged to ~2e-11), and
  // Expansion<2> against Expansion<3> on the same grid
  JuliaParams p;
  p.cx = -0.2;
  p.nr = 40; p.nt = 121;
  const auto rd = julia_area<double>(p);
  const auto r2 = julia_area<Expansion<2>>(p);
  const auto r3 = julia_area<Expansion<3>>(p);
  const double d23 = double(Expansion<3>(r2.area.x[0], r2.area.x[1], 0.0, nonoverlap) - r3.area);
  print("  c = -0.2: area %.15f, residual %.3g, %d refinements, %d GMRES iterations, %.2f s; dd - td %.3g",
        rd.area, rd.residual, rd.refinements, rd.gmres_iters, rd.secs, d23);
  ASSERT_LE(abs(rd.area - 3.0673885651), 1e-9);
  ASSERT_LE(abs(d23), 1e-28);
}

TEST(grade) {
  // Angular grading changes the discretization, not the answer: near the cusp it needs far fewer points
  JuliaParams p;
  p.cx = 0.25 - std::ldexp(1, -9);
  p.nr = 110; p.nt = 331;
  const double uniform = julia_area<double>(p).area;
  p.nr = 60; p.nt = 121; p.grade = 0.6;
  const double graded = julia_area<double>(p).area;
  print("  c = 1/4 - 2^-9: uniform 110x331 %.16f, graded 60x121 %.16f", uniform, graded);
  ASSERT_LE(abs(uniform - graded), 1e-13);
}

TEST(parabolic_preconditioner) {
  // Near c = 1/4, (1 - L+)^-1 over the branch fixing q leaves the answer alone and cuts GMRES iterations
  JuliaParams p;
  p.cx = 0.25 - std::ldexp(1, -12);
  p.nr = 80; p.nt = 181; p.grade = 0.744;
  p.pre = 0;
  const auto plain = julia_area<double>(p);
  p.pre = 1;
  const auto pre = julia_area<double>(p);
  print("  c = 1/4 - 2^-12: plain %.16f (%d GMRES iterations, %.2f s), preconditioned %.16f (%d, %.2f s)",
        plain.area, plain.gmres_iters, plain.secs, pre.area, pre.gmres_iters, pre.secs);
  ASSERT_LE(abs(plain.area - pre.area), 1e-13);
  ASSERT_LE(4 * pre.gmres_iters, plain.gmres_iters);
}

TEST(eigenvalue) {
  // L's leading eigenvalue e^{P(2)}: exactly 1/2 at c = 0 (L^n 1 = 2^-n |z|^-2 (1 - 2^-n)), and at real c ≥ the
  // variational bound |f'(q)|^-2 from the delta measure at the repelling fixed point q = (1 + √(1 - 4c))/2
  JuliaParams p;
  p.eig = true;
  p.nr = 40; p.nt = 121;
  ASSERT_LE(abs(julia_area<double>(p).rho - 0.5), 1e-12);
  for (const double c : {0.1, 0.2, 0.24}) {
    p.cx = c;
    const double rho = julia_area<double>(p).rho, q = (1 + sqrt(1 - 4 * c)) / 2, bound = 1 / (4 * q * q);
    print("  c = %g: rho %.12f, |f'(q)|^-2 %.12f", c, rho, bound);
    ASSERT_LE(bound, rho);
    ASSERT_LE(rho, 1);
  }
}

}  // namespace
}  // namespace mandelbrot
