// Areas of the satellite bulbs of one hyperbolic component, for the renormalization cascade.
//
// The parent W has period P and multiplier map c_W(λ) (λ the multiplier of its attracting P-cycle).  Its p/q
// satellite is the period qP component whose root is c_W(e^{2πip/q}).  For each requested p/q we
//   - continue W's cycle from its center to λ0 = e^{2πip/q}, giving the root c_r and c_W'(λ0);
//   - find the child's center by Newton on f_c^{qP}(0) = 0 from c_r + λ0 c_W'(λ0) / q^2 (exact for the cardioid's
//     1/2 bulb, and the right scale in general), and check it is the child: exact period qP, and continuing its own
//     multiplier to 1 lands on c_r;
//   - trace the child's boundary c(μ), |μ| = 1, by continuation around the circle, and integrate
//     area = 1/2 ∫ Re((c - c0)‾ μ c'(μ)) dθ (relative to the center c0, against cancellation) by the trapezoid
//     rule (spectral: c is analytic), at N and 2N points.
// Normalized area F = area q^4 / (π |c_W'(λ0)|^2): about 1 for every bulb if bulbs are universal in the parent's
// multiplier coordinate.
// For the cardioid (P = 1) it works in cusp coordinates ζ = z - 1/2, δ = c - 1/4, where bulbs near c = 1/4 keep
// full relative precision (in c itself, a bulb of radius 1e-9 near 0.25 loses most of its digits).
// Failures print the stage ("failed parent/center/period/area").  $BULB_TOL (default 1e-11) is the final Newton
// step accepted when roundoff stops it short of 1e-15: bulbs with large interior digits (hundreds of near-parabolic
// passes per cycle) need 1e-9.  $BULB_STEPS, $BULB_SUBSTEPS set the continuation steps (64, 4).
//   ./build/release/bulb_areas P c_re c_im N < list
// reads lines "p q" and prints "p q center_re center_im area F |c_W'|^2 (2N - N)/area root_err/size".  With
// $BULB_EXP=1 boundary points are polished in Expansion<2>, conv is the N/2 subrule's relative change, and two more
// columns give the low parts of area and F (area = col5 + col10, F = col6 + col11).  With $BULB_PROFILE=1 (double
// mode) ten more columns: the parent's κ = λ c_W''/c_W' at the root (re, im) and the child's Fourier coefficients
// h_1..h_4 of log|c'(e^{iθ})|^2 on its boundary (re, im each).  The
// parent's (c_re, c_im) is its center (0 0 for the cardioid).  CPU threads: $MANDELBROT_THREADS.
#include "complex.h"
#include "expansion_arith.h"
#include "nearest.h"
#include <atomic>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <vector>
typedef std::complex<double> Cx;
using mandelbrot::Complex;
typedef mandelbrot::Expansion<2> E;
typedef Complex<E> CE;

// The map is x ↦ x^2 + s x + c, with s = 0 (z ↦ z^2 + c) or s = 1: cusp coordinates ζ = z - 1/2, δ = c - 1/4 for
// the cardioid's bulbs, where ζ ↦ ζ^2 + ζ + δ keeps full relative precision in bulbs near the cusp (tiny δ and
// orbits lingering near ζ = 0).  The critical point is x = -s/2.
static double shift = 0;
// Continuation steps: radial to the boundary, and substeps between boundary points ($BULB_STEPS, $BULB_SUBSTEPS)
static int radial_steps = 64, substeps = 4;
// Newton accepts a final step below this when it can't reach 1e-15 ($BULB_TOL)
static double accept = 1e-11;

// Newton for (z, c) with f_c^n(z) = z, (f_c^n)'(z) = mu; also returns dc/dμ.  True if converged.
static bool solve(const int n, const Cx mu, Cx& z, Cx& c, Cx* dcdmu = nullptr) {
  double last = INFINITY;
  for (int it = 0; it < 40; it++) {
    Cx x = z, xz = 1, xc = 0, xzz = 0, xzc = 0;
    for (int i = 0; i < n; i++) {
      const Cx nxzz = 2.0 * xz * xz + (2.0 * x + shift) * xzz, nxzc = 2.0 * xc * xz + (2.0 * x + shift) * xzc;
      xzz = nxzz; xzc = nxzc;
      const Cx df = 2.0 * x + shift;
      xc = df * xc + 1.0;
      xz = df * xz;
      x = x * x + shift * x + c;
    }
    const Cx F1 = x - z, F2 = xz - mu, a = xz - 1.0, b = xc, d = xzz, e = xzc, det = a * e - b * d;
    if (dcdmu) *dcdmu = a / det;
    const Cx dz = (F1 * e - b * F2) / det, dc = (a * F2 - d * F1) / det;
    z -= dz; c -= dc;
    last = std::abs(dz) + std::abs(dc);
    if (!(last < 1)) return false;
    if (last < 1e-15 * (1 + std::abs(c))) return true;
  }
  return last < accept;  // Roundoff limits the final steps (near-parabolic passes amplify it)
}

// Continue (z, c) radially from μ = 0 (z at the critical point, c at the center) to mu
static bool radial(const int n, const Cx mu, Cx& z, Cx& c, const int steps, Cx* dcdmu = nullptr) {
  for (int s = 1; s <= steps; s++)
    if (!solve(n, mu * (double(s) / steps), z, c, s == steps ? dcdmu : nullptr)) return false;
  return true;
}

// Double-double (Expansion<2>) polish, for areas to ~1e-25 relative ($BULB_EXP=1)
static CE to_e(const Cx z) { return CE(E(z.real()), E(z.imag())); }
static Cx to_d(const CE z) { return Cx(double(z.r), double(z.i)); }
static CE cdiv(const CE a, const CE b) {
  const E d = sqr(b.r) + sqr(b.i);
  const CE n = a * conj(b);
  return CE(n.r / d, n.i / d);
}
static const CE shift_e() { return CE(E(shift), E(0.0)); }

// One Newton step in E for f^n(z) = z, (f^n)'(z) = mu; returns dc/dμ
static CE step_e(const int n, const CE mu, CE& z, CE& c) {
  const CE s = shift_e(), one(1), two(2);
  CE x = z, xz(1), xc(0), xzz(0), xzc(0);
  for (int i = 0; i < n; i++) {
    const CE df = twice(x) + s;
    const CE nxzz = twice(sqr(xz)) + df * xzz, nxzc = twice(xc * xz) + df * xzc;
    xzz = nxzz; xzc = nxzc;
    xc = df * xc + one;
    xz = df * xz;
    x = sqr(x) + s * x + c;
  }
  const CE F1 = x - z, F2 = xz - mu, a = xz - one, b = xc, d = xzz, e = xzc, det = a * e - b * d;
  z -= cdiv(F1 * e - b * F2, det);
  c -= cdiv(a * F2 - d * F1, det);
  return cdiv(a, det);
}

// Area of the period n component with center c0 in E, at N boundary points (double continuation, E polish), and
// the relative change from the N/2 subrule (every other point: a rotated trapezoid rule)
static E area_e(const int n, const CE c0, const int N, double& conv, bool& ok) {
  const Cx c0d = to_d(c0);
  Cx z = -shift / 2, c = c0d;
  ok = radial(n, std::polar(1.0, M_PI / N), z, c, radial_steps);
  E sum(0.0), half_sum(0.0);
  for (int j = 0; j < N && ok; j++) {
    const Cx mu = std::polar(1.0, M_PI * (2 * j + 1) / N);
    if (j) {
      const Cx prev = std::polar(1.0, M_PI * (2 * j - 1) / N);
      for (int s = 1; s <= substeps && ok; s++)
        ok = solve(n, prev * std::polar(1.0, 2 * M_PI * s / (double(substeps) * N)), z, c);
    }
    ok = ok && solve(n, mu, z, c);
    if (!ok) break;
    const CE mu_e = mandelbrot::nearest_twiddle<E>(2 * j + 1, 2 * N);
    CE ze = to_e(z), ce = to_e(c), dc;
    for (int it = 0; it < 3; it++) dc = step_e(n, mu_e, ze, ce);
    const E t = (conj(ce - c0) * mu_e * dc).r;
    sum += t;
    if (j % 2 == 0) half_sum += t;
  }
  const E A = sum / E(int64_t(N)), A2 = half_sum / E(int64_t(N / 2));
  conv = double((A - A2) / A);
  return A;  // Times π by the caller
}

// Area of the period n component with center c0, at N boundary points; optionally log|c'(μ_j)|^2 at the points
static double area(const int n, const Cx c0, const int N, bool& ok, std::vector<double>* prof = nullptr) {
  Cx z = -shift / 2, c = c0;
  ok = radial(n, std::polar(1.0, M_PI / N), z, c, radial_steps);
  double sum = 0;
  for (int j = 0; j < N && ok; j++) {
    const Cx mu = std::polar(1.0, M_PI * (2 * j + 1) / N);
    Cx dc;
    // March around the circle in substeps
    if (j) {
      const Cx prev = std::polar(1.0, M_PI * (2 * j - 1) / N);
      for (int s = 1; s <= substeps && ok; s++)
        ok = solve(n, prev * std::polar(1.0, 2 * M_PI * s / (double(substeps) * N)), z, c);
    }
    ok = ok && solve(n, mu, z, c, &dc);
    sum += (std::conj(c - c0) * mu * dc).real();  // Relative to the center: bulbs are tiny
    if (prof) prof->push_back(std::log(std::norm(dc)));
  }
  return M_PI * sum / N;
}

int main(int argc, char** argv) {
  if (argc < 5) { fprintf(stderr, "usage: bulb_areas P c_re c_im N < list\n"); return 1; }
  const int P = atoi(argv[1]), N = atoi(argv[4]);
  if (getenv("BULB_STEPS")) radial_steps = atoi(getenv("BULB_STEPS"));
  if (getenv("BULB_SUBSTEPS")) substeps = atoi(getenv("BULB_SUBSTEPS"));
  if (getenv("BULB_TOL")) accept = atof(getenv("BULB_TOL"));
  const bool exp2 = getenv("BULB_EXP") && atoi(getenv("BULB_EXP"));
  const E pi_e = mandelbrot::nearest_pi<E>();
  Cx center(atof(argv[2]), atof(argv[3]));
  if (P == 1) {  // Cusp coordinates for the cardioid
    shift = 1;
    center -= 0.25;
  }
  const Cx crit = -shift / 2;
  std::vector<std::pair<int, int>> jobs;
  int p, q;
  while (scanf("%d %d", &p, &q) == 2) jobs.push_back({p, q});
  struct Out { Cx cc; double A, F, w, conv, root_err; bool ok; const char* why; E Ae, Fe; Cx kappa, h[4]; };
  const bool profile = getenv("BULB_PROFILE") && atoi(getenv("BULB_PROFILE"));
  std::vector<Out> out(jobs.size());
  std::atomic<size_t> next(0);
  const char* te = getenv("MANDELBROT_THREADS");
  const int T = te ? atoi(te) : int(std::thread::hardware_concurrency());
  std::vector<std::thread> pool;
  for (int t = 0; t < T; t++)
    pool.emplace_back([&]() {
      for (size_t i; (i = next++) < jobs.size();) {
        const int p = jobs[i].first, q = jobs[i].second, n = q * P;
        Out& o = out[i];
        o.ok = false;
        o.why = "parent";
        const Cx l0 = std::polar(1.0, 2 * M_PI * p / q);
        Cx zr = crit, cr = center, dW;
        if (!radial(P, l0, zr, cr, 128, &dW)) continue;
        // Child center: Newton on g(c) = f_c^n(0)
        Cx c = cr + l0 * dW / double(q * q);
        bool conv = false;
        for (int it = 0; it < 200 && !conv; it++) {
          Cx x = crit, dx = 0;
          for (int k = 0; k < n; k++) { dx = (2.0 * x + shift) * dx + 1.0; x = x * x + shift * x + c; }
          const Cx step = (x - crit) / dx;
          c -= step;
          conv = std::abs(step) < 1e-16 * (1 + std::abs(c));
          if (!(std::abs(step) < 1)) break;
        }
        o.why = "center";
        if (!conv) continue;
        bool exact = true;
        Cx x = crit;
        for (int k = 1; k < n && exact; k++) {
          x = x * x + shift * x + c;
          if (n % k == 0 && std::abs(x - crit) < 1e-8) exact = false;
        }
        o.why = "period";
        if (!exact) continue;
        // Root check: continue the child's own multiplier to 1 - 1e-4 (closer, the system degenerates)
        Cx zc = crit, cc = c;
        const double size = std::abs(dW) / (q * q);
        o.root_err = radial(n, 1 - 1e-4, zc, cc, 128) ? std::abs(cc - cr) / size : INFINITY;
        if (exp2) {
          // Polish the center in E, then the area at N points
          CE ce = to_e(c), crit_e = to_e(crit);
          const CE sh = shift_e(), one(1);
          for (int it = 0; it < 3; it++) {
            CE x = crit_e, dx(0);
            for (int k = 0; k < n; k++) { dx = (twice(x) + sh) * dx + one; x = sqr(x) + sh * x + ce; }
            ce -= cdiv(x - crit_e, dx);
          }
          bool ok;
          double conv;
          const E A = pi_e * area_e(n, ce, N, conv, ok);
          o.why = "area";
          if (!ok) continue;
          // |c_W'(λ0)|^2: exact for the cardioid, |1 - λ0|^2 / 4; else from the double continuation
          E w(std::norm(dW));
          if (P == 1) {
            const CE l0e = mandelbrot::nearest_twiddle<E>(p, q), d = CE(1) - l0e;
            w = (sqr(d.r) + sqr(d.i)) / E(int64_t(4));
          }
          const E q4 = E(int64_t(q) * q * q * q);
          o.cc = to_d(ce) + 0.25 * shift; o.Ae = A; o.Fe = A * q4 / (pi_e * w);
          o.A = double(A); o.F = double(o.Fe); o.w = double(w); o.conv = conv; o.ok = true;
          continue;
        }
        bool ok1, ok2;
        std::vector<double> prof;
        const double A1 = area(n, c, N, ok1), A2 = area(n, c, 2 * N, ok2, profile ? &prof : nullptr);
        if (profile && ok2) {
          // Parent: κ = λ c_W''(λ0) / c_W'(λ0) by central differences along the ray; child: harmonics of log|c'|^2
          const double h = 1e-4;
          Cx zp = crit, cp = center, dp, zm = crit, cm = center, dm;
          radial(P, l0 * (1 + h), zp, cp, 128, &dp);
          radial(P, l0 * (1 - h), zm, cm, 128, &dm);
          o.kappa = (dp - dm) / (2 * h * dW);
          const int M = int(prof.size());
          for (int k = 1; k <= 4; k++) {
            Cx sk = 0;
            for (int j = 0; j < M; j++) sk += prof[j] * std::polar(1.0, -k * M_PI * (2 * j + 1) / M);
            o.h[k - 1] = sk / double(M);
          }
        }
        o.why = "area";
        if (!ok1 || !ok2) continue;
        o.cc = c + 0.25 * shift; o.A = A2; o.w = std::norm(dW);
        o.F = A2 * double(q) * q * q * q / (M_PI * o.w);
        o.conv = (A2 - A1) / A2;
        o.ok = true;
      }
    });
  for (auto& th : pool) th.join();
  for (size_t i = 0; i < jobs.size(); i++) {
    const auto& o = out[i];
    if (o.ok && exp2)  // Two more columns: the E low parts of area and F
      printf("%d %d %.17g %.17g %.17g %.17g %.17g %.1e %.1e %.17g %.17g\n", jobs[i].first, jobs[i].second,
             o.cc.real(), o.cc.imag(), o.Ae.x[0], o.Fe.x[0], o.w, o.conv, o.root_err, o.Ae.x[1], o.Fe.x[1]);
    else if (o.ok && profile) {
      printf("%d %d %.17g %.17g %.17g %.15f %.17g %.1e %.1e %.12g %.12g", jobs[i].first, jobs[i].second, o.cc.real(),
             o.cc.imag(), o.A, o.F, o.w, o.conv, o.root_err, o.kappa.real(), o.kappa.imag());
      for (int k = 0; k < 4; k++) printf(" %.12g %.12g", o.h[k].real(), o.h[k].imag());
      printf("\n");
    } else if (o.ok)
      printf("%d %d %.17g %.17g %.17g %.15f %.17g %.1e %.1e\n", jobs[i].first, jobs[i].second, o.cc.real(),
             o.cc.imag(), o.A, o.F, o.w, o.conv, o.root_err);
    else
      printf("%d %d failed %s\n", jobs[i].first, jobs[i].second, o.why);
  }
}
