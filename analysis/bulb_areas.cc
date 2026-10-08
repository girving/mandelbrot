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
//   ./build/release/bulb_areas P c_re c_im N < list
// reads lines "p q" and prints "p q center_re center_im area F |c_W'|^2 (2N - N)/area root_err/size".  The
// parent's (c_re, c_im) is its center (0 0 for the cardioid).  CPU threads: $MANDELBROT_THREADS.
#include <atomic>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <vector>
typedef std::complex<double> Cx;

// The map is x ↦ x^2 + s x + c, with s = 0 (z ↦ z^2 + c) or s = 1: cusp coordinates ζ = z - 1/2, δ = c - 1/4 for
// the cardioid's bulbs, where ζ ↦ ζ^2 + ζ + δ keeps full relative precision in bulbs near the cusp (tiny δ and
// orbits lingering near ζ = 0).  The critical point is x = -s/2.
static double shift = 0;

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
  return last < 1e-11;
}

// Continue (z, c) radially from μ = 0 (z at the critical point, c at the center) to mu
static bool radial(const int n, const Cx mu, Cx& z, Cx& c, const int steps, Cx* dcdmu = nullptr) {
  for (int s = 1; s <= steps; s++)
    if (!solve(n, mu * (double(s) / steps), z, c, s == steps ? dcdmu : nullptr)) return false;
  return true;
}

// Area of the period n component with center c0, at N boundary points
static double area(const int n, const Cx c0, const int N, bool& ok) {
  Cx z = -shift / 2, c = c0;
  ok = radial(n, std::polar(1.0, M_PI / N), z, c, 64);
  double sum = 0;
  for (int j = 0; j < N && ok; j++) {
    const Cx mu = std::polar(1.0, M_PI * (2 * j + 1) / N);
    Cx dc;
    // March around the circle in substeps
    if (j) {
      const Cx prev = std::polar(1.0, M_PI * (2 * j - 1) / N);
      for (int s = 1; s <= 4 && ok; s++) ok = solve(n, prev * std::polar(1.0, 2 * M_PI * s / (4.0 * N)), z, c);
    }
    ok = ok && solve(n, mu, z, c, &dc);
    sum += (std::conj(c - c0) * mu * dc).real();  // Relative to the center: bulbs are tiny
  }
  return M_PI * sum / N;
}

int main(int argc, char** argv) {
  if (argc < 5) { fprintf(stderr, "usage: bulb_areas P c_re c_im N < list\n"); return 1; }
  const int P = atoi(argv[1]), N = atoi(argv[4]);
  Cx center(atof(argv[2]), atof(argv[3]));
  if (P == 1) {  // Cusp coordinates for the cardioid
    shift = 1;
    center -= 0.25;
  }
  const Cx crit = -shift / 2;
  std::vector<std::pair<int, int>> jobs;
  int p, q;
  while (scanf("%d %d", &p, &q) == 2) jobs.push_back({p, q});
  struct Out { Cx cc; double A, F, w, conv, root_err; bool ok; };
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
        if (!conv) continue;
        bool exact = true;
        Cx x = crit;
        for (int k = 1; k < n && exact; k++) {
          x = x * x + shift * x + c;
          if (n % k == 0 && std::abs(x - crit) < 1e-8) exact = false;
        }
        if (!exact) continue;
        // Root check: continue the child's own multiplier to 1 - 1e-4 (closer, the system degenerates)
        Cx zc = crit, cc = c;
        const double size = std::abs(dW) / (q * q);
        o.root_err = radial(n, 1 - 1e-4, zc, cc, 128) ? std::abs(cc - cr) / size : INFINITY;
        bool ok1, ok2;
        const double A1 = area(n, c, N, ok1), A2 = area(n, c, 2 * N, ok2);
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
    if (o.ok)
      printf("%d %d %.17g %.17g %.17g %.15f %.17g %.1e %.1e\n", jobs[i].first, jobs[i].second, o.cc.real(),
             o.cc.imag(), o.A, o.F, o.w, o.conv, o.root_err);
    else
      printf("%d %d failed\n", jobs[i].first, jobs[i].second);
  }
}
