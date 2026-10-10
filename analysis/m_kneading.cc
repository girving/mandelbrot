// Exact kneading of Mandelbrot centres by their own rays: the first kneading 0 (the first entry after 1 of the internal
// address: a component in the wake of the cardioid's p/q bulb has it at q).
//
//   ./build/release/m_kneading [threads] < "name c_re c_im P"   →   "name c_re c_im P first0 margin" (or "name failed")
//
// For the centre c of period P (polished by Newton): the root ρ of the critical point's Fatou component is the landing
// point of its internal ray 0 (f^P ≈ A z² near 0, Böttcher φ ≈ A z, so the ray at radius t is near t/A; walked out by
// inverse iteration of f^P, then Newton on f^P(z) = z); the external ray landing at ρ is traced as a field line of the
// Green's function upward from just outside ρ; with the internal rays 0 → ±ρ and the mirror ray (-z: the ray of angle
// θ/2 + 1/2) it is the kneading partition.  ν_i = 1 iff f^i(c) is on the side of c; first0 = the first i < P with
// ν_i = 0 (P if none), margin = the smallest distance of an orbit point f^i(c), i < first0, to the partition.
#include <atomic>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <thread>
#include <vector>
using std::complex;
using std::vector;
typedef complex<double> C;

static bool center(C& c, const int P) {
  for (int it = 0; it < 100; it++) {
    C z = 0, dz = 0;
    for (int j = 0; j < P; j++) { dz = 2.0 * z * dz + 1.0; z = z * z + c; }
    const C step = z / dz;
    c -= step;
    if (std::abs(step) < 1e-15 * (1 + std::abs(c))) return true;
    if (!(std::abs(step) < 1)) return false;
  }
  return false;
}
// f^P and its derivative at z
static void fP(const C c, const int P, const C z0, C& z, C& dz) {
  z = z0; dz = 1;
  for (int j = 0; j < P; j++) { dz = 2.0 * z * dz; z = z * z + c; }
}
// Green's function gradient: G = 2^-n log|f^n|, returns false if z does not escape in nmax steps
// The Green's function near the Julia set underflows (G = 2^-n log|f^n| with n in the thousands), so return log G and
// the gradient step scale: G/|∇G| = |f^n| log|f^n| / |(f^n)'| (2^-n cancels) and the gradient's direction
static bool green(const C c, const C z0, double& logG, double& scale, C& dir, const int nmax) {
  C z = z0, dz = 1;
  for (int n = 1; n <= nmax; n++) {
    dz = 2.0 * z * dz; z = z * z + c;
    if (std::norm(dz) > 1e200) return false;
    if (std::norm(z) > 1e20) {
      const double L = std::log(std::abs(z));
      logG = std::log(L) - n * std::log(2.0);
      scale = std::abs(z) * L / std::abs(dz);
      const C g = std::conj(dz / z);
      dir = g / std::abs(g);
      return true;
    }
  }
  return false;
}

int main(int argc, char** argv) {
  const int threads = argc > 1 ? atoi(argv[1]) : 2;
  struct Job { std::string name; C c; int P; };
  vector<Job> jobs;
  char name[256];
  double cr, ci;
  int P;
  while (scanf("%255s %lf %lf %d", name, &cr, &ci, &P) == 4) jobs.push_back({name, C(cr, ci), P});
  vector<std::string> out(jobs.size());
  std::atomic<int> next(0);
  vector<std::thread> pool;
  for (int t = 0; t < threads; t++)
    pool.emplace_back([&]() {
      for (int i; (i = next++) < int(jobs.size());) {
        const auto& J = jobs[i];
        char line[512];
        snprintf(line, sizeof(line), "%s failed", J.name.c_str());
        out[i] = line;
        C c = J.c;
        const int P = J.P;
        if (!center(c, P)) continue;
        // A = (f^P)''(0)/2
        C z = 0, a = 1, b = 0;
        for (int j = 0; j < P; j++) { b = 2.0 * (a * a + z * b); a = 2.0 * z * a; z = z * z + c; }
        const C A = b / 2.0;
        // internal ray 0: points at radius t → √t, each solving f^P(w) = previous
        double rad = 1e-3, rprev = 0;
        C prev = 0, cur = rad / A;
        vector<C> iray = {C(0), cur};
        bool ok = true;
        for (int k = 0; k < 60 && rad < 1 - 1e-12; k++) {
          const double rn = std::sqrt(rad);
          C w = cur + (cur - prev) * ((rn - rad) / (rad - rprev));
          bool cv = false;
          for (int it = 0; it < 60; it++) {
            C fz, dfz;
            fP(c, P, w, fz, dfz);
            const C res = fz - cur;
            w -= res / dfz;
            if (std::abs(res) < 1e-13 * (1 + std::abs(cur))) { cv = true; break; }
          }
          if (!cv) { ok = false; break; }
          prev = cur; rprev = rad; cur = w; rad = rn; iray.push_back(cur);
        }
        if (!ok) continue;
        C rho = cur;
        for (int it = 0; it < 60; it++) {
          C fz, dfz;
          fP(c, P, rho, fz, dfz);
          const C res = fz - rho;
          rho -= res / (dfz - 1.0);
          if (std::abs(res) < 1e-14 * (1 + std::abs(rho))) break;
        }
        iray.push_back(rho);
        // the external ray landing at ρ: a field line of G up from ρ + ε d, d = outward along the internal ray
        const C d = (rho - iray[iray.size() - 3]) / std::abs(rho - iray[iray.size() - 3]);
        const int nmax = 200000;
        vector<C> ray;
        bool found = false;
        for (const double eps : {1e-9, 1e-8, 1e-7, 1e-6})
          for (const double ang : {0.0, 0.3, -0.3, 0.7, -0.7, 1.2, -1.2}) {
            if (found) break;
            const C z0 = rho + eps * std::abs(rho) * d * std::polar(1.0, ang);
            double lg, sc;
            C dr;
            if (!green(c, z0, lg, sc, dr, nmax)) continue;
            ray.assign(1, z0);
            C zz = z0;
            bool good = true;
            double eta = 0.02;   // relative increase of G per step, halved when a trial point does not escape
            for (int st = 0; st < 400000 && std::abs(zz) < 12; st++) {
              double l1, l2, l3, s1, s2, s3;
              C d1, d2, d3;
              if (!green(c, zz, l1, s1, d1, nmax)) { good = false; break; }
              const C h1 = eta * s1 * d1;
              if (!green(c, zz + 0.5 * h1, l2, s2, d2, nmax)) { eta *= 0.5; if (eta < 1e-8) { good = false; break; } continue; }
              const C znew = zz + eta * s2 * d2 * std::exp(l1 - l2);   // G1/|∇G(mid)| = s2 G1/G2
              if (!green(c, znew, l3, s3, d3, nmax) || !(l3 > l1)) { eta *= 0.5; if (eta < 1e-8) { good = false; break; } continue; }
              zz = znew;
              ray.push_back(zz);
              eta = std::min(0.05, eta * 1.5);
            }
            if (getenv("MK_DEBUG")) fprintf(stderr, "  trace from eps %g ang %g: %zu points, |z| %.3g, good %d\n", eps, ang, ray.size(), std::abs(zz), good);
            if (good && std::abs(zz) >= 12) found = true;
          }
        if (getenv("MK_DEBUG")) {
          C fz, dfz;
          fP(c, P, rho, fz, dfz);
          fprintf(stderr, "%s: c %.12f%+.12fi A %.4g%+.4gi |ρ| %.6f ρ %.8f%+.8fi |f^P'(ρ)| %.4g, internal ray %zu pts, found %d\n", J.name.c_str(),
                  c.real(), c.imag(), A.real(), A.imag(), std::abs(rho), rho.real(), rho.imag(), std::abs(dfz), iray.size(), found);
          for (const double eps : {1e-9, 1e-7, 1e-5})
            for (int q = 0; q < 8; q++) {
              const C z0 = rho + eps * std::abs(rho) * d * std::polar(1.0, 2 * M_PI * q / 8);
              double lg, sc;
              C dr;
              fprintf(stderr, "  eps %g ang %d/8 escapes %d\n", eps, q, green(c, z0, lg, sc, dr, nmax));
            }
        }
        if (!found) { snprintf(line, sizeof(line), "%s failed ray", J.name.c_str()); out[i] = line; continue; }
        // partition polygon: far end of the ray → … → ρ → internal ray → 0 → -internal ray → -ρ → -ray → far, closed by
        // an arc of radius 200
        vector<C> poly;
        const C far1 = ray.back() * (200.0 / std::abs(ray.back()));
        poly.push_back(far1);
        for (int k = int(ray.size()) - 1; k >= 0; k--) poly.push_back(ray[k]);
        for (int k = int(iray.size()) - 1; k >= 0; k--) poly.push_back(iray[k]);
        for (size_t k = 1; k < iray.size(); k++) poly.push_back(-iray[k]);
        for (size_t k = 0; k < ray.size(); k++) poly.push_back(-ray[k]);
        poly.push_back(-far1);
        const double a0 = std::arg(-far1);
        for (int k = 1; k < 360; k++) poly.push_back(std::polar(200.0, a0 + M_PI * k / 360));
        const auto inside = [&](const C q) {
          bool in = false;
          for (size_t e = 0; e < poly.size(); e++) {
            const C u = poly[e], v = poly[(e + 1) % poly.size()];
            if ((u.imag() > q.imag()) != (v.imag() > q.imag())) {
              const double x = u.real() + (q.imag() - u.imag()) * (v.real() - u.real()) / (v.imag() - u.imag());
              if (x > q.real()) in = !in;
            }
          }
          return in;
        };
        const auto dist = [&](const C q) {
          double dm = INFINITY;
          for (size_t e = 0; e < poly.size(); e++) {
            const C u = poly[e], v = poly[(e + 1) % poly.size()], uv = v - u;
            const double t = std::max(0.0, std::min(1.0, std::real((q - u) * std::conj(uv)) / std::max(std::norm(uv), 1e-300)));
            dm = std::min(dm, std::abs(q - (u + t * uv)));
          }
          return dm;
        };
        const bool vside = inside(c);
        int first0 = P;
        double margin = dist(c);
        C w = c;
        for (int k = 1; k < P; k++) {
          w = w * w + c;
          margin = std::min(margin, dist(w));
          if (inside(w) != vside) { first0 = k + 1; break; }   // ν index: f^k(c) is the orbit point at time k + 1
        }
        snprintf(line, sizeof(line), "%s %.17g %.17g %d %d %.3g", J.name.c_str(), c.real(), c.imag(), P, first0, margin);
        out[i] = line;
      }
    });
  for (auto& t : pool) t.join();
  for (const auto& l : out) printf("%s\n", l.c_str());
}
