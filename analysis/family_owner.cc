// The family expansion, measured: M's escape-time tail split into per-root families.
//
// Near a p/q satellite root (period p = qP, child of a period-P parent) an escaping parameter's critical orbit
// passes through the root's gate: it lingers near a near-parabolic p-cycle, whose multiplier m satisfies
// m - 1 ≈ -q^2 (λ/λ0 - 1) with λ the parent's multiplier.  Each escaping parameter (step ≥ 64) is owned by the
// outermost gate it passes: the smallest Q ≤ Qscan such that the orbit returns, |z_i - z_{i-Q}| < eps, and Newton
// from its closest return finds a Q-cycle with |m - 1| < μ0; the owner is the catalog root of period Q nearest to c.
// Ownership is nested: deeper roots inside a root's gate disk (seahorse-valley copies, whose escapes take many
// passes through the gate) belong to its family, so the families partition the exterior.  The family expansion
// says root r's owned tail is A_r G_q(n / P_r), one profile per type q, with A_r = |dc/dλ_parent|^2 at the root.
//
// Roots come from component centers (hyperbolic's dump: "p c_re c_im area" per line): continue (z, c) along the
// component's own multiplier from 0 to 1 - 1e-6; a period P | p cycle nearby with multiplier a primitive (p/P)-th
// root of unity makes it a satellite, whose exact root and A_r come from Newton in the parent's coordinate.
// Primitive roots (cusps) get A_r = |κ_r / κ_cardioid|^2, with c ≈ c_r + κ_r (λ - 1)^2.
//
// Parameters are sampled uniformly in [-2, 0.5] × [0, 1.2] and classified like escape_tree's leaves (Orbit<double>,
// with Newton certificates).  Areas are doubled for the lower half plane: the total is right, but an off-axis root's
// row also counts its conjugate's family.
//   ./build/release/family_owner components.txt max_period samples_log2 max_iter_log2 [seed] [eps] [Qscan] [mu0]
// (defaults eps 0.1, Qscan 64, μ0 0.5).  CPU threads: $MANDELBROT_THREADS; per-sample dump: $FAMILY_DUMP.
#include "engine.h"
#include "orbit.h"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <thread>
#include <vector>
using namespace mandelbrot;
typedef std::complex<double> Cx;

// Newton for (z, c) with f_c^P(z) = z, (f_c^P)'(z) = lam; also returns dc/dλ there
static bool solve(const int P, const Cx lam, Cx& z, Cx& c, Cx* dcdl = nullptr) {
  double last = INFINITY;
  for (int it = 0; it < 60; it++) {
    Cx x = z, xz = 1, xc = 0, xzz = 0, xzc = 0;
    for (int i = 0; i < P; i++) {
      const Cx nxzz = 2.0 * (xz * xz + x * xzz), nxzc = 2.0 * (xc * xz + x * xzc);
      xzz = nxzz; xzc = nxzc;
      xc = 2.0 * x * xc + 1.0;
      xz = 2.0 * x * xz;
      x = x * x + c;
    }
    const Cx F1 = x - z, F2 = xz - lam, a = xz - 1.0, b = xc, d = xzz, e = xzc, det = a * e - b * d;
    if (dcdl) *dcdl = a / det;
    const Cx dz = (F1 * e - b * F2) / det, dc = (a * F2 - d * F1) / det;
    z -= dz; c -= dc;
    last = std::abs(dz) + std::abs(dc);
    if (last < 1e-14 * (1 + std::abs(z))) return true;
  }
  return last < 1e-9;  // Near a root the system degenerates, and roundoff limits the final steps
}

struct RootInfo {
  int p, P, q;  // Period, parent period (p for primitive), satellite denominator (1 for primitive)
  Cx c;
  double A;     // Area scale in the multiplier coordinate
};

int main(int argc, char** argv) {
  if (argc < 5) { fprintf(stderr, "usage: family_owner components.txt max_period samples_log2 max_iter_log2 [seed] [eps] [Qscan] [mu0]\n"); return 1; }
  const int Pmax = atoi(argv[2]);
  const int64_t S = int64_t(1) << atoi(argv[3]), max_iter = int64_t(1) << atoi(argv[4]);
  const uint64_t seed = argc > 5 ? atoll(argv[5]) : 1;
  const double eps = argc > 6 ? atof(argv[6]) : 0.1;
  const int Qscan = argc > 7 ? atoi(argv[7]) : 64;
  const double mu0 = argc > 8 ? atof(argv[8]) : 0.5;
  if (Qscan < Pmax || Qscan > 255) { fprintf(stderr, "need max_period ≤ Qscan ≤ 255\n"); return 1; }

  // Roots of components up to period Pmax (upper half plane, including the real axis)
  std::vector<RootInfo> roots;
  roots.push_back({1, 1, 1, Cx(0.25, 0), 1.0});  // The main cardioid's cusp, κ = -1/4
  FILE* fp = fopen(argv[1], "r");
  if (!fp) { fprintf(stderr, "can't open %s\n", argv[1]); return 1; }
  int p, fails = 0;
  double cr, ci, area;
  while (fscanf(fp, "%d %lf %lf %lf", &p, &cr, &ci, &area) == 4) {
    if (p < 2 || p > Pmax || ci < -1e-12) continue;
    Cx z = 0, c(cr, ci);
    bool ok = true;
    const double end = 1 - 1e-6;
    for (int s = 1; s <= 200 && ok; s++) ok = solve(p, Cx(end * s / 200, 0), z, c);
    if (!ok) { fails++; continue; }
    // Satellite if a period P | p cycle near the root's cycle has multiplier a primitive (p/P)-th root of unity
    int P = p;
    Cx zP = z, mP = 1;
    for (int d = 1; d < p && P == p; d++) {
      if (p % d) continue;
      Cx w = z;
      for (int it = 0; it < 50; it++) {
        Cx x = w, m = 1;
        for (int i = 0; i < d; i++) { m = 2.0 * x * m; x = x * x + c; }
        w -= (x - w) / (m - 1.0);
      }
      Cx x = w, m = 1;
      for (int i = 0; i < d; i++) { m = 2.0 * x * m; x = x * x + c; }
      if (std::abs(x - w) < 1e-10 && std::abs(w - z) < 0.1 && std::abs(std::pow(m, p / d) - 1.0) < 1e-3) {
        P = d; zP = w; mP = m;
      }
    }
    RootInfo r{p, P, p / P, c, 0};
    if (P < p) {
      // |dc/dλ|^2 in the parent's multiplier coordinate, at the parent's cycle
      const Cx lam0 = std::polar(1.0, std::round(std::arg(mP) * (p / P) / (2 * M_PI)) * 2 * M_PI / (p / P));
      Cx dcdl;
      if (!solve(P, lam0, zP, c, &dcdl)) { fails++; continue; }
      r.c = c;
      r.A = std::norm(dcdl);
    } else {
      // Cusp: c ≈ c_r + κ (λ - 1)^2, from the continuation at λ = 1 - h; relative to the cardioid's κ = -1/4
      const double h = 1e-3;
      if (!solve(p, 1.0, z, c)) { fails++; continue; }
      r.c = c;
      Cx z2 = z, c2 = c;
      if (!solve(p, Cx(1 - h, 0), z2, c2)) { fails++; continue; }
      const Cx kappa = (c2 - c) / (h * h);
      r.A = std::norm(kappa / -0.25);
    }
    roots.push_back(r);
  }
  fclose(fp);
  // Roots by period, for nearest-root lookup
  std::vector<std::vector<int>> by_period(Pmax + 1);
  for (int i = 0; i < int(roots.size()); i++) by_period[roots[i].p].push_back(i);
  printf("%zu roots up to period %d (%d failed)\n", roots.size(), Pmax, fails);

  // Sampling: per thread histograms [root or deep (index = roots.size()) or unattributed (+1)][octave]
  const int R = int(roots.size()), NJ = 40;
  const int T = cpu_threads();
  std::vector<std::vector<int64_t>> hist(T, std::vector<int64_t>((R + 2) * NJ, 0));
  // $FAMILY_DUMP: write "owner gate_period steps c_re c_im" for every escape at step ≥ 128
  const char* dump_path = getenv("FAMILY_DUMP");
  FILE* dump = dump_path ? fopen(dump_path, "w") : nullptr;
  std::mutex dump_mutex;
  std::vector<std::thread> pool;
  for (int t = 0; t < T; t++)
    pool.emplace_back([&, t]() {
      auto& h = hist[t];
      for (int64_t s = t; s < S; s += T) {
        const double x = -2 + 2.5 * uniform(seed, s, 0), y = 1.2 * uniform(seed, s, 1);
        Orbit<double> o;
        if (o.start(x, y, 64)) continue;
        o.finish(max_iter, 256);
        const Escape e = o.result();
        if (e.steps < 64) continue;
        const int j = std::min(NJ - 1, 63 - __builtin_clzll(uint64_t(e.steps)));
        // The owner is the outermost gate the orbit passes: the smallest Q ≤ Qscan such that the orbit returns,
        // |z_i - z_{i-Q}| < eps, and Newton from its closest return finds a near-parabolic Q-cycle, |m - 1| < μ0.
        // Near a p/q satellite root the child's multiplier is m - 1 ≈ -q^2 (λ/λ0 - 1), so this is a disk around the
        // root in its parent's multiplier coordinate (radius ≈ μ0/q^2), and deeper roots inside it (seahorse-valley
        // copies, whose escapes take many passes through the gate) belong to its family.  No gate ≤ Qscan: deep.
        const Cx c(x, y);
        Cx ring[256], zbest[256];
        double dbest[256];
        for (int Q = 0; Q <= Qscan; Q++) dbest[Q] = INFINITY;
        Cx z = 0;
        for (int64_t i = 0; i < e.steps; i++) {
          if (i >= Qscan)
            for (int Q = 1; Q <= Qscan; Q++) {
              const double d = std::norm(z - ring[(i - Q) & 255]);
              if (d < dbest[Q]) { dbest[Q] = d; zbest[Q] = z; }
            }
          ring[i & 255] = z;
          z = z * z + c;
        }
        int Qs = 0;
        for (int Q = 1; Q <= Qscan && !Qs; Q++) {
          if (!(dbest[Q] < eps * eps)) continue;
          Cx w = zbest[Q], m = 1;
          bool ok = false;
          for (int it = 0; it < 60; it++) {
            Cx x = w;
            m = 1;
            for (int k = 0; k < Q; k++) { m = 2.0 * x * m; x = x * x + c; }
            const Cx step = (x - w) / (m - 1.0);
            w -= step;
            if (!(std::norm(step) < 1)) break;
            if (std::norm(step) < 1e-18) { ok = true; break; }  // A near-double root: ~sqrt(ε) accuracy
          }
          if (ok && std::norm(w - zbest[Q]) < eps * eps && std::norm(m - 1.0) < mu0 * mu0) Qs = Q;
        }
        int owner = R + 1;  // No gate ≤ Qscan
        if (Qs) {
          owner = R;  // A gate of period in (Pmax, Qscan], or no listed root
          double dmin = INFINITY;
          if (Qs <= Pmax)
            for (const int i : by_period[Qs]) {
              const double d = std::norm(c - roots[i].c);
              if (d < dmin) { dmin = d; owner = i; }
            }
        }
        h[owner * NJ + j]++;
        if (dump && e.steps >= 128) {
          const std::lock_guard<std::mutex> lock(dump_mutex);
          fprintf(dump, "%d %d %lld %.17g %.17g\n", owner, Qs, (long long)e.steps, x, y);
        }
      }
    });
  for (auto& th : pool) th.join();
  if (dump) fclose(dump);
  std::vector<int64_t> H((R + 2) * NJ, 0);
  for (auto& h : hist) for (size_t i = 0; i < H.size(); i++) H[i] += h[i];
  const double w = 2 * 2.5 * 1.2 / double(S);  // Area per sample, doubled
  printf("area per sample %.6e; columns: octaves j = 6 … %d (area escaping in [2^j, 2^{j+1}))\n", w, NJ - 1);
  printf("total     ");
  for (int j = 6; j < NJ; j++) { int64_t s = 0; for (int r = 0; r < R + 2; r++) s += H[r * NJ + j]; printf(" %.4e", s * w); }
  printf("\n");
  for (int r = 0; r < R + 2; r++) {
    int64_t s = 0;
    for (int j = 0; j < NJ; j++) s += H[r * NJ + j];
    if (!s) continue;
    if (r < R) printf("root p %2d P %2d q %2d c %+.10f %+.10fi A %.6e:", roots[r].p, roots[r].P, roots[r].q,
                      roots[r].c.real(), roots[r].c.imag(), roots[r].A);
    else printf(r == R ? "deep (gate period > max_period, or unlisted):" : "no gate (period > Qscan):");
    for (int j = 6; j < NJ; j++) printf(" %lld", (long long)H[r * NJ + j]);
    printf("\n");
  }
}
