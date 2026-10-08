// How much of M's hyperbolic area lies in the satellite tree: the components reached from the main cardioid by
// satellite bifurcations alone (bulbs of bulbs of ...), as opposed to those inside primitive copies.
//
// Reads hyperbolic's dump ("p c_re c_im area" per component).  For each component: continue its own multiplier from
// the center towards 1 to find the root; a period P | p cycle there whose multiplier is a primitive (p/P)-th root
// of unity makes it a satellite, and continuing that cycle's multiplier back to 0 lands on the parent's center,
// matched against the catalog.  A component is in the satellite tree if it is the cardioid or a satellite of a
// component in the tree; otherwise its primitive ancestor (the first primitive component up its chain of parents)
// has period > 1, and it lies in that copy.
//   ./build/release/satellite_tree components.txt
// CPU threads: $MANDELBROT_THREADS.
#include <atomic>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <thread>
#include <vector>
typedef std::complex<double> Cx;

// Newton for (z, c) with f_c^n(z) = z, (f_c^n)'(z) = mu
static bool solve(const int n, const Cx mu, Cx& z, Cx& c) {
  double last = INFINITY;
  for (int it = 0; it < 40; it++) {
    Cx x = z, xz = 1, xc = 0, xzz = 0, xzc = 0;
    for (int i = 0; i < n; i++) {
      const Cx nxzz = 2.0 * (xz * xz + x * xzz), nxzc = 2.0 * (xc * xz + x * xzc);
      xzz = nxzz; xzc = nxzc;
      xc = 2.0 * x * xc + 1.0;
      xz = 2.0 * x * xz;
      x = x * x + c;
    }
    const Cx F1 = x - z, F2 = xz - mu, a = xz - 1.0, b = xc, d = xzz, e = xzc, det = a * e - b * d;
    const Cx dz = (F1 * e - b * F2) / det, dc = (a * F2 - d * F1) / det;
    z -= dz; c -= dc;
    last = std::abs(dz) + std::abs(dc);
    if (!(last < 1)) return false;
    if (last < 1e-15 * (1 + std::abs(c))) return true;
  }
  return last < 1e-8;  // Near a root the system degenerates
}

struct Comp {
  int p;
  Cx c;
  double area;
  int parent = -1;     // Satellite parent's index, or -1 for primitive (or unresolved)
  bool resolved = true;
};

int main(int argc, char** argv) {
  if (argc < 2) { fprintf(stderr, "usage: satellite_tree components.txt\n"); return 1; }
  std::vector<Comp> comps;
  FILE* fp = fopen(argv[1], "r");
  if (!fp) { fprintf(stderr, "can't open %s\n", argv[1]); return 1; }
  int p;
  double cr, ci, a;
  while (fscanf(fp, "%d %lf %lf %lf", &p, &cr, &ci, &a) == 4) comps.push_back({p, Cx(cr, ci), a});
  fclose(fp);
  int pmax = 0;
  for (const auto& w : comps) pmax = std::max(pmax, w.p);
  std::vector<std::vector<int>> by_period(pmax + 1);
  for (int i = 0; i < int(comps.size()); i++) by_period[comps[i].p].push_back(i);

  std::atomic<size_t> next(0);
  const char* te = getenv("MANDELBROT_THREADS");
  const int T = te ? atoi(te) : int(std::thread::hardware_concurrency());
  std::vector<std::thread> pool;
  for (int t = 0; t < T; t++)
    pool.emplace_back([&]() {
      for (size_t i; (i = next++) < comps.size();) {
        Comp& w = comps[i];
        if (w.p == 1) continue;
        // Root: own multiplier to 1 - 1e-6
        Cx z = 0, c = w.c;
        bool ok = true;
        for (int s = 1; s <= 200 && ok; s++) ok = solve(w.p, Cx((1 - 1e-6) * s / 200, 0), z, c);
        if (!ok) { w.resolved = false; continue; }
        for (int d = 1; d < w.p && w.parent < 0; d++) {
          if (w.p % d) continue;
          // A period d cycle near z with multiplier a primitive (p/d)-th root of unity
          Cx u = z;
          for (int it = 0; it < 50; it++) {
            Cx x = u, m = 1;
            for (int k = 0; k < d; k++) { m = 2.0 * x * m; x = x * x + c; }
            u -= (x - u) / (m - 1.0);
          }
          Cx x = u, m = 1;
          for (int k = 0; k < d; k++) { m = 2.0 * x * m; x = x * x + c; }
          if (!(std::abs(x - u) < 1e-9 && std::abs(u - z) < 1 && std::abs(std::pow(m, w.p / d) - 1.0) < 1e-3))
            continue;
          // Primitive (p/d)-th root of unity: m^k ≠ 1 for proper divisors k of p/d
          const int qq = w.p / d;
          bool primitive_root = true;
          for (int k = 1; k < qq; k++)
            if (qq % k == 0 && std::abs(std::pow(m, k) - 1.0) < 1e-3) primitive_root = false;
          if (!primitive_root) continue;
          // Continue the parent's multiplier from m back to 0: lands on the parent's center
          Cx zp = u, cp = c;
          bool okp = true;
          for (int s = 199; s >= 0 && okp; s--) okp = solve(d, m * (s / 200.0), zp, cp);
          if (!okp) { w.resolved = false; break; }
          double best = INFINITY;
          int bi = -1;
          for (const int j : by_period[d]) {
            const double dist = std::abs(comps[j].c - cp);
            if (dist < best) { best = dist; bi = j; }
          }
          if (best < 1e-8) w.parent = bi;
          else w.resolved = false;
        }
      }
    });
  for (auto& th : pool) th.join();

  // Primitive ancestor of each component: follow satellite parents (periods decrease, so recurse in order)
  std::vector<int> anc(comps.size(), -1);
  for (int per = 1; per <= pmax; per++)
    for (const int i : by_period[per]) anc[i] = comps[i].parent < 0 ? i : anc[comps[i].parent];
  // Per period: total area, satellite-tree area, area in copies by primitive ancestor period, unresolved
  printf("period  components  total area       satellite tree   in copies        unresolved  (copy area by ancestor period)\n");
  double T_all = 0, T_sat = 0, T_copy = 0;
  for (int per = 1; per <= pmax; per++) {
    double all = 0, sat = 0, copy = 0;
    int unres = 0;
    std::map<int, double> by_anc;
    for (const int i : by_period[per]) {
      all += comps[i].area;
      if (!comps[i].resolved) { unres++; continue; }
      const int r = anc[i];
      if (comps[r].p == 1) sat += comps[i].area;
      else { copy += comps[i].area; by_anc[comps[r].p] += comps[i].area; }
    }
    T_all += all; T_sat += sat; T_copy += copy;
    printf("  %2d    %6zu     %.10e  %.10e  %.10e  %d   ", per, by_period[per].size(), all, sat, copy, unres);
    for (const auto& kv : by_anc) printf(" %d:%.2e", kv.first, kv.second);
    printf("\n");
  }
  printf("total   %.12f  satellite tree %.12f (%.4f%%)  copies %.12f (%.4f%%)\n", T_all, T_sat, 100 * T_sat / T_all,
         T_copy, 100 * T_copy / T_all);
}
