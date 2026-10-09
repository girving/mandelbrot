// Seahorse-valley families: non-renormalizable primitive components of the cardioid limbs k/(2k+1), up to k.
//
//   ./build/release/limb_families stats J                      # families with extra period j ≤ J, k = 2, 3, 4
//   ./build/release/limb_families jobs J k1,k2,... [threads]   # bulb_batch P = 0 jobs for each family at each k
//   ./build/release/limb_families keys J                       # j, index, key words per family
//   ./build/release/limb_families limb J k                     # every NRP of limb k: j and key words
//   ./build/release/limb_families custom threads < lines "name lo hi k1,k2,..."   # jobs for given keys
//   ./build/release/limb_families size threads < lines "name lo hi k1,k2,..."     # centers and size estimates only
//
// Every angle in the limb k/(2k+1) starts with (01)^{k-1}0 (the common prefix of its wake words (01)^{k-1}001 and
// (01)^{k-1}010), and stripping (01)^{k-1} from both angle words of a component gives a key independent of k: the
// family.  `stats` checks that with lavaurs restricted to the wake, for the non-renormalizable primitives
// (maximal_tuning finds no copy containing them; the limb's own bulb tunes its satellite tree, which is not a
// family).  A family's extra period is j = p - (2k + 1).  `jobs` takes the families from limb K0 = 4 (periods ≤ 9 + J; limb 2 is preasymptotic for j ≥ 5),
// builds each member's angle words at the requested k, traces both parameter rays to near the root, and runs
// Newton for the center from each endpoint; the two must agree.  Output lines "j<j>_<i>_<k> 0 c_re c_im 0 p" on
// stdout, failures on stderr.
#include "angles.h"
#include "complex.h"
#include "expansion.h"
#include "expansion_arith.h"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <mutex>
#include <set>
#include <string>
#include <thread>
#include <vector>
using namespace mandelbrot;
using std::string;
using std::vector;
typedef std::complex<double> C;
typedef Expansion<2> E2;
typedef Complex<E2> Ce;

static string bits(const uint64_t k, const int q) {
  string s(q, '0');
  for (int i = 0; i < q; i++) s[q - 1 - i] = '0' + (k >> i & 1);
  return s;
}

struct Family { string lo, hi; };  // Angle words with (01)^{k-1} stripped

static vector<Family> families(const int k, const int J) {
  const int q = 2 * k + 1;
  const auto roots = lavaurs(q + J, cardioid_wake(k, q));
  const auto parent = maximal_tuning(roots);
  const string pre = [&] { string s; for (int i = 1; i < k; i++) s += "01"; return s; }();
  vector<Family> fs;
  for (size_t i = 0; i < roots.size(); i++) {
    const auto& w = roots[i].w;
    if (roots[i].satellite || parent[i] >= 0 || w.lo.q == q) continue;
    const string lo = bits(w.lo.k, w.lo.q), hi = bits(w.hi.k, w.hi.q);
    if (lo.compare(0, pre.size(), pre) || hi.compare(0, pre.size(), pre)) {
      fprintf(stderr, "limb %d: angles %s %s lack the prefix\n", k, lo.c_str(), hi.c_str());
      exit(1);
    }
    fs.push_back({lo.substr(pre.size()), hi.substr(pre.size())});
  }
  return fs;
}

// Parameter ray of the periodic angle with word w, from potential log(er) inward through `levels` levels
static C ray_in(const string& w, const int levels, const int S = 8) {
  const double er = 65536;
  const int p = w.size();
  const auto angle = [&](const int shift) {  // 2^shift θ mod 1 from its first 60 bits
    double t = 0, f = 0.5;
    for (int i = 0; i < 60; i++, f /= 2) t += f * (w[(shift + i) % p] - '0');
    return t;
  };
  C c = std::polar(er, 2 * M_PI * angle(0));
  for (int n = 1; n <= levels; n++) {
    const double t = angle((n - 1) % p);
    for (int j = 0; j < S; j++) {
      const C target = std::polar(std::pow(er, std::pow(0.5, (j + 1.0) / S)), 2 * M_PI * t);
      for (int it = 0; it < 64; it++) {
        C z = 0, dz = 0;
        for (int m = 0; m < n; m++) { dz = 2.0 * z * dz + 1.0; z = z * z + c; }
        const C step = (z - target) / dz;
        c -= step;
        if (std::abs(step) < 1e-15 * std::abs(c)) break;
      }
    }
  }
  return c;
}

// Newton for the center from c; converged if the last step is below 1e-11 |δ| (roundoff stops it near 1e-14 at
// periods ~100)
static bool center(C& c, const int p) {
  double last = INFINITY;
  for (int it = 0; it < 100; it++) {
    C z = 0, dz = 0;
    for (int i = 0; i < p; i++) { dz = 2.0 * z * dz + 1.0; z = z * z + c; }
    const C step = z / dz;
    c -= step;
    last = std::abs(step);
    if (last < 1e-15 * std::abs(c + 0.75)) break;
  }
  return last < 1e-11 * std::abs(c + 0.75);
}

// Center in double-double, from a guess: (optionally) double Newton as far as it goes, then simplified Newton with the
// residual f^p(0) in Complex<E2> and the Jacobian in double.  Then the size estimate from the E2 orbit rounded
// pointwise (each factor 2 z_i to 1e-16 relative).  Converged if the last step is below 1e-10 |s|.
struct Center { Ce c; double s2, lam, beta; };
static bool center_e2(Ce c, const int p, Center& out, const bool double_first) {
  if (double_first) {
    C g(double(c.r), double(c.i));
    center(g, p);
    c = Ce(E2(g.real()), E2(g.imag()));
  }
  double step = INFINITY;
  const auto orbit = [&](auto&& visit) {  // visit(z_i) for i = 1..p, z_i in E2
    Ce z(E2(0.0), E2(0.0));
    for (int i = 1; i <= p; i++) { z = sqr(z) + c; visit(i, z); }
  };
  for (int it = 0; it < 12; it++) {
    C dz = 0, zd = 0, zp = 0;
    orbit([&](const int i, const Ce& z) {  // dz_i = 2 z_{i-1} dz_{i-1} + 1 with z_0 = 0
      dz = 2.0 * zd * dz + 1.0;
      zd = C(double(z.r), double(z.i));
      if (i == p) zp = zd;
    });
    const C d = zp / dz;
    c = c - Ce(E2(d.real()), E2(d.imag()));
    const double last = step;
    step = std::abs(d);
    if (step < 1e-31 || step > 0.5 * last) break;
  }
  C prod = 1, beta = 0;
  orbit([&](const int i, const Ce& z) {
    if (i < p) { prod *= 2.0 * C(double(z.r), double(z.i)); beta += 1.0 / prod; }
  });
  out = {c, 1 / std::norm(beta * prod * prod), std::abs(prod), std::abs(beta)};
  return step < 1e-10 * std::sqrt(out.s2);
}

int main(int argc, char** argv) {
  if (argc < 3) { fprintf(stderr, "usage: limb_families stats J | jobs J k1,k2,... [threads]\n"); return 1; }
  const string mode = argv[1];
  const int J = atoi(argv[2]);
  if (mode == "stats") {
    std::map<int, std::set<std::pair<string, string>>> keys;
    for (const int k : {2, 3, 4, 5}) {
      vector<int> count(J + 1);
      for (const auto& f : families(k, J)) {
        const int j = int(f.lo.size()) + 2 * (k - 1) - (2 * k + 1);
        count[j]++;
        keys[k].insert({f.lo, f.hi});
      }
      printf("limb %d/%d: families per extra period j = 0..%d:", k, 2 * k + 1, J);
      for (int j = 0; j <= J; j++) printf(" %d", count[j]);
      printf("\n");
    }
    for (const auto [k0, k] : {std::pair(2, 3), std::pair(3, 4), std::pair(4, 5), std::pair(3, 5)}) {
      int same = 0, missing = 0;
      for (const auto& f : keys[k0]) (keys[k].count(f) ? same : missing)++;
      int extra = 0;
      for (const auto& f : keys[k]) extra += !keys[k0].count(f);
      printf("limb %d vs %d: %d keys shared, %d missing, %d new\n", k0, k, same, missing, extra);
      vector<int> miss(J + 1), shown(J + 1);
      for (const auto& f : keys[k0])
        if (!keys[k].count(f)) {
          const int j = f.first.size() - 3;  // key length = p - 2(k-1) = j + 3
          miss[j]++;
          if (getenv("SHOW") && shown[j]++ < 3) printf("  missing j=%d: %s %s\n", j, f.first.c_str(), f.second.c_str());
        }
      for (const auto& f : keys[k])
        if (!keys[k0].count(f) && getenv("SHOW") && shown[f.first.size() - 3]++ < 6)
          printf("  new j=%zu: %s %s\n", f.first.size() - 3, f.first.c_str(), f.second.c_str());
      printf("  missing per j:");
      for (int j = 0; j <= J; j++) printf(" %d", miss[j]);
      printf("\n");
    }
    return 0;
  }
  if (mode == "limb") {  // every NRP of limb k with extra period ≤ J: "j lo hi" (key words, (01)^{k-1} stripped)
    const int k = atoi(argv[3]);
    for (const auto& f : families(k, J))
      printf("%d %s %s\n", int(f.lo.size()) + 2 * (k - 1) - (2 * k + 1), f.lo.c_str(), f.hi.c_str());
    return 0;
  }
  if (mode == "keys") {  // j, index, key words (limb K0's families, as `jobs` numbers them)
    const int K0 = getenv("K0") ? atoi(getenv("K0")) : 4;
    vector<int> index(J + 1);
    for (const auto& f : families(K0, J)) {
      const int j = f.lo.size() - 3;
      printf("%d %d %s %s\n", j, index[j]++, f.lo.c_str(), f.hi.c_str());
    }
    return 0;
  }
  // Jobs: a family (name and key words) at one k
  struct Job { string name, lo, hi; int k; };
  vector<Job> jobs;
  const auto parse_ks = [](const char* s) {
    vector<int> ks;
    while (*s) { ks.push_back(strtol(s, (char**)&s, 10)); if (*s == ',') s++; }
    return ks;
  };
  int threads = 2, nfam = 0;
  if (mode == "size") {  // stdin lines "name lo hi k1,k2,...": centers and size estimates only (no areas)
    // Per family, k runs through every integer from its least to its largest requested k (only requested k are
    // printed).  The first 7 by rays (both rays, double-double centers that must agree); later ones by Newton in
    // double-double from the 7-point extrapolation of kδ (≈ iπ/2, slowly varying; integer weights, exact in E2),
    // accepted only if Newton moves less than 1% of the component's size |s| (a different root of f^p(0) is at
    // least ~|s| away), else rays again.  Output, as each family finishes: "name_k δre_hi δre_lo δim_hi δim_lo p
    // |s|^2 |Λ| |β|" with δ = c + 3/4, s = 1/(β Λ²), Λ = Π_{i<p} 2 z_i, β = Σ_{i<p} 1/Π_{j≤i} 2 z_j.
    threads = atoi(argv[2]);
    struct Fam { string name, lo, hi; vector<int> ks; };
    vector<Fam> fams;
    char name[256], lo[4096], hi[4096], ks[65536];
    while (scanf("%255s %4095s %4095s %65535s", name, lo, hi, ks) == 4) {
      auto k = parse_ks(ks);
      std::sort(k.begin(), k.end());
      fams.push_back({name, lo, hi, k});
    }
    std::atomic<int64_t> next(0), rays(0), failed(0);
    std::mutex io;
    vector<std::thread> pool;
    const auto delta = [](const Ce& c) { return Ce(c.r + E2(0.75), c.i); };
    const auto to_c = [](const Ce& z) { return C(double(z.r), double(z.i)); };
    for (int t = 0; t < threads; t++)
      pool.emplace_back([&]() {
        for (int64_t n; (n = next.fetch_add(1)) < int64_t(fams.size());) {
          const auto& f = fams[n];
          std::set<int> want(f.ks.begin(), f.ks.end());
          vector<Ce> kd;  // k δ for consecutive k
          string text;
          string pre;
          for (int i = 1; i < f.ks.front(); i++) pre += "01";
          for (int k = f.ks.front(); k <= f.ks.back(); k++, pre += "01") {
            const string lo = pre + f.lo, hi = pre + f.hi;
            const int p = lo.size();
            Center cen;
            bool ok = false;
            if (kd.size() >= 7) {
              // Extrapolate w = kδ to k from k-1..k-7: w = Σ_{i=1}^7 (-1)^{i+1} C(7,i) w_{k-i}
              static const int64_t binom[8] = {1, 7, 21, 35, 35, 21, 7, 1};
              const size_t m = kd.size();
              Ce w(E2(0.0), E2(0.0));
              for (int i = 1; i <= 7; i++) {
                const E2 b(i % 2 ? binom[i] : -binom[i]);
                w = w + Ce(b * kd[m - i].r, b * kd[m - i].i);
              }
              // δ = w / k in double-double: d0 = w/k in double, then the residual's quotient
              const auto div = [k](const E2 x) {
                const double d0 = double(x) / k;
                const E2 r = x - E2(int64_t(k)) * E2(d0);
                return E2(d0) + E2(double(r) / k);
              };
              const Ce pred(div(w.r) - E2(0.75), div(w.i));
              ok = center_e2(pred, p, cen, false) && std::abs(to_c(cen.c - pred)) < 1e-2 * std::sqrt(cen.s2);
            }
            if (!ok) {
              rays++;
              Center a, b;
              const auto seed = [](const C c) { return Ce(E2(c.real()), E2(c.imag())); };
              const bool oka = center_e2(seed(ray_in(lo, 2 * p + 4)), p, a, true);
              const bool okb = center_e2(seed(ray_in(hi, 2 * p + 4)), p, b, true);
              ok = oka && okb && std::abs(to_c(a.c - b.c)) < 1e-3 * std::sqrt(a.s2);
              cen = a;
            }
            if (!ok) {
              failed++;
              fprintf(stderr, "%s_%d: failed; stopping this family\n", f.name.c_str(), k);
              break;
            }
            const Ce d = delta(cen.c);
            kd.push_back(Ce(E2(int64_t(k)) * d.r, E2(int64_t(k)) * d.i));
            if (want.count(k)) {
              char line[512];
              snprintf(line, sizeof(line), "%s_%d %.17g %.17g %.17g %.17g %d %.17g %.17g %.17g\n", f.name.c_str(), k,
                       d.r.x[0], d.r.x[1], d.i.x[0], d.i.x[1], p, cen.s2, cen.lam, cen.beta);
              text += line;
            }
          }
          std::lock_guard<std::mutex> lock(io);
          fputs(text.c_str(), stdout);
          fflush(stdout);
        }
      });
    for (auto& t : pool) t.join();
    fprintf(stderr, "limb_families size: %zu families, %lld ray pairs, %lld failed\n", fams.size(), (long long)rays,
            (long long)failed);
    return 0;
  }
  if (mode == "custom") {  // stdin lines "name lo hi k1,k2,...": the key words need not come from a census
    threads = atoi(argv[2]);
    char name[256], lo[1024], hi[1024], ks[1024];
    while (scanf("%255s %1023s %1023s %1023s", name, lo, hi, ks) == 4) {
      nfam++;
      for (const int k : parse_ks(ks)) jobs.push_back({name, lo, hi, k});
    }
  } else {
    const auto ks = parse_ks(argv[3]);
    threads = argc > 4 ? atoi(argv[4]) : 2;
    const int K0 = getenv("K0") ? atoi(getenv("K0")) : 4;
    vector<int> index(J + 1);
    for (const auto& f : families(K0, J)) {
      const int j = f.lo.size() - 3;
      const string name = "j" + std::to_string(j) + "_" + std::to_string(index[j]++);
      nfam++;
      for (const int k : ks) jobs.push_back({name, f.lo, f.hi, k});
    }
  }
  vector<string> out(jobs.size());
  std::atomic<int64_t> next(0), failed(0);
  vector<std::thread> pool;
  for (int t = 0; t < threads; t++)
    pool.emplace_back([&]() {
      for (int64_t n; (n = next.fetch_add(1)) < int64_t(jobs.size());) {
        const auto& job = jobs[n];
        string pre;
        for (int i = 1; i < job.k; i++) pre += "01";
        const string lo = pre + job.lo, hi = pre + job.hi;
        const int p = lo.size();
        C a = ray_in(lo, 2 * p + 4), b = ray_in(hi, 2 * p + 4);
        const bool oka = center(a, p), okb = center(b, p), ok = oka && okb;
        char line[512];
        if (ok && std::abs(a - b) < 1e-9 * std::abs(a + 0.75)) {
          snprintf(line, sizeof(line), "%s_%d 0 %.17g %.17g 0 %d", job.name.c_str(), job.k, a.real(), a.imag(), p);
          out[n] = line;
        } else {
          fprintf(stderr, "%s_%d: rays disagree (%s %s): %.6g%+.6gi vs %.6g%+.6gi\n", job.name.c_str(), job.k,
                  lo.c_str(), hi.c_str(), a.real(), a.imag(), b.real(), b.imag());
          failed++;
        }
      }
    });
  for (auto& t : pool) t.join();
  for (const auto& s : out) if (s.size()) printf("%s\n", s.c_str());
  fprintf(stderr, "limb_families: %d families, %zu jobs, %lld failed\n", nfam, jobs.size(), (long long)failed);
}
