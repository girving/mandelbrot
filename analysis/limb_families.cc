// Seahorse-valley families: non-renormalizable primitive components of the cardioid limbs k/(2k+1), up to k.
//
//   ./build/release/limb_families stats J                      # families with extra period j ≤ J, k = 2, 3, 4
//   ./build/release/limb_families jobs J k1,k2,... [threads]   # bulb_batch P = 0 jobs for each family at each k
//   ./build/release/limb_families keys J                       # j, index, key words per family
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
#include <algorithm>
#include <atomic>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
#include <string>
#include <thread>
#include <vector>
using namespace mandelbrot;
using std::string;
using std::vector;
typedef std::complex<double> C;

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
    // Per family, rays at the first three k (increasing); later k by Newton from a quadratic extrapolation of
    // 1/δ in k, accepted if it lands within 5% of the last spacing from the prediction (else rays again).  Output
    // "name_k c_re c_im p |s|^2 |Λ| |β|" with s = 1/(β Λ²), Λ = Π_{i<p} 2 z_i, β = Σ_{i<p} 1/Π_{j≤i} 2 z_j.
    threads = atoi(argv[2]);
    struct Fam { string name, lo, hi; vector<int> ks; };
    vector<Fam> fams;
    char name[256], lo[4096], hi[4096], ks[65536];
    while (scanf("%255s %4095s %4095s %65535s", name, lo, hi, ks) == 4) {
      auto k = parse_ks(ks);
      std::sort(k.begin(), k.end());
      fams.push_back({name, lo, hi, k});
    }
    vector<string> out(fams.size());
    std::atomic<int64_t> next(0), rays(0), failed(0);
    vector<std::thread> pool;
    for (int t = 0; t < threads; t++)
      pool.emplace_back([&]() {
        for (int64_t n; (n = next.fetch_add(1)) < int64_t(fams.size());) {
          const auto& f = fams[n];
          vector<C> cs;
          vector<int> done;
          string text;
          for (const int k : f.ks) {
            string pre;
            for (int i = 1; i < k; i++) pre += "01";
            const string lo = pre + f.lo, hi = pre + f.hi;
            const int p = lo.size();
            C c;
            bool ok = false;
            const int m = cs.size();
            if (m >= 3) {
              // Lagrange extrapolation of u = 1/δ through the last three (k, u)
              const double k0 = done[m-3], k1 = done[m-2], k2 = done[m-1];
              const C u0 = 1.0 / (cs[m-3] + 0.75), u1 = 1.0 / (cs[m-2] + 0.75), u2 = 1.0 / (cs[m-1] + 0.75);
              const double x = k;
              const C u = u0 * ((x - k1) * (x - k2) / ((k0 - k1) * (k0 - k2))) +
                          u1 * ((x - k0) * (x - k2) / ((k1 - k0) * (k1 - k2))) +
                          u2 * ((x - k0) * (x - k1) / ((k2 - k0) * (k2 - k1)));
              const C pred = 1.0 / u - 0.75;
              c = pred;
              ok = center(c, p) && std::abs(c - pred) < 0.05 * std::abs(cs[m-1] - cs[m-2]);
            }
            if (!ok) {
              rays++;
              C a = ray_in(lo, 2 * p + 4), b = ray_in(hi, 2 * p + 4);
              const bool oka = center(a, p), okb = center(b, p);
              ok = oka && okb && std::abs(a - b) < 1e-9 * std::abs(a + 0.75);
              c = a;
            }
            if (!ok) { failed++; fprintf(stderr, "%s_%d: failed\n", f.name.c_str(), k); continue; }
            C z = 0, prod = 1, beta = 0;
            for (int i = 1; i < p; i++) { z = z * z + c; prod *= 2.0 * z; beta += 1.0 / prod; }
            const double s2 = 1 / std::norm(beta * prod * prod);
            char line[512];
            snprintf(line, sizeof(line), "%s_%d %.17g %.17g %d %.17g %.17g %.17g\n", f.name.c_str(), k, c.real(),
                     c.imag(), p, s2, std::abs(prod), std::abs(beta));
            text += line;
            cs.push_back(c);
            done.push_back(k);
          }
          out[n] = text;
        }
      });
    for (auto& t : pool) t.join();
    for (const auto& s : out) fputs(s.c_str(), stdout);
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
