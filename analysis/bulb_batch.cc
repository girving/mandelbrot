// Batched bulb areas (bulb.h) on CPU threads or the GPU.
//
//   ./build/release/bulb_batch [--cuda] [--N 64] [--polish 2] < jobs > out
//   ./build/release/bulb_batch [--cuda] --all Q [--parent P c_re c_im] > out   (every child p/q with q ≤ Q, p ≤ q/2
//                                                                                for the cardioid; default parent)
// reads lines "key P c_re c_im p q" (key: any token naming the parent; P, (c_re, c_im): the parent's period and
// center; p/q the child) and prints "key p q center_re center_im area_hi F_hi w conv 0 area_lo F_lo" per line (the
// columns of bulb_areas with BULB_EXP=1, after the key), or "key p q failed <stage>".
#include "bulb.h"
#include "wall_time.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <string>
using namespace mandelbrot;

int main(int argc, char** argv) {
  BulbParams params;
  int all = 0, Pp = 1;
  double pcr = 0, pci = 0;
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--cuda")) params.cuda = true;
    else if (!strcmp(argv[i], "--N") && i + 1 < argc) params.N = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--polish") && i + 1 < argc) params.polish = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--substeps") && i + 1 < argc) params.substeps = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--all") && i + 1 < argc) all = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--parent") && i + 3 < argc) {
      Pp = atoi(argv[++i]); pcr = atof(argv[++i]); pci = atof(argv[++i]);
    }
    else { fprintf(stderr, "usage: bulb_batch [--cuda] [--N 64] [--polish 2] [--substeps 4] < jobs\n"); return 1; }
  }
  if (getenv("BULB_TOL")) params.accept = atof(getenv("BULB_TOL"));
  vector<std::string> keys;
  vector<BulbJob> jobs;
  char key[256];
  int P, p, q;
  double cr, ci;
  if (all)
    for (q = 2; q <= all; q++)
      for (p = 1; p < q; p++)
        if (std::gcd(p, q) == 1 && (Pp > 1 || 2 * p <= q)) {
          keys.push_back("all");
          jobs.push_back(bulb_job(Pp, Complex<double>(pcr, pci), p, q));
        }
  while (!all && scanf("%255s %d %lf %lf %d %d", key, &P, &cr, &ci, &p, &q) == 6) {
    keys.push_back(key);
    jobs.push_back(bulb_job(P, Complex<double>(cr, ci), p, q));
  }
  const auto t0 = wall_time();
  const auto res = bulb_areas(jobs, params);
  fprintf(stderr, "bulb_batch: %zu bulbs in %.3f s (%s)\n", jobs.size(), (wall_time() - t0).seconds(),
          params.cuda ? "gpu" : "cpu");
  static const char* why[] = {"ok", "parent", "center", "period", "area"};
  for (size_t i = 0; i < jobs.size(); i++) {
    const auto& r = res[i];
    if (r.status == bulb_ok)
      printf("%s %d %d %.17g %.17g %.17g %.17g %.17g %.1e 0 %.17g %.17g\n", keys[i].c_str(), jobs[i].p, jobs[i].q,
             r.center.r, r.center.i, r.area.x[0], r.F.x[0], r.w, r.conv, r.area.x[1], r.F.x[1]);
    else
      printf("%s %d %d failed %s\n", keys[i].c_str(), jobs[i].p, jobs[i].q, why[r.status]);
  }
}
