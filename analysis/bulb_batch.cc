// Batched bulb areas (bulb.h) on CPU threads or the GPU.
//
//   ./build/release/bulb_batch [--cuda] [--N 64] [--polish 2] [--out file] < jobs
//   ./build/release/bulb_batch [--cuda] --all Q [--parent P c_re c_im] > out   (every child p/q with q ≤ Q, p ≤ q/2
//                                                                                for the cardioid; default parent)
// reads lines "key P c_re c_im p q" (key: any token naming the parent; P, (c_re, c_im): the parent's period and
// center; p/q the child) and prints "key p q center_re center_im area_hi F_hi w conv 0 area_lo F_lo" per line (the
// columns of bulb_areas with BULB_EXP=1, after the key), or "key p q failed <stage>".  With --local, reads lines
// "key p c_re_hi c_re_lo c_im_hi c_im_lo" instead: deep components of period p whose centers are known to double-double
// precision (bulb_job_local), printed as "key 0 p ..." like P = 0 jobs.
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
  bool local = false;
  const char* out_path = nullptr;
  double pcr = 0, pci = 0;
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--cuda")) params.cuda = true;
    else if (!strcmp(argv[i], "--N") && i + 1 < argc) params.N = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--polish") && i + 1 < argc) params.polish = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--substeps") && i + 1 < argc) params.substeps = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--all") && i + 1 < argc) all = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--out") && i + 1 < argc) out_path = argv[++i];
    else if (!strcmp(argv[i], "--local")) local = true;
    else if (!strcmp(argv[i], "--parent") && i + 3 < argc) {
      Pp = atoi(argv[++i]); pcr = atof(argv[++i]); pci = atof(argv[++i]);
    }
    else { fprintf(stderr, "usage: bulb_batch [--cuda] [--N 64] [--polish 2] [--substeps 4] [--out file] < jobs\n"); return 1; }
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
  double crl, cil;
  while (!all && local && scanf("%255s %d %lf %lf %lf %lf", key, &q, &cr, &crl, &ci, &cil) == 6) {
    keys.push_back(key);
    jobs.push_back(bulb_job_local(Complex<double>(cr, ci), Complex<double>(crl, cil), q));
  }
  while (!all && !local && scanf("%255s %d %lf %lf %d %d", key, &P, &cr, &ci, &p, &q) == 6) {
    keys.push_back(key);
    jobs.push_back(bulb_job(P, Complex<double>(cr, ci), p, q));
  }
  const auto t0 = wall_time();
  const auto res = bulb_areas(jobs, params);
  fprintf(stderr, "bulb_batch: %zu bulbs in %.3f s (%s)\n", jobs.size(), (wall_time() - t0).seconds(),
          params.cuda ? "gpu" : "cpu");
  static const char* why[] = {"ok", "parent", "center", "period", "area"};
  FILE* out = out_path ? fopen(out_path, "w") : stdout;
  if (!out) { fprintf(stderr, "can't open %s\n", out_path); return 1; }
  int failed = 0;
  for (size_t i = 0; i < jobs.size(); i++) {
    const auto& r = res[i];
    if (r.status == bulb_ok)
    {
      // area and F parts with 60 digits: %.17g identifies a double but is up to ~5e-18 relative off its value, which
      // the lower parts cannot repair (the 3e-17 floor of the old M-side family fits; %.25g still leaves 5e-26)
      fprintf(out, "%s %d %d %.17g %.17g %.60g %.60g %.17g %.1e 0 %.60g %.60g", keys[i].c_str(), jobs[i].p, jobs[i].q,
             r.center.r, r.center.i, r.area.x[0], r.F.x[0], r.w, r.conv, r.area.x[1], r.F.x[1]);
      for (int k = 2; k < BULB_E; k++) fprintf(out, " %.60g %.60g", r.area.x[k], r.F.x[k]);   // bulb_batch3's third parts
      fprintf(out, "\n");
    }
    else {
      fprintf(out, "%s %d %d failed %s\n", keys[i].c_str(), jobs[i].p, jobs[i].q, why[r.status]);
      failed++;
    }
  }
  if (out != stdout) fclose(out);
  fprintf(stderr, "bulb_batch: %d failed\n", failed);
}
