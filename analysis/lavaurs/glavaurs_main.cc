// The general p/q-root Lavaurs model (glavaurs.h): components, areas, and children.
//
//   glavaurs p q single nmax re0 re1 im0 im1 grid threads   # single-transit centers by grid Newton, with areas
//   glavaurs p q area threads < "name r n re im"              # areas: "name r n cre cim area_σ C conv cusp"
//   glavaurs p q tree cmin dmax [rloc]                        # single-transit centers as backward paths
//   glavaurs p q children threads < "name r n_u u_re u_im radius n_c c_re c_im jmin jmax"   # r-transit children
//   glavaurs p q locate < "name r re im"                      # Θ_r(σ), Θ_r', Π H'
//   glavaurs p q consist                                      # cross-petal branch conventions
// C = (4π² sin²(πp/q)/q⁴) area_σ is the family constant lim k⁴ area_M (limbs [CF(p/q), k]).
#include "glavaurs.h"
#include <atomic>
#include <cmath>
#include <complex>
#include <functional>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <string>
#include <thread>
#include <vector>
using namespace mandelbrot;
typedef Complex<double> Cd;
typedef std::complex<double> SCd;

int main(int argc, char** argv) {
  if (argc < 4) { fprintf(stderr, "usage: glavaurs p q mode ...\n"); return 1; }
  const int p = atoi(argv[1]), q = atoi(argv[2]);
  const std::string mode = argv[3];
  const GeneralLavaurs L(p, q, getenv("GL_SIDE") ? atoi(getenv("GL_SIDE")) : 1);
  const double K = 4 * M_PI * M_PI * std::pow(std::sin(M_PI * p / q), 2) / std::pow(q, 4);
  fprintf(stderr, "glavaurs %d/%d: a_{-q} = %.12g%+.12gi, β = %.12g%+.12gi, ζ0 = %.12g%+.12gi, K = %.10g\n", p, q,
          L.a[0].r, L.a[0].i, L.beta.r, L.beta.i, L.zeta0().r, L.zeta0().i, K);
  if (mode == "single") {
    const int nmax = atoi(argv[4]), grid = atoi(argv[9]), threads = atoi(argv[10]);
    const double r0 = atof(argv[5]), r1 = atof(argv[6]), i0 = atof(argv[7]), i1 = atof(argv[8]);
    struct Job { int n; Cd g; };
    std::vector<Job> jobs;
    for (int n = 0; n <= nmax; n++)
      for (int a = 0; a < grid; a++)
        for (int b = 0; b < grid; b++)
          jobs.push_back({n, Cd(r0 + (r1 - r0) * (a + 0.5) / grid, i0 + (i1 - i0) * (b + 0.5) / grid)});
    std::mutex mu;
    std::vector<std::pair<int, Cd>> found;
    std::atomic<int64_t> next(0);
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; t++)
      pool.emplace_back([&]() {
        for (int64_t i; (i = next++) < int64_t(jobs.size());) {
          Cd s = jobs[i].g;
          if (!L.center(1, jobs[i].n, s)) continue;
          if (!(s.r >= r0 - 1 && s.r <= r1 + 1 && s.i >= i0 - 1 && s.i <= i1 + 1)) continue;
          std::lock_guard<std::mutex> lock(mu);
          bool dup = false;
          for (const auto& f : found) dup |= f.first == jobs[i].n && std::hypot(f.second.r - s.r, f.second.i - s.i) < 1e-8;
          if (!dup) found.push_back({jobs[i].n, s});
        }
      });
    for (auto& t : pool) t.join();
    for (const auto& f : found) {
      Cd cen, a1;
      double ar, conv, cusp;
      if (L.area(1, f.first, f.second, cen, ar, conv, cusp, a1))
        printf("S 1 %d %.17g %.17g %.10e %.10e %.1e %.3e\n", f.first, cen.r, cen.i, ar, K * ar, conv, cusp);
      else
        printf("S 1 %d %.17g %.17g area-failed\n", f.first, f.second.r, f.second.i);
    }
    return 0;
  }
  if (mode == "area") {
    const int threads = argc > 4 ? atoi(argv[4]) : 2;
    struct Job { std::string name; int r, n; Cd g; };
    std::vector<Job> jobs;
    char name[256];
    int r, n;
    double sr, si;
    while (scanf("%255s %d %d %lf %lf", name, &r, &n, &sr, &si) == 5) jobs.push_back({name, r, n, Cd(sr, si)});
    std::vector<std::string> out(jobs.size());
    std::atomic<int64_t> next(0);
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; t++)
      pool.emplace_back([&]() {
        for (int64_t i; (i = next++) < int64_t(jobs.size());) {
          const auto& j = jobs[i];
          Cd cen, a1;
          double ar, conv, cusp;
          char line[512];
          if (L.area(j.r, j.n, j.g, cen, ar, conv, cusp, a1))
            snprintf(line, sizeof(line), "%s %d %d %.17g %.17g %.10e %.10e %.1e %.3e", j.name.c_str(), j.r, j.n, cen.r, cen.i, ar, K * ar, conv, cusp);
          else
            snprintf(line, sizeof(line), "%s %d %d failed", j.name.c_str(), j.r, j.n);
          out[i] = line;
        }
      });
    for (auto& t : pool) t.join();
    for (const auto& s : out) printf("%s\n", s.c_str());
    return 0;
  }
  if (mode == "tree") {
    // Single-transit centers as backward paths of the critical point into a repelling petal: y_0 = crit,
    // y_{i+1} ∈ f^{-1}(y_i) until |y| < r_loc in a repelling petal (inside r_switch the local inverse branch is followed
    // deterministically, the other preimage of each point on the way explored as usual); then step forward t < q times
    // into the transit's exit petal s, so that with L backward steps left, n = L mod q and m = (L - n)/q:
    // Ψ_s(Φ_s(y) + m) = f^{qm}(y), σ = Φ_s(y) + m - ζ0, D = (f^L)'(y)/Φ_s'(y), and the size estimate
    // C ≈ K A_card / |D² Φ_a'(v)|².  Paths are pruned when that estimate with the current derivative (an upper bound
    // once |y| is small) falls below cmin, or past dmax branching steps.  Prints "T 1 n re im C_est L".
    const double cmin = atof(argv[4]);
    const int dmax = atoi(argv[5]);
    // default r_loc: where psi trusts the series (|a_{-q}| r^{-q} ≈ 400)
    const double rloc = argc > 6 ? atof(argv[6]) : std::pow(std::hypot(L.a[0].r, L.a[0].i) / 400, 1.0 / q);
    Cd z0, dz0, ddz0;
    int s0;
    L.phi_a(L.v, z0, dz0, ddz0, s0);
    const int sx = L.exit_petal(s0);
    const double Acard = 3 * M_PI / 8, cv2 = std::norm(std::complex<double>(dz0.r, dz0.i));
    int64_t found = 0, pruned = 0;
    const double rswitch = getenv("GL_RSWITCH") ? atof(getenv("GL_RSWITCH")) : 0.3;
    // an upper bound for |Φ_rep'| at entry (|u| ≈ rswitch), with a safety factor 10
    const double logB = std::log(10 * q * std::hypot(L.a[0].r, L.a[0].i) * std::pow(rswitch, -q - 1));
    std::function<void(Cd, int, int, double)> dfs = [&](Cd y, int depth, int bdepth, double logP) {
      // depth = backward steps, bdepth = branching steps (bounded by dmax)
      // logP = log |(f^depth)'(y)|; near 0 in a repelling sector, follow the local inverse branch into the petal
      // (deterministic), and also explore the other preimage below (the local one too if the follow did not land)
      bool landed = false;
      if (std::hypot(y.r, y.i) < rswitch && L.petal(y, 1) >= 0) {
        Cd u = y;
        int dd = depth;
        double lp = logP;
        bool ok = true;
        while (std::hypot(u.r, u.i) >= rloc) {
          const SCd disc = std::sqrt(SCd(L.lam.r, L.lam.i) * SCd(L.lam.r, L.lam.i) + 4.0 * SCd(u.r, u.i));
          SCd a1 = (-SCd(L.lam.r, L.lam.i) + disc) / 2.0, a2 = (-SCd(L.lam.r, L.lam.i) - disc) / 2.0;
          const SCd nu = std::abs(a1) < std::abs(a2) ? a1 : a2, far = std::abs(a1) < std::abs(a2) ? a2 : a1;
          // paths that follow the local branch for a while and then leave by the other preimage
          if (dd > depth) dfs(Cd(far.real(), far.imag()), dd + 1, bdepth + 1, lp + std::log(std::abs(SCd(L.lam.r, L.lam.i) + 2.0 * far)));
          lp += std::log(std::abs(SCd(L.lam.r, L.lam.i) + 2.0 * nu));
          u = Cd(nu.real(), nu.imag());
          if (++dd > depth + 200000) { ok = false; break; }
        }
        int k = ok ? L.petal(u, 1) : -1;
        if (k >= 0) {
          landed = true;
          // step forward into the exit petal sx, where Ψ_sx(Φ_sx(u) + m) = f^{qm}(u) (psi's branch convention)
          int t = 0;
          for (; t < q && k != sx; t++) {
            lp -= std::log(std::abs(SCd(L.lam.r, L.lam.i) + 2.0 * SCd(u.r, u.i)));
            u = L.lam * u + sqr(u);
            k = L.petal(u, 1);
          }
          if (k != sx) return;
          Cd sr, sd, sdd;
          L.series(u, L.axis(1, k), sr, sd, sdd);
          const int n = (((dd - t) % q) + q) % q;
          const Cd sigma = sr + Cd(double((dd - t - n) / q)) - z0;
          const double logD = lp - std::log(std::hypot(sd.r, sd.i));
          const double C = K * Acard / (std::exp(4 * logD) * cv2);
          if (C > cmin) { printf("T 1 %d %.17g %.17g %.6e %d\n", n, sigma.r, sigma.i, C, dd - t); found++; }
        }
      }
      if (bdepth >= dmax) { pruned++; return; }
      if (std::log(K * Acard / cv2) + 4 * (logB - logP) < std::log(cmin)) { pruned++; return; }
      const SCd disc = std::sqrt(SCd(L.lam.r, L.lam.i) * SCd(L.lam.r, L.lam.i) + 4.0 * SCd(y.r, y.i));
      const bool near = landed;
      const SCd b1 = (-SCd(L.lam.r, L.lam.i) + disc) / 2.0, b2 = (-SCd(L.lam.r, L.lam.i) - disc) / 2.0;
      const SCd local = std::abs(b1) < std::abs(b2) ? b1 : b2;
      for (const int sg : {1, -1}) {
        const SCd yy = (-SCd(L.lam.r, L.lam.i) + double(sg) * disc) / 2.0;
        if (near && yy == local) continue;   // the local branch was followed into the petal above
        const Cd yc(yy.real(), yy.imag());
        const double fp = std::abs(SCd(L.lam.r, L.lam.i) + 2.0 * yy);
        dfs(yc, depth + 1, bdepth + 1, logP + std::log(fp));
      }
    };
    dfs(L.crit, 0, 0, 0.0);
    fprintf(stderr, "tree: %lld centers, %lld paths cut at depth %d\n", (long long)found, (long long)pruned, dmax);
    return 0;
  }
  if (mode == "locate") {
    char name[256];
    int r;
    double sr, si;
    while (scanf("%255s %d %lf %lf", name, &r, &sr, &si) == 4) {
      Cd t, d, dd, hp;
      if (!L.theta(Cd(sr, si), r, t, d, dd, hp)) { printf("%s %d failed\n", name, r); continue; }
      printf("%s %d %.15g %.15g %.15g %.15g %.15g %.15g\n", name, r, t.r, t.i, d.r, d.i, hp.r, hp.i);
    }
    return 0;
  }
  if (mode == "children") {
    // The r-transit children of the (r-1)-transit source U (center, excursion, radius) over the single-transit target
    // c (n_c < q, σ_c mod 1): Newton on Θ_r(σ) = σ_c + j/q from both roots of U's local quadratic (plus CHILDREN_STARTS
    // points on each of three circles), each verified as a center with excursion n_c - j (this also checks that the
    // last transit's petal matches j mod q).  Prints "name|±|j re im w sat", w = |Θ_r' Π H'|^-2 (C ≈ C_c w), sat = 1
    // within 4 radii of U.
    const int threads = argc > 4 ? atoi(argv[4]) : 2;
    const int m = getenv("CHILDREN_STARTS") ? atoi(getenv("CHILDREN_STARTS")) : 0;
    struct Q { std::string name; int r, nu, nc, jlo, jhi; Cd u, c; double rad; };
    std::vector<Q> qs;
    char name[256];
    int r, nu, nc, jlo, jhi;
    double ur, ui, rad, cr, ci;
    while (scanf("%255s %d %d %lf %lf %lf %d %lf %lf %d %d", name, &r, &nu, &ur, &ui, &rad, &nc, &cr, &ci, &jlo, &jhi) == 11)
      qs.push_back({name, r, nu, nc, jlo, jhi, Cd(ur, ui), Cd(cr, ci), rad});
    std::vector<std::string> out(qs.size());
    std::atomic<int64_t> next(0);
    std::atomic<int64_t> done(0);
    std::vector<std::thread> pool;
    for (int th = 0; th < threads; th++)
      pool.emplace_back([&]() {
        for (int64_t qi; (qi = next++) < int64_t(qs.size());) {
          const auto& Q = qs[qi];
          std::string text;
          Cd t0, d0, dd0, hp0;
          if (!L.theta(Q.u, Q.r, t0, d0, dd0, hp0)) continue;
          for (int j = Q.jlo; j <= Q.jhi; j++) {
            if (Q.nc - j < 0) continue;
            const SCd y(Q.c.r + double(j) / q, Q.c.i);
            const SCd a = 0.5 * SCd(dd0.r, dd0.i), b(d0.r, d0.i), cc = SCd(t0.r, t0.i) - y;
            const SCd sq = std::sqrt(b * b - 4.0 * a * cc);
            std::vector<SCd> starts;
            const SCd u(Q.u.r, Q.u.i);
            for (const int sgn : {1, -1}) starts.push_back(u + (-b + double(sgn) * sq) / (2.0 * a));
            const double rho = std::abs(starts[0] - u);
            for (const double f : {0.5, 1.0, 2.0})
              for (int k = 0; k < m; k++) starts.push_back(u + std::polar(f * rho, 2 * M_PI * (k + 0.5) / m));
            std::vector<SCd> found;
            for (size_t si = 0; si < starts.size(); si++) {
              SCd s = starts[si];
              bool ok = false;
              for (int it = 0; it < 60; it++) {
                Cd t, d, dd, hp;
                if (!L.theta(Cd(s.real(), s.imag()), Q.r, t, d, dd, hp)) break;
                const SCd step = (SCd(t.r, t.i) - y) / SCd(d.r, d.i);
                s -= step;
                if (!(std::abs(step) < 10)) break;
                if (std::abs(step) < 1e-11 * (1 + std::abs(s))) { ok = true; break; }
              }
              if (!ok) continue;
              bool dup = false;
              for (const auto& fk : found) dup |= std::abs(fk - s) < 1e-9;
              if (dup) continue;
              found.push_back(s);
              Cd cen(s.real(), s.imag());
              if (!L.center(Q.r, Q.nc - j, cen) || std::hypot(cen.r - s.real(), cen.i - s.imag()) > 1e-8) continue;
              Cd t, d, dd, hp;
              L.theta(Cd(s.real(), s.imag()), Q.r, t, d, dd, hp);
              const double w = 1 / std::norm(SCd(d.r, d.i) * SCd(hp.r, hp.i));
              const int sat = std::abs(s - u) < 4 * Q.rad;
              char line[512];
              snprintf(line, sizeof(line), "%s|%c|%d %.17g %.17g %.6e %d\n", Q.name.c_str(), si == 0 ? '+' : '-', j,
                       s.real(), s.imag(), w, sat);
              text += line;
            }
          }
          out[qi] = text;
          if (++done % 1000 == 0) fprintf(stderr, "children: %lld / %zu\n", (long long)done, qs.size());
        }
      });
    for (auto& t : pool) t.join();
    for (const auto& t : out) fputs(t.c_str(), stdout);
    return 0;
  }
  if (mode == "debug") {
    const Cd sig(atof(argv[4]), atof(argv[5]));
    Cd s0, d0, dd0;
    int k;
    const bool ok = L.phi_a(L.v, s0, d0, dd0, k);
    printf("phi_a(v) ok %d = %.12g%+.12gi petal %d (d %.6g%+.6gi)\n", ok, s0.r, s0.i, k, d0.r, d0.i);
    for (int kk = 0; kk < q; kk++) {
      Cd w, d, dd;
      const bool ok2 = L.psi(s0 + sig, kk, w, d, dd);
      if (kk == L.exit_petal(k)) printf("(exit petal) ");
      Cd fw = L.lam * w + w * w;
      printf("psi(ζ0+σ, %d) ok %d = %.12g%+.12gi; f of it %.6g%+.6gi (crit %.6g%+.6gi)\n", kk, ok2, w.r, w.i, fw.r, fw.i, L.crit.r, L.crit.i);
    }
    Cd sg = sig;
    printf("center(1,1) ok %d -> %.12g%+.12gi\n", L.center(1, 1, sg), sg.r, sg.i);
    return 0;
  }
  if (mode == "consist") {
    // f-equivariance of the per-petal conventions: Ψ_{k'}(ζ + 1/q) = f(Ψ_k(ζ)) (k' the petal of f(Ψ_k(ζ))), and the
    // attracting series Φ(f(w)) = Φ(w) + 1/q across petals
    for (int k = 0; k < q; k++) {
      const Cd z(-60, 0.7);
      Cd w, d, dd, w2, d2, dd2;
      L.psi(z, k, w, d, dd);
      const Cd fw = L.lam * w + w * w;
      const int k2 = L.petal(fw, 1);
      // the shift j 2πiβ/q (constant on the petal) relating the branches
      int bj = 0;
      double be = 1e300;
      for (int j = -2 * q; j <= 2 * q; j++) {
        L.psi(z + Cd(1.0 / q) + Cd(0, 2 * M_PI / q) * L.beta * Cd(double(j)), k2, w2, d2, dd2);
        const double e = std::hypot((fw - w2).r, (fw - w2).i);
        if (e < be) { be = e; bj = j; }
      }
      printf("repelling %d -> %d: f(Ψ_k(ζ)) = Ψ_k'(ζ + 1/q + %d 2πiβ/q) to %.3g (|w| %.3g)\n", k, k2, bj, be, std::hypot(w.r, w.i));
    }
    for (int k = 0; k < q; k++) {
      const double base = -std::arg(std::complex<double>(-L.A.r, -L.A.i)) / q + 2 * M_PI * k / q;
      const Cd w(0.01 * cos(base + 0.1), 0.01 * sin(base + 0.1));
      const int kk = L.petal(w, -1);
      const Cd fw = L.lam * w + w * w;
      const int k2 = L.petal(fw, -1);
      Cd s1, d1, e1, s2, d2, e2;
      L.series(w, L.axis(-1, kk), s1, d1, e1);
      L.series(fw, L.axis(-1, k2), s2, d2, e2);
      const Cd e = s2 - s1 - Cd(1.0 / q);
      const auto jb = std::complex<double>(e.r, e.i) / (std::complex<double>(0, 2 * M_PI / q) * std::complex<double>(L.beta.r, L.beta.i));
      printf("attracting %d -> %d: Φ(f(w)) - Φ(w) - 1/q = %.3g%+.3gi = %.6g%+.2gi 2πiβ/q; exit petals %d -> %d\n", kk, k2, e.r, e.i, jb.real(), jb.imag(), L.exit_petal(kk), L.exit_petal(k2));
    }
    return 0;
  }
  fprintf(stderr, "unknown mode %s\n", mode.c_str());
  return 1;
}
