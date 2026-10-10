// The general p/q-root Lavaurs model (glavaurs.h): components, areas, and children.
//
//   glavaurs p q single nmax re0 re1 im0 im1 grid threads   # single-transit centers by grid Newton, with areas
//   glavaurs p q area threads < "name r n re im"              # areas: "name r n cre cim area_σ C conv cusp C_nf"
//   glavaurs p q tree cmin dmax [rloc]                        # single-transit centers as backward paths
//   glavaurs p q children threads < "name r n_u u_re u_im radius n_c c_re c_im jmin jmax"   # r-transit children
//   glavaurs p q locate < "name r re im"                      # Θ_r(σ), Θ_r', Π H'
//   glavaurs p q walk r n count < "re im" seeds               # a family σ_k ≈ σ_0 + k (n fixed): centers and C_nf
//   glavaurs p q dchildren threads < (as children)            # children by the one-step map (linearized source)
//   glavaurs p q fchildren threads < (as children)            # children by the frozen-σ horn-map dynamics
//   glavaurs p q frefine threads < "name r y_re y_im sX_re sX_im sW_re sW_im"   # frozen-σ child from the exact one
//   glavaurs p q consist                                      # cross-petal branch conventions
//   glavaurs p q jetarea threads < "name r n re im"           # C_nf and the first-order shape-corrected area
//   glavaurs p q tunelabel X:p:cre:cim,... < "U r n re im"     # predicted tunings U*X from U's weight-4 jet
//   glavaurs p q hp r n re im [N Nb]                          # center, C_nf, area in double, Expansion<2>, Expansion<3>
//   glavaurs p q arc sgn                                      # the critical arc: Φ_a(w) = Φ_a(crit) + i sgn t
//   glavaurs p q orbit J < "name r n re im"                   # explicit critical-orbit points, for kneading
//   glavaurs p q tuned threads < "W r n re im U rU nU ure uim"   # little-Julia-set tuning test of W by U
//   glavaurs p q classify theta threads < "name r n re im"    # limbs by kneading: "name m offset" | "name bulb r"
// C = (4π² sin²(πp/q)/q⁴) area_σ is the family constant lim k⁴ area_M (limbs [CF(p/q), k]).
#include "glavaurs.h"
#include "expansion_arith.h"
#include <atomic>
#include <iostream>
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
static Cd cd_mul(const Cd a, const Cd b) { const SCd c = SCd(a.r, a.i) * SCd(b.r, b.i); return Cd(c.real(), c.imag()); }
static Cd cd_div(const Cd a, const Cd b) { const SCd c = SCd(a.r, a.i) / SCd(b.r, b.i); return Cd(c.real(), c.imag()); }

// Explicit critical-orbit points of an r-transit component (gate interiors omitted): per transit T = 1..r the steps
// until the attracting petal ('e', T - 1, s: T = 1 from the critical value, s = 1 being v; T > 1 after the landing point
// x_{T-1}), then the J exit steps before the landing point x_T and x_T itself ('x', T, -J..0), then the final
// excursion up to (not including) the critical point ('f', r, s).
struct OrbitPoint { char kind; int T, s; Cd w; };
static bool orbit_points(const GeneralLavaurs& L, const int r, const int n, const Cd sigma, const int J,
                         std::vector<OrbitPoint>& pts) {
  const int q = L.q;
  int dp;
  {
    // f maps repelling petal k to petal k + dp
    Cd w, d, dd;
    L.psi(Cd(-60, 0), 0, w, d, dd);
    dp = L.petal(L.lam * w + w * w, 1);
  }
  const Cd tau = cd_mul(Cd(0, 2 * M_PI / q), L.beta);
  Cd x = L.v;
  for (int T = 1; T <= r; T++) {
    Cd w = T == 1 ? x : L.lam * x + x * x;
    for (int s = 1; s < 1000000; s++) {
      if (std::hypot(w.r, w.i) < 0.05 && L.petal(w, -1) >= 0) break;
      pts.push_back({'e', T - 1, s, w});
      w = L.lam * w + w * w;
      if (std::hypot(w.r, w.i) > 10) return false;
    }
    Cd p0, d0, dd0;
    int pet;
    if (!L.phi_a(x, p0, d0, dd0, pet)) return false;
    const Cd zeta = p0 + sigma + L.transit_shift(pet);
    const int e = L.exit_petal(pet);
    // the predecessor in petal b = a - dp of Ψ_a(ζ) is Ψ_b(ζ - 1/q - (a - b) τ)
    std::vector<Cd> ex;
    Cd z = zeta;
    int a = e;
    for (int j = 1; j <= J; j++) {
      const int b = ((a - dp) % q + q) % q;
      z = z - Cd(1.0 / q) - cd_mul(Cd(double(a - b)), tau);
      Cd y, d, dd;
      if (!L.psi(z, b, y, d, dd)) return false;
      ex.push_back(y);
      a = b;
    }
    for (int j = J; j >= 1; j--) pts.push_back({'x', T, -j, ex[j - 1]});
    Cd d, dd;
    if (!L.psi(zeta, e, x, d, dd)) return false;
    pts.push_back({'x', T, 0, x});
  }
  Cd w = x;
  for (int s = 1; s < n; s++) {
    w = L.lam * w + w * w;
    pts.push_back({'f', r, s, w});
  }
  return true;
}

// Truncated Taylor series for the shape of a component (host, double): univariate T4 (degree 4) and bivariate B4 in
// (u, δ) with weight i + 2j ≤ 4 (u = w - crit, δ = σ - σ_c: near a center δ ~ u², so these are the consistent orders)
struct GLT4 { std::complex<double> c[5]; };
static inline GLT4 t4_mul(const GLT4& a, const GLT4& b) {
  GLT4 r{};
  for (int i = 0; i < 5; i++) for (int j = 0; i + j < 5; j++) r.c[i + j] += a.c[i] * b.c[j];
  return r;
}
// f ∘ g for g with zero constant term: Σ_k f_k g^k
static inline GLT4 t4_compose(const GLT4& f, const GLT4& g) {
  GLT4 r{}, p{};
  p.c[0] = 1;
  for (int k = 0; k < 5; k++) {
    for (int i = 0; i < 5; i++) r.c[i] += f.c[k] * p.c[i];
    p = t4_mul(p, g);
  }
  return r;
}
struct GLB4 {
  // monomials u^i δ^j with i + 2j ≤ 4
  static constexpr int M = 9;
  static constexpr int I[M] = {0, 1, 2, 3, 4, 0, 1, 2, 0}, J[M] = {0, 0, 0, 0, 0, 1, 1, 1, 2};
  std::complex<double> c[M];
  static int index(const int i, const int j) {
    for (int m = 0; m < M; m++) if (I[m] == i && J[m] == j) return m;
    return -1;
  }
};
static inline GLB4 b4_mul(const GLB4& a, const GLB4& b) {
  GLB4 r{};
  for (int x = 0; x < GLB4::M; x++)
    for (int y = 0; y < GLB4::M; y++) {
      const int k = GLB4::index(GLB4::I[x] + GLB4::I[y], GLB4::J[x] + GLB4::J[y]);
      if (k >= 0) r.c[k] += a.c[x] * b.c[y];
    }
  return r;
}
// g(x) for the univariate Taylor series g about x's constant term
static inline GLB4 b4_apply(const GLT4& g, const GLB4& x) {
  GLB4 d = x, r{}, p{};
  d.c[0] = 0;
  p.c[0] = 1;
  for (int k = 0; k < 5; k++) {
    for (int m = 0; m < GLB4::M; m++) r.c[m] += g.c[k] * p.c[m];
    p = b4_mul(p, d);
  }
  return r;
}

// Taylor coefficients (Φ^{(k)}/k!, k ≤ 4) of the Fatou series at w on the branch of axis ax
static inline GLT4 gl_series_t4(const GLCore& L, const std::complex<double> w, const double ax) {
  typedef std::complex<double> SC;
  GLT4 t{};
  Complex<double> s, d, dd;
  L.series(Complex<double>(w.real(), w.imag()), ax, s, d, dd);
  t.c[0] = SC(s.r, s.i);
  const SC beta(L.beta.r, L.beta.i);
  // β log: derivatives (-1)^{k-1} (k-1)!/w^k; w^j: j (j-1)..(j-k+1) w^{j-k}
  for (int k = 1; k <= 4; k++) {
    double fact = 1; for (int i = 2; i < k; i++) fact *= i;   // (k-1)!
    SC v = beta * ((k % 2 ? 1.0 : -1.0) * fact) / std::pow(w, k);
    for (int j = -L.q; j <= L.N; j++) {
      if (!j) continue;
      double ff = 1; for (int i = 0; i < k; i++) ff *= (j - i);
      v += SC(L.a[j + L.q].r, L.a[j + L.q].i) * ff * std::pow(w, j - k);
    }
    double kf = 1; for (int i = 2; i <= k; i++) kf *= i;
    t.c[k] = v / kf;
  }
  return t;
}
static inline GLT4 t4_f(const GLCore& L, const GLT4& g) {   // λ g + g²
  GLT4 r = t4_mul(g, g);
  const std::complex<double> lam(L.lam.r, L.lam.i);
  for (int i = 0; i < 5; i++) r.c[i] += lam * g.c[i];
  return r;
}
// Φ_a(w0 + h) as a series in h, and the entering petal
static inline bool gl_phi_a_t4(const GLCore& L, const std::complex<double> w0, GLT4& out, int& pet) {
  GLT4 g{};
  g.c[0] = w0; g.c[1] = 1;
  for (int n = 0; n < (1 << 20); n++) {
    const std::complex<double> w = g.c[0];
    if (std::abs(w) < L.r0) {
      const int k = L.petal(Complex<double>(w.real(), w.imag()), -1);
      if (k >= 0) {
        const GLT4 P = gl_series_t4(L, w, L.axis(-1, k));
        GLT4 dlt = g; dlt.c[0] = 0;
        out = t4_compose(P, dlt);
        out.c[0] -= double(n) / L.q;
        const Complex<double> off = L.hit_offset(k);
        out.c[0] += std::complex<double>(off.r, off.i);
        pet = L.label(k, n);
        return true;
      }
    }
    g = t4_f(L, g);
    if (std::abs(g.c[0]) > 10) return false;
  }
  return false;
}
// Ψ_k(ζ0 + h) as a series in h
static inline bool gl_psi_t4(const GLCore& L, const std::complex<double> z0, const int k, GLT4& out) {
  typedef std::complex<double> SC;
  Complex<double> w, d, dd;
  if (!L.psi(Complex<double>(z0.real(), z0.imag()), k, w, d, dd)) return false;   // validity and the same branch
  const double mm = std::ceil(z0.real() + std::max(GLCore::kRepel, 2 * std::fabs(z0.imag())));
  const long long m = mm > 0 ? (long long)mm : 0;
  const SC zl = z0 - double(m);
  // u with Φ(u) = zl: Newton from psi's own choice (rerun the local solve)
  const double ax = L.axis(1, k);
  const double base = -L.argA / L.q + 2 * M_PI * k / L.q;
  const SC ratio = SC(L.a[0].r, L.a[0].i) / zl;
  const SC rr = std::polar(std::pow(std::abs(ratio), 1.0 / L.q), std::arg(ratio) / L.q);
  SC u = rr;
  double bd = 10;
  for (int t = 0; t < L.q; t++) {
    const SC uu = rr * std::polar(1.0, 2 * M_PI * t / L.q);
    const double dl = std::fabs(std::remainder(std::arg(uu) - base, 2 * M_PI));
    if (dl < bd) { bd = dl; u = uu; }
  }
  for (int it = 0; it < 80; it++) {
    Complex<double> s, sd, sdd;
    L.series(Complex<double>(u.real(), u.imag()), ax, s, sd, sdd);
    const SC step = (SC(s.r, s.i) - zl) / SC(sd.r, sd.i);
    u -= step;
    if (std::abs(step) < 1e-16 * std::abs(u)) break;
  }
  const GLT4 P = gl_series_t4(L, u, ax);
  // reversion: Φ(u + Δ) = zl + h, Δ = Σ b_k h^k
  const SC p1 = P.c[1], p2 = P.c[2], p3 = P.c[3], p4 = P.c[4];
  const SC b1 = 1.0 / p1, b2 = -p2 * b1 * b1 / p1, b3 = -(2.0 * p2 * b1 * b2 + p3 * b1 * b1 * b1) / p1;
  const SC b4 = -(p2 * (b2 * b2 + 2.0 * b1 * b3) + 3.0 * p3 * b1 * b1 * b2 + p4 * b1 * b1 * b1 * b1) / p1;
  GLT4 g{};
  g.c[0] = u; g.c[1] = b1; g.c[2] = b2; g.c[3] = b3; g.c[4] = b4;
  for (long long i = 0; i < L.q * m; i++) {
    g = t4_f(L, g);
    if (std::abs(g.c[0]) > 10) return false;
  }
  out = g;
  return true;
}
// The return map R_{σ_c + δ}(crit + u) as a bivariate series
static inline bool gl_return_b4(const GLCore& L, const int r, const int n, const std::complex<double> sc, GLB4& x) {
  typedef std::complex<double> SC;
  const SC lam(L.lam.r, L.lam.i);
  x = GLB4{};
  x.c[0] = SC(L.crit.r, L.crit.i); x.c[GLB4::index(1, 0)] = 1;
  const auto F = [&]() { GLB4 y = b4_mul(x, x); for (int m = 0; m < GLB4::M; m++) y.c[m] += lam * x.c[m]; x = y; };
  F();
  for (int t = 0; t < r; t++) {
    GLT4 P;
    int pet;
    if (!gl_phi_a_t4(L, x.c[0], P, pet)) return false;
    x = b4_apply(P, x);
    x.c[0] += sc;
    x.c[GLB4::index(0, 1)] += 1.0;
    const Complex<double> sh = L.transit_shift(pet);
    GLT4 Q;
    if (!gl_psi_t4(L, x.c[0] + SC(sh.r, sh.i), L.exit_petal(pet), Q)) return false;
    x = b4_apply(Q, x);
  }
  for (int i = 0; i < n; i++) F();
  return true;
}

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
          if (L.area(j.r, j.n, j.g, cen, ar, conv, cusp, a1)) {
            // the normal-form constant at the center: R(w) ≈ crit + D δσ + A (w - crit)², a cardioid of area 3π/8 in
            // c = A D δσ, so C_nf = K (3π/8)/|A D|²: what the component's children inherit (not its exact area)
            GLJet x;
            double cnf = NAN;
            if (L.core.return_map(j.r, j.n, L.crit, cen, x))
              cnf = K * (3 * M_PI / 8) / std::norm(SCd(x.ww.r / 2, x.ww.i / 2) * SCd(x.s.r, x.s.i));
            snprintf(line, sizeof(line), "%s %d %d %.17g %.17g %.10e %.10e %.1e %.3e %.10e", j.name.c_str(), j.r, j.n,
                     cen.r, cen.i, ar, K * ar, conv, cusp, cnf);
          }
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
  if (mode == "walk") {
    // Continue a family of r-transit components with excursion n whose members sit at σ_k + k (long excursions: the
    // representative (σ - k, n + q r k)): from two members, predict by polynomial extrapolation, Newton for the center,
    // print "k re im C_nf" (the normal-form constant: the family's exact masses to the shape factor)
    const int r = atoi(argv[4]), n = atoi(argv[5]), count = atoi(argv[6]);
    std::vector<Cd> sig;
    double sr, si;
    while (scanf("%lf %lf", &sr, &si) == 2) { Cd c(sr, si); if (L.center(r, n, c)) sig.push_back(c); else break; }
    if (sig.size() < 2) { fprintf(stderr, "walk: need two seeds\n"); return 1; }
    double lastA = NAN, lastD = NAN;
    const auto cnf_at = [&](const Cd c) {
      GLJet x;
      if (!L.core.return_map(r, n, L.crit, c, x)) return double(NAN);
      lastA = std::abs(SCd(x.ww.r / 2, x.ww.i / 2)); lastD = std::abs(SCd(x.s.r, x.s.i));
      return K * (3 * M_PI / 8) / std::norm(SCd(x.ww.r / 2, x.ww.i / 2) * SCd(x.s.r, x.s.i));
    };
    for (int k = 0; k < int(sig.size()); k++) { const double cn = cnf_at(sig[k]); printf("%d %.17g %.17g %.10e %.6e %.6e\n", k, sig[k].r, sig[k].i, cn, lastA, lastD); }
    for (int k = int(sig.size()); k < count; k++) {
      // extrapolate the drift d_k = σ_k - k from up to 4 previous members (Lagrange in k)
      const int m = std::min(4, k);
      SCd pred = 0;
      for (int a = 0; a < m; a++) {
        const int ka = k - 1 - a;
        SCd l = 1;
        for (int b = 0; b < m; b++) if (b != a) { const int kb = k - 1 - b; l *= double(k - kb) / double(ka - kb); }
        pred += l * SCd(sig[ka].r - ka, sig[ka].i);
      }
      Cd c(pred.real() + k, pred.imag());
      if (!L.center(r, n, c)) { fprintf(stderr, "walk: no center at k = %d\n", k); break; }
      if (std::hypot((c - sig[k - 1]).r - 1, (c - sig[k - 1]).i) > 0.5) { fprintf(stderr, "walk: jumped at k = %d\n", k); break; }
      sig.push_back(c);
      const double cn = cnf_at(c);
      printf("%d %.17g %.17g %.10e %.6e %.6e\n", k, c.r, c.i, cn, lastA, lastD);
      fflush(stdout);
    }
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
  if (mode == "frefine") {
    // For an exact child σ_W of source σ_X (Θ_r(σ_W) = y): Newton on the frozen-σ Θ_r (σ frozen at σ_X after the first
    // transit) from σ_W.  Prints "name sf_re sf_im w_exact wT wmult" (w = |Θ' Π H'|^-2 exact at σ_W; the frozen child's
    // weight with the parameter derivative along the frozen orbit, and |Π H'|^-4), or "name failed".
    const int threads = argc > 4 ? atoi(argv[4]) : 2;
    // FR_T=1: the slice crosses the fiber at the source's transversality T_{r-1}(σ_X) = Θ'_{r-1}/Π H'_{r-1}
    const bool useT = getenv("FR_T") && atoi(getenv("FR_T"));
    struct Q { std::string name; int r; Cd y, sx, sw; };
    std::vector<Q> qs;
    char name[256];
    int r;
    double a1, a2, b1, b2, c1, c2;
    while (scanf("%255s %d %lf %lf %lf %lf %lf %lf", name, &r, &a1, &a2, &b1, &b2, &c1, &c2) == 8)
      qs.push_back({name, r, Cd(a1, a2), Cd(b1, b2), Cd(c1, c2)});
    std::vector<std::string> out(qs.size());
    std::atomic<int64_t> next(0);
    std::vector<std::thread> pool;
    for (int th = 0; th < threads; th++)
      pool.emplace_back([&]() {
        for (int64_t i; (i = next++) < int64_t(qs.size());) {
          const auto& Q = qs[i];
          char line[512];
          Cd t, d, dd, hp;
          snprintf(line, sizeof(line), "%s failed", Q.name.c_str());
          if (L.theta(Q.sw, Q.r, t, d, dd, hp)) {
            const double we = 1 / std::norm(SCd(d.r, d.i) * SCd(hp.r, hp.i));
            Cd speed(1);
            if (useT) {
              Cd t1, d1, dd1, h1;
              if (L.theta(Q.sx, Q.r - 1, t1, d1, dd1, h1)) speed = cd_div(d1, h1);
            }
            SCd s(Q.sw.r, Q.sw.i);
            bool ok = false;
            Cd tf, dz, dq, hq, pr;
            int pet;
            for (int it = 0; it < 60; it++) {
              if (!L.core.theta_frozen(Cd(s.real(), s.imag()), Q.sx, Q.r, tf, dz, dq, hq, pr, pet, speed)) break;
              const SCd step = (SCd(tf.r, tf.i) - SCd(Q.y.r, Q.y.i)) / SCd(dz.r, dz.i);
              s -= step;
              if (!(std::abs(step) < 10)) break;
              if (std::abs(step) < 1e-11 * (1 + std::abs(s))) { ok = true; break; }
            }
            if (ok) {
              L.core.theta_frozen(Cd(s.real(), s.imag()), Q.sx, Q.r, tf, dz, dq, hq, pr, pet, speed);
              const SCd H(hq.r, hq.i);
              // The holomorphic motion of the fiber chain in δ (F_δ = H + σ_X + δ): z(δ) = z0 + z1 δ + z2 δ², with the
              // target fixed and z_k^{(i)} from the next point's by H' z1 + 1 = z1', H' z2 + H''/2 z1² = z2'; the child is
              // the crossing ζ0 + σ_X + δ = z^{(1)}(δ), solved to first and second order (FR_T unset: z0 = ζ0 + σ_f)
              Cd h1[kMaxChain], h2[kMaxChain], prr;
              const Cd z0 = L.zeta0() + Cd(s.real(), s.imag());
              SCd d1 = NAN, d2 = NAN;
              if (!useT && gl_fiber_chain(L.core, z0, Q.sx, Q.r, h1, h2, prr)) {
                SCd z1 = 0, z2 = 0;
                for (int i = Q.r - 1; i >= 1; i--) {
                  const SCd a(h1[i].r, h1[i].i), b(h2[i].r, h2[i].i);
                  const SCd n1 = (z1 - 1.0) / a, n2 = (z2 - 0.5 * b * n1 * n1) / a;
                  z1 = n1; z2 = n2;
                }
                const SCd rhs = SCd(z0.r, z0.i) - SCd(L.zeta0().r + Q.sx.r, L.zeta0().i + Q.sx.i);
                d1 = rhs / (1.0 - z1);
                // second order: (1 - z1) δ - z2 δ² = rhs, by Newton from d1
                SCd dd = d1;
                for (int it = 0; it < 20; it++) dd -= ((1.0 - z1) * dd - z2 * dd * dd - rhs) / ((1.0 - z1) - 2.0 * z2 * dd);
                d2 = dd;
              }
              const SCd sx(Q.sx.r, Q.sx.i);
              snprintf(line, sizeof(line), "%s %.17g %.17g %.6e %.6e %.6e %.17g %.17g %.17g %.17g", Q.name.c_str(), s.real(),
                       s.imag(), we, 1 / std::norm(SCd(dq.r, dq.i) * H), 1 / std::norm(H * H), (sx + d1).real(),
                       (sx + d1).imag(), (sx + d2).real(), (sx + d2).imag());
            }
          }
          out[i] = line;
        }
      });
    for (auto& t : pool) t.join();
    for (const auto& l : out) printf("%s\n", l.c_str());
    return 0;
  }
  if (mode == "fchildren") {
    // Children by the frozen-σ dynamics: with σ frozen at the source's σ_X in every transit after the first, the
    // r-transit condition is F^{r-1}(p_1) ∈ targets for the single map F(p) = H(p) + σ_X, p_1 = ζ0 + σ exactly linear
    // in σ (no linearization of the source's parameter map).  Newton on the frozen Θ_r from the same starts as
    // children; kept when the last landing point reaches the critical point in n_c - j steps.  Prints "name|±|j re im w
    // 0 sat" with w = |dpar Π H'|^-2 (FCH_T=1, default: the parameter derivative along the frozen orbit) or |Π H'|^-4
    // (FCH_T=0: purely multiplicative).
    const int threads = argc > 4 ? atoi(argv[4]) : 2;
    const int m = getenv("CHILDREN_STARTS") ? atoi(getenv("CHILDREN_STARTS")) : 0;
    const bool useT = !getenv("FCH_T") || atoi(getenv("FCH_T"));
    struct Q { std::string name; int r, nu, nc, jlo, jhi; Cd u, c; double rad; };
    std::vector<Q> qs;
    char name[256];
    int r, nu, nc, jlo, jhi;
    double ur, ui, rad, cr, ci;
    while (scanf("%255s %d %d %lf %lf %lf %d %lf %lf %d %d", name, &r, &nu, &ur, &ui, &rad, &nc, &cr, &ci, &jlo, &jhi) == 11)
      qs.push_back({name, r, nu, nc, jlo, jhi, Cd(ur, ui), Cd(cr, ci), rad});
    std::vector<std::string> out(qs.size());
    std::atomic<int64_t> next(0);
    std::vector<std::thread> pool;
    for (int th = 0; th < threads; th++)
      pool.emplace_back([&]() {
        for (int64_t qi; (qi = next++) < int64_t(qs.size());) {
          const auto& Q = qs[qi];
          std::string text;
          Cd t0, d0, dd0, hp0;
          if (!L.theta(Q.u, Q.r, t0, d0, dd0, hp0)) continue;   // the starts, as in children
          const SCd u(Q.u.r, Q.u.i);
          for (int j = Q.jlo; j <= Q.jhi; j++) {
            if (Q.nc - j < 0) continue;
            const Cd ts = L.core.target_shift(j);
            const SCd y = SCd(Q.c.r + double(j) / q, Q.c.i) + SCd(ts.r, ts.i);
            const SCd a = 0.5 * SCd(dd0.r, dd0.i), b(d0.r, d0.i), cc = SCd(t0.r, t0.i) - y;
            const SCd sq = std::sqrt(b * b - 4.0 * a * cc);
            std::vector<SCd> starts;
            for (const int sgn : {1, -1}) starts.push_back(u + (-b + double(sgn) * sq) / (2.0 * a));
            const double rho = std::abs(starts[0] - u);
            for (const double f : {0.5, 1.0, 2.0})
              for (int k = 0; k < m; k++) starts.push_back(u + std::polar(f * rho, 2 * M_PI * (k + 0.5) / m));
            std::vector<SCd> found;
            for (size_t si = 0; si < starts.size(); si++) {
              SCd s = starts[si];
              bool ok = false;
              Cd tf, dz, dq, hp, pr;
              int pet = -1;
              for (int it = 0; it < 60; it++) {
                if (!L.core.theta_frozen(Cd(s.real(), s.imag()), Q.u, Q.r, tf, dz, dq, hp, pr, pet)) break;
                const SCd step = (SCd(tf.r, tf.i) - y) / SCd(dz.r, dz.i);
                s -= step;
                if (!(std::abs(step) < 10)) break;
                if (std::abs(step) < 1e-11 * (1 + std::abs(s))) { ok = true; break; }
              }
              if (!ok) continue;
              bool dup = false;
              for (const auto& fk : found) dup |= std::abs(fk - s) < 1e-9;
              if (dup) continue;
              found.push_back(s);
              L.core.theta_frozen(Cd(s.real(), s.imag()), Q.u, Q.r, tf, dz, dq, hp, pr, pet);
              Cd x, d1, d2;
              if (!L.psi(pr, L.exit_petal(pet), x, d1, d2)) continue;
              for (int i = 0; i < Q.nc - j; i++) x = L.lam * x + x * x;
              if (!(std::hypot((x - L.crit).r, (x - L.crit).i) < 1e-6)) continue;
              const SCd H(hp.r, hp.i);
              const double w = useT ? 1 / std::norm(SCd(dq.r, dq.i) * H) : 1 / std::norm(H * H);
              const int sat = std::abs(s - u) < 4 * Q.rad;
              char line[512];
              snprintf(line, sizeof(line), "%s|%c|%d %.17g %.17g %.6e 0 %d\n", Q.name.c_str(), si == 0 ? '+' : '-', j,
                       s.real(), s.imag(), w, sat);
              text += line;
            }
          }
          out[qi] = text;
        }
      });
    for (auto& t : pool) t.join();
    for (const auto& t : out) fputs(t.c_str(), stdout);
    return 0;
  }
  if (mode == "dchildren") {
    // The one-step (dynamical) version of children: at the source X, p* = p_{r-1}(σ_X) is a critical point of the horn
    // map and to first order p_{r-1}(σ) = p* + Θ'_{r-1}(σ - σ_X), so the children solve F(p) = H(p) + ε (p - p*) + σ_X -
    // ζ0 = σ_c + j/q with ε = 1/Θ'_{r-1}(σ_X) (exact at r = 2, where ε = 1: the ζ formulation), and Θ_r' = Θ'_{r-1}
    // (H'(p) + ε), Π H'_r = Π H'_{r-1} H'(p): one transfer-operator step of p ↦ F(p) with the extra state ε, which
    // evolves as ε ↦ ε/(H' + ε).  DCH_EPS=0: ε = 0 (frozen σ, the purely multiplicative weight).  Solutions are kept
    // when the next landing point reaches the critical point in n_c - j steps.  Prints "name|±|j sw_re sw_im w eps' sat"
    // with σ_W = σ_X + ε (p - p*), w = |Θ_r' Π H'_r|^-2, eps' the child's ε, sat within 4 radii.
    const int threads = argc > 4 ? atoi(argv[4]) : 2;
    const int m = getenv("CHILDREN_STARTS") ? atoi(getenv("CHILDREN_STARTS")) : 0;
    const Cd z0 = L.zeta0();
    const bool use_eps = !getenv("DCH_EPS") || atoi(getenv("DCH_EPS"));
    struct Q { std::string name; int r, nu, nc, jlo, jhi; Cd u, c; double rad; };
    std::vector<Q> qs;
    char name[256];
    int r, nu, nc, jlo, jhi;
    double ur, ui, rad, cr, ci;
    while (scanf("%255s %d %d %lf %lf %lf %d %lf %lf %d %d", name, &r, &nu, &ur, &ui, &rad, &nc, &cr, &ci, &jlo, &jhi) == 11)
      qs.push_back({name, r, nu, nc, jlo, jhi, Cd(ur, ui), Cd(cr, ci), rad});
    std::vector<std::string> out(qs.size());
    std::atomic<int64_t> next(0);
    std::vector<std::thread> pool;
    for (int th = 0; th < threads; th++)
      pool.emplace_back([&]() {
        for (int64_t qi; (qi = next++) < int64_t(qs.size());) {
          const auto& Q = qs[qi];
          std::string text;
          Cd t0, d0, dd0, hp0;
          int pet;
          if (!L.theta(Q.u, Q.r - 1, t0, d0, dd0, hp0, &pet)) continue;
          const SCd eps = use_eps ? 1.0 / SCd(d0.r, d0.i) : SCd(0);
          const SCd ps = SCd(t0.r + z0.r, t0.i + z0.i), Tp(d0.r, d0.i), Hp(hp0.r, hp0.i), u(Q.u.r, Q.u.i);
          Cd h0, dh0, ddh0;
          int pet2;
          if (!L.horn(Cd(ps.real(), ps.imag()), pet, h0, dh0, ddh0, pet2)) continue;
          for (int j = Q.jlo; j <= Q.jhi; j++) {
            if (Q.nc - j < 0) continue;
            const Cd ts = L.core.target_shift(j);
            const SCd target = SCd(Q.c.r + double(j) / q, Q.c.i) + SCd(ts.r, ts.i) + SCd(z0.r, z0.i) - u;
            // local quadratic H(p*) + ε u + a u² = target (H'(p*) = 0)
            const SCd a = 0.5 * SCd(ddh0.r, ddh0.i), cc = SCd(h0.r, h0.i) - target;
            const SCd sq = std::sqrt(eps * eps - 4.0 * a * cc);
            std::vector<SCd> starts = {ps + (-eps + sq) / (2.0 * a), ps + (-eps - sq) / (2.0 * a)};
            const double rho = std::abs(starts[0] - ps);
            for (const double f : {0.5, 1.0, 2.0})
              for (int k = 0; k < m; k++) starts.push_back(ps + std::polar(f * rho, 2 * M_PI * (k + 0.5) / m));
            std::vector<SCd> found;
            for (size_t si = 0; si < starts.size(); si++) {
              SCd pp = starts[si];
              bool ok = false;
              Cd h, dh, ddh;
              int pe2 = -1;
              for (int it = 0; it < 60; it++) {
                if (!L.horn(Cd(pp.real(), pp.imag()), pet, h, dh, ddh, pe2)) break;
                const SCd step = (SCd(h.r, h.i) + eps * (pp - ps) - target) / (SCd(dh.r, dh.i) + eps);
                pp -= step;
                if (!(std::abs(step) < 10)) break;
                if (std::abs(step) < 1e-11 * (1 + std::abs(pp))) { ok = true; break; }
              }
              if (!ok) continue;
              bool dup = false;
              for (const auto& fk : found) dup |= std::abs(fk - pp) < 1e-9;
              if (dup) continue;
              found.push_back(pp);
              L.horn(Cd(pp.real(), pp.imag()), pet, h, dh, ddh, pe2);
              // the next landing point must reach the critical point in n_c - j steps (the petal combinatorics)
              Cd x, d1, d2;
              const SCd sm = u + eps * (pp - ps), sw = u + (pp - ps) / Tp;   // the model's σ; the predicted child
              if (!L.psi(Cd(h.r + sm.real(), h.i + sm.imag()), L.exit_petal(pe2), x, d1, d2)) continue;
              for (int i = 0; i < Q.nc - j; i++) x = L.lam * x + x * x;
              if (!(std::hypot((x - L.crit).r, (x - L.crit).i) < 1e-6)) continue;
              const SCd H1(dh.r, dh.i);
              const double w = 1 / std::norm(Tp * (H1 + eps) * Hp * H1);
              const SCd eps2 = eps / (H1 + eps);
              const int sat = std::abs(sw - u) < 4 * Q.rad;
              char line[512];
              snprintf(line, sizeof(line), "%s|%c|%d %.17g %.17g %.6e %.6e %d\n", Q.name.c_str(), si == 0 ? '+' : '-',
                       j, sw.real(), sw.imag(), w, std::abs(eps2), sat);
              text += line;
            }
          }
          out[qi] = text;
        }
      });
    for (auto& t : pool) t.join();
    for (const auto& t : out) fputs(t.c_str(), stdout);
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
    // On the GPU when there is one (GL_GPU=0: CPU threads); same algorithm and output
    if (!getenv("GL_GPU") || atoi(getenv("GL_GPU"))) {
      std::vector<GLChildJob> jobs;
      for (const auto& Q : qs) jobs.push_back(GLChildJob{Q.r, Q.nc, Q.jlo, Q.jhi, Q.u, Q.c, Q.rad});
      std::vector<GLChild> kids;
      if (glavaurs_gpu_children(L.core, jobs, m, kids)) {
        for (const auto& k : kids)
          printf("%s|%c|%d %.17g %.17g %.6e %d\n", qs[k.job].name.c_str(), k.start == 0 ? '+' : '-', k.j, k.s.r, k.s.i,
                 k.w, k.sat);
        return 0;
      }
    }
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
            const Cd ts = L.core.target_shift(j);
            const SCd y = SCd(Q.c.r + double(j) / q, Q.c.i) + SCd(ts.r, ts.i);
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
  if (mode == "arc") {
    // The lift under Φ_a (2:1 at the critical point) of the vertical half-line {Φ_a(crit) + i sgn t, t > 0}: two
    // branches from the critical point, ending at α (w = 0) and -α (w = -λ).  Prints "A branch re im".
    const double sgn = atof(argv[4]);
    Cd p0, d0, dd0;
    int k;
    if (!L.phi_a(L.crit, p0, d0, dd0, k)) return 1;
    for (const int br : {1, -1}) {
      double t = 1e-4;
      const SCd w0 = SCd(L.crit.r, L.crit.i) + double(br) * std::sqrt(2.0 * SCd(0, sgn * t) / SCd(dd0.r, dd0.i));
      Cd w(w0.real(), w0.imag());
      for (int it = 0; it < 4000 && t < 1e4; it++) {
        bool ok = false;
        for (int nt = 0; nt < 60; nt++) {
          Cd s, d, dd;
          int kk;
          if (!L.phi_a(w, s, d, dd, kk)) break;
          const SCd step = (SCd(s.r, s.i) - SCd(p0.r, p0.i + sgn * t)) / SCd(d.r, d.i);
          w = w - Cd(step.real(), step.imag());
          // phi_a runs ~1/kR0^q iterations: residuals bottom out near 1e-12
          if (std::abs(step) < 1e-9 && std::abs(step * SCd(d.r, d.i)) < 1e-10) { ok = true; break; }
        }
        if (!ok) break;
        printf("A %d %.17g %.17g\n", br, w.r, w.i);
        if (std::hypot(w.r, w.i) < 0.01 || std::hypot(w.r + L.lam.r, w.i + L.lam.i) < 0.01) break;
        t *= 1.03;
      }
    }
    return 0;
  }
  if (mode == "orbit") {
    // The critical orbit of an r-transit component as explicit points (orbit_points); prints "name kind T s re im"
    const int J = atoi(argv[4]);
    char name[256];
    int r, n;
    double sr, si;
    while (scanf("%255s %d %d %lf %lf", name, &r, &n, &sr, &si) == 5) {
      std::vector<OrbitPoint> pts;
      if (!orbit_points(L, r, n, Cd(sr, si), J, pts)) { printf("%s failed\n", name); continue; }
      for (const auto& o : pts) printf("%s %c %d %d %.17g %.17g\n", name, o.kind, o.T, o.s, o.w.r, o.w.i);
    }
    return 0;
  }
  if (mode == "tuned") {
    // W ∈ U*M iff at σ_W the critical orbit of U's return map R_U (r_U transits, excursion n_U) stays in U's little
    // filled Julia set through all p = r/r_U returns (the p-th return is W's own and closes the cycle; requires
    // n + 1 = p (n_U + 1) with σ_W in U's representative).  In the normal form R_U(w) ≈ crit + D + a (w - crit)² the
    // little K lies in |w - crit| ≤ 2/|a|; prints "W U ratio" with ratio = max_{j<p} |R_U^j(crit) - crit| |a|/2
    // (≤ 1 for tunings; ∞ if the orbit fails).
    const int threads = argc > 4 ? atoi(argv[4]) : 2;
    struct Job { std::string w, u; int r, n, ru, nu; Cd sw, su; };
    std::vector<Job> jobs;
    char wn[256], un[256];
    int r, n, ru, nu;
    double a, b, c, d;
    while (scanf("%255s %d %d %lf %lf %255s %d %d %lf %lf", wn, &r, &n, &a, &b, un, &ru, &nu, &c, &d) == 10)
      jobs.push_back({wn, un, r, n, ru, nu, Cd(a, b), Cd(c, d)});
    std::vector<std::string> out(jobs.size());
    std::atomic<int64_t> next(0);
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; t++)
      pool.emplace_back([&]() {
        for (int64_t i; (i = next++) < int64_t(jobs.size());) {
          const auto& jb = jobs[i];
          char line[600];
          double ratio = INFINITY;
          if (jb.r % jb.ru == 0 && jb.n + 1 == (jb.r / jb.ru) * (jb.nu + 1)) {
            const int pp = jb.r / jb.ru;
            Cd v, dw, dww;
            if (L.return_map(jb.ru, jb.nu, L.crit, jb.sw, v, dw, dww)) {
              const double R = 2 / std::hypot(dww.r / 2, dww.i / 2);
              double mx = std::hypot((v - L.crit).r, (v - L.crit).i);
              Cd x = v;
              bool ok = true;
              for (int j = 2; j < pp && ok; j++) {
                ok = L.return_map(jb.ru, jb.nu, x, jb.sw, v, dw, dww);
                x = v;
                mx = std::max(mx, std::hypot((x - L.crit).r, (x - L.crit).i));
                if (!(mx < 100 * R)) ok = false;
              }
              if (ok) ratio = mx / R;
            }
          }
          snprintf(line, sizeof(line), "%s %s %.4g", jb.w.c_str(), jb.u.c_str(), ratio);
          out[i] = line;
        }
      });
    for (auto& t : pool) t.join();
    for (const auto& l : out) printf("%s\n", l.c_str());
    return 0;
  }
  if (mode == "classify") {
    // The limb of each component (gl_limb.py's kneading classifier in C++): partition R_{θ/2} ∪ A ∪ R_{(θ+1)/2} with
    // θ the root's parameter angle on the gate's side (gate +1: the lower angle, gate -1: the upper) and A the arc
    // Φ_a = Φ_a(crit) + i gate t; ν = 1 on the side of the angles (θ/2, (θ+1)/2) and within 0.15 of α.  Prints
    // "name m offset" at the first 0 (offset from the landing point x_m; ≡ -1 mod q, b = (offset + 1)/q mod m), or
    // "name bulb r" (no 0 before the critical point), or "name failed".
    const double theta = atof(argv[4]);
    const int threads = argc > 5 ? atoi(argv[5]) : 2;
    const bool knead = getenv("KNEAD") && atoi(getenv("KNEAD"));   // print the whole kneading word instead
    const int J = 40;
    const SCd lam(L.lam.r, L.lam.i), c0 = lam / 2.0 - lam * lam / 4.0, alpha = lam / 2.0;
    const auto ray = [&](const double t) {
      std::vector<SCd> pts;
      const double ER = 1e4;
      SCd z = std::polar(std::sqrt(ER), 2 * M_PI * t);
      for (int n = 1; n <= 40; n++) {
        const double ang = 2 * M_PI * std::fmod(std::ldexp(t, n), 1.0);
        for (int j = 1; j <= 16; j++) {
          const SCd target = std::polar(std::pow(ER, std::pow(2.0, -j / 16.0)), ang);
          for (int it = 0; it < 60; it++) {
            SCd w = z, dw = 1;
            for (int m = 0; m < n; m++) { dw = 2.0 * w * dw; w = w * w + c0; }
            const SCd step = (w - target) / dw;
            z -= step;
            if (std::abs(step) < 1e-15 * (1 + std::abs(z))) break;
          }
          pts.push_back(z);
        }
      }
      return pts;
    };
    // the arc's branches, in z coordinates
    std::vector<SCd> br[2];
    {
      const double sgn = L.side;
      Cd p0, d0, dd0;
      int k;
      if (!L.phi_a(L.crit, p0, d0, dd0, k)) return 1;
      for (int b = 0; b < 2; b++) {
        double t = 1e-4;
        const SCd w0 = SCd(L.crit.r, L.crit.i) + (b ? -1.0 : 1.0) * std::sqrt(2.0 * SCd(0, sgn * t) / SCd(dd0.r, dd0.i));
        Cd w(w0.real(), w0.imag());
        for (int it = 0; it < 4000 && t < 1e4; it++) {
          bool ok = false;
          for (int nt = 0; nt < 60; nt++) {
            Cd s, d, dd;
            int kk;
            if (!L.phi_a(w, s, d, dd, kk)) break;
            const SCd step = (SCd(s.r, s.i) - SCd(p0.r, p0.i + sgn * t)) / SCd(d.r, d.i);
            w = w - Cd(step.real(), step.imag());
            if (std::abs(step) < 1e-9 && std::abs(step * SCd(d.r, d.i)) < 1e-10) { ok = true; break; }
          }
          if (!ok) break;
          br[b].push_back(SCd(w.r, w.i) + alpha);
          if (std::abs(br[b].back() - alpha) < 0.01 || std::abs(br[b].back() + alpha) < 0.01) break;
          t *= 1.03;
        }
      }
    }
    const int ba = std::abs(br[0].back() - alpha) < std::abs(br[1].back() - alpha) ? 0 : 1;
    if (!(std::abs(br[ba].back() - alpha) < 0.05 && std::abs(br[1 - ba].back() + alpha) < 0.05)) {
      fprintf(stderr, "classify: the arc does not reach ±α\n");
      return 1;
    }
    const double a1 = theta / 2, a2 = (theta + 1) / 2;
    const auto r1 = ray(a1), r2 = ray(a2);
    const bool r1a = std::abs(r1.back() - alpha) < std::abs(r2.back() - alpha);
    const auto& arc1 = r1a ? br[ba] : br[1 - ba];
    const auto& arc2 = r1a ? br[1 - ba] : br[ba];
    std::vector<SCd> poly;
    const double far = 200;
    poly.push_back(std::polar(far, 2 * M_PI * a1));
    poly.insert(poly.end(), r1.begin(), r1.end());
    poly.insert(poly.end(), arc1.rbegin(), arc1.rend());
    poly.push_back(SCd(L.crit.r, L.crit.i) + alpha);
    poly.insert(poly.end(), arc2.begin(), arc2.end());
    poly.insert(poly.end(), r2.rbegin(), r2.rend());
    for (int d = 0; d < 720; d++) poly.push_back(std::polar(far, 2 * M_PI * (a2 - (a2 - a1) * d / 720)));
    const auto inside = [&](const SCd z) {
      bool in = false;
      for (size_t i = 0; i < poly.size(); i++) {
        const SCd a = poly[i], b = poly[(i + 1) % poly.size()];
        if ((a.imag() > z.imag()) != (b.imag() > z.imag())) {
          const double x = a.real() + (z.imag() - a.imag()) * (b.real() - a.real()) / (b.imag() - a.imag());
          if (x > z.real()) in = !in;
        }
      }
      return in;
    };
    struct Job { std::string name; int r, n; Cd s; };
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
          const auto& jb = jobs[i];
          std::vector<OrbitPoint> pts;
          char line[512];
          if (!orbit_points(L, jb.r, jb.n, jb.s, J, pts)) { snprintf(line, sizeof(line), "%s failed", jb.name.c_str()); out[i] = line; continue; }
          if (knead) {
            // the whole explicit kneading word: e-points, 'G' at each gate passage (its interior is all 1s), the exit
            // steps and landing point, ..., the final excursion (the critical point itself omitted)
            // points within 0.15 of α (gate steps, always 1) and gate passages collapse into one marker '|', so only
            // the symbols of points away from α remain: a word independent of how many steps a passage takes
            std::string w = jb.name + " ";
            for (const auto& o : pts) {
              const bool near = std::hypot(o.w.r, o.w.i) < 0.15 || (o.kind == 'x' && o.s < 0);
              if (near) { if (w.back() != '|') w += '|'; }
              else w += inside(SCd(o.w.r, o.w.i) + alpha) ? '1' : '0';
            }
            out[i] = w;
            continue;
          }
          snprintf(line, sizeof(line), "%s bulb %d", jb.name.c_str(), jb.r);
          for (const auto& o : pts) {
            if (std::hypot(o.w.r, o.w.i) < 0.15) continue;
            if (!inside(SCd(o.w.r, o.w.i) + alpha)) {
              if (o.kind == 'e' && o.T == 0) snprintf(line, sizeof(line), "%s pre %d", jb.name.c_str(), o.s);
              else snprintf(line, sizeof(line), "%s %d %d", jb.name.c_str(), o.T, o.s);
              break;
            }
          }
          out[i] = line;
        }
      });
    for (auto& t : pool) t.join();
    for (const auto& l : out) printf("%s\n", l.c_str());
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
  if (mode == "jetarea") {
    // At a center the return map is R = crit + D δ + A u² + E u δ + F δ² + B u³ + ... (u = w - crit, δ = σ - σ_c,
    // R_w = 0).  The attracting fixed point with multiplier μ: to zeroth order u0 = μ/(2A), δ0 = (μ/2 - μ²/4)/(AD)
    // (the cardioid: area_σ = (3π/8)/|AD|², C_nf); to first order in (B, E, F) u1 = -(3B u0² + E δ0)/(2A) and
    // D δ1 = (1 - μ) u1 - B u0³ - E δ0 u0 - F δ0², a polynomial in μ of degree 4, whose area is π Σ k |a_k|² exactly.
    // Prints "name C_nf C_jet1 |B|/|A|^2·|...| scale" (the first-order area constant and a nonlinearity measure).
    const int threads = argc > 4 ? atoi(argv[4]) : 2;
    struct Job { std::string name; int r, n; Cd s; };
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
          const auto& jb = jobs[i];
          char line[512];
          snprintf(line, sizeof(line), "%s failed", jb.name.c_str());
          Cd c = jb.s;
          GLJet3<double> x;
          if (L.center(jb.r, jb.n, c) && gl_return_map3(L.core, jb.r, jb.n, L.crit, c, x)) {
            typedef std::vector<SCd> P;   // polynomial coefficients in μ
            const auto mul = [](const P& a, const P& b) { P c(a.size() + b.size() - 1, 0.0); for (size_t i = 0; i < a.size(); i++) for (size_t j = 0; j < b.size(); j++) c[i + j] += a[i] * b[j]; return c; };
            const auto add = [](P a, const P& b, const SCd f) { if (a.size() < b.size()) a.resize(b.size(), 0.0); for (size_t i = 0; i < b.size(); i++) a[i] += f * b[i]; return a; };
            const SCd A = 0.5 * SCd(x.ww.r, x.ww.i), D(x.s.r, x.s.i), E(x.ws.r, x.ws.i), F = 0.5 * SCd(x.ss.r, x.ss.i),
                      B = SCd(x.www.r, x.www.i) / 6.0;
            const P mu = {0.0, 1.0};
            const P u0 = {0.0, 1.0 / (2.0 * A)};
            const P d0 = {0.0, 0.5 / (A * D), -0.25 / (A * D)};
            P u1 = add(add(P{}, mul(u0, u0), -3.0 * B / (2.0 * A)), d0, -E / (2.0 * A));
            P Dd1 = mul(P{1.0, -1.0}, u1);
            Dd1 = add(Dd1, mul(mul(u0, u0), u0), -B);
            Dd1 = add(Dd1, mul(d0, u0), -E);
            Dd1 = add(Dd1, mul(d0, d0), -F);
            const P dl = add(d0, Dd1, 1.0 / D);
            const auto area = [](const P& a) { double s = 0; for (size_t k = 1; k < a.size(); k++) s += k * std::norm(a[k]); return M_PI * s; };
            const double cnf = K * area(d0), cj1 = K * area(dl);
            // nonlinearity: the first-order correction's size relative to the cardioid's
            double num = 0, den = 0;
            for (size_t k = 1; k < dl.size(); k++) { num += std::norm(dl[k] - (k < d0.size() ? d0[k] : 0.0)); den += std::norm(k < d0.size() ? d0[k] : 0.0); }
            // the weight-4 truncated model: solve R(u, δ) = crit + u, ∂_u R = μ on |μ| = 1 (radially, then around)
            double cb4 = NAN;
            GLB4 X;
            const bool b4ok = gl_return_b4(L.core, jb.r, jb.n, SCd(c.r, c.i), X);
            if (getenv("JET_DEBUG")) {
              fprintf(stderr, "b4 %d:", b4ok);
              for (int m = 0; m < GLB4::M; m++) fprintf(stderr, " (%d,%d) %.3e%+.3ei", GLB4::I[m], GLB4::J[m], X.c[m].real(), X.c[m].imag());
              fprintf(stderr, "\n jet3: A %.3e D %.3e E %.3e F %.3e B %.3e\n", std::abs(A), std::abs(D), std::abs(E), std::abs(F), std::abs(B));
            }
            if (b4ok) {
              const auto ev = [&](const SCd u, const SCd d, SCd& R, SCd& Ru, SCd& Rd, SCd& Ruu, SCd& Rud) {
                R = Ru = Rd = Ruu = Rud = 0;
                // integer powers by multiplication (std::pow on complex 0 gives nan)
                SCd up[6] = {1, u, u * u, u * u * u, u * u * u * u, 0}, dp[4] = {1, d, d * d, 0};
                const auto pw = [&](const SCd* t, const int e) { return e >= 0 ? t[e] : SCd(0); };
                for (int m = 0; m < GLB4::M; m++) {
                  const int i = GLB4::I[m], j = GLB4::J[m];
                  const SCd cm = X.c[m];
                  R += cm * pw(up, i) * pw(dp, j);
                  if (i >= 1) Ru += cm * double(i) * pw(up, i - 1) * pw(dp, j);
                  if (j >= 1) Rd += cm * double(j) * pw(up, i) * pw(dp, j - 1);
                  if (i >= 2) Ruu += cm * double(i * (i - 1)) * pw(up, i - 2) * pw(dp, j);
                  if (i >= 1 && j >= 1) Rud += cm * double(i * j) * pw(up, i - 1) * pw(dp, j - 1);
                }
              };
              const SCd cr(L.crit.r, L.crit.i);
              const int Nb = 64;
              std::vector<SCd> pts(Nb);
              SCd u = 0, d = 0;
              bool ok = true;
              const auto solve = [&](const SCd mu) {
                for (int it = 0; it < 50; it++) {
                  SCd R, Ru, Rd, Ruu, Rud;
                  ev(u, d, R, Ru, Rd, Ruu, Rud);
                  const SCd F1 = R - cr - u, F2 = Ru - mu;
                  const SCd a11 = Ru - 1.0, a12 = Rd, a21 = Ruu, a22 = Rud, det = a11 * a22 - a12 * a21;
                  const SCd du = (F1 * a22 - a12 * F2) / det, dd = (a11 * F2 - a21 * F1) / det;
                  u -= du; d -= dd;
                  if (std::abs(du) + std::abs(dd) < 1e-12 * (std::abs(u) + std::abs(d))) return true;
                }
                return false;
              };
              for (int i = 1; i <= 16 && ok; i++) ok = solve(std::polar(double(i) / 16, M_PI / Nb));
              for (int j = 0; j < Nb && ok; j++) {
                if (j) for (int t = 1; t < 4 && ok; t++) ok = solve(std::polar(1.0, M_PI * (2 * j - 1 + 2.0 * t / 4) / Nb));
                if (ok) ok = solve(std::polar(1.0, M_PI * (2 * j + 1) / Nb));
                pts[j] = d;
              }
              if (ok) {
                double sum = 0;
                for (int k = 1; k < Nb; k++) {
                  SCd ak = 0;
                  for (int j = 0; j < Nb; j++) ak += pts[j] * std::polar(1.0, -M_PI * double((int64_t(k) * (2 * j + 1)) % (2 * Nb)) / Nb);
                  sum += k * std::norm(ak);
                }
                cb4 = K * M_PI * sum / (double(Nb) * Nb);
              }
            }
            snprintf(line, sizeof(line), "%s %.10e %.10e %.3e %.10e", jb.name.c_str(), cnf, cj1, std::sqrt(num / den), cb4);
          }
          out[i] = line;
        }
      });
    for (auto& t : pool) t.join();
    for (const auto& l : out) printf("%s\n", l.c_str());
    return 0;
  }
  if (mode == "tunelabel") {
    // U*X for primitive X of period p (center c_X in M): in U's weight-4 jet model R(crit + u, σ_U + δ), solve for the
    // δ at which the model's critical point is p-periodic, from the normal-form guess δ = c_X/(A D), by Newton on
    // F(δ) = u_p(δ) (u_{i+1} = R(crit + u_i, σ_U + δ) - crit, u_0 = 0); then refine the center of the r_U p-transit,
    // excursion p (n_U + 1) - 1 component exactly from σ_U + δ.  Prints "U*X model_re model_im center_re center_im"
    // (or "failed").
    struct TX { std::string name; int p; SCd c; };
    std::vector<TX> xs;
    {
      const std::string t = argv[4];
      size_t i = 0;
      while (i < t.size()) {
        size_t j = t.find(',', i);
        if (j == std::string::npos) j = t.size();
        char nm[64];
        int pp;
        double cr, ci;
        if (sscanf(t.substr(i, j - i).c_str(), "%63[^:]:%d:%lf:%lf", nm, &pp, &cr, &ci) == 4) xs.push_back({nm, pp, SCd(cr, ci)});
        i = j + 1;
      }
    }
    char name[256];
    int r, n;
    double sr, si;
    while (scanf("%255s %d %d %lf %lf", name, &r, &n, &sr, &si) == 5) {
      Cd c(sr, si);
      if (!L.center(r, n, c)) { printf("%s center failed\n", name); continue; }
      GLB4 X;
      if (!gl_return_b4(L.core, r, n, SCd(c.r, c.i), X)) { printf("%s jet failed\n", name); continue; }
      const SCd cr(L.crit.r, L.crit.i), A = X.c[GLB4::index(2, 0)], D = X.c[GLB4::index(0, 1)];
      // the model's map u ↦ R(crit + u, σ_U + δ) - crit and its δ-derivative
      const auto step = [&](const SCd u, const SCd du, const SCd d, SCd& un, SCd& dun) {
        SCd R = 0, Ru = 0, Rd = 0;
        SCd up[5] = {1, u, u * u, u * u * u, u * u * u * u}, dp[3] = {1, d, d * d};
        for (int m = 0; m < GLB4::M; m++) {
          const int i = GLB4::I[m], j = GLB4::J[m];
          R += X.c[m] * up[i] * dp[j];
          if (i >= 1) Ru += X.c[m] * double(i) * up[i - 1] * dp[j];
          if (j >= 1) Rd += X.c[m] * double(j) * up[i] * dp[j - 1];
        }
        un = R - cr; dun = Ru * du + Rd;
      };
      for (const auto& tx : xs) {
        SCd d = tx.c / (A * D);
        bool ok = false;
        for (int it = 0; it < 60; it++) {
          // Newton on G = u_p / Π_{d | p, d < p} u_d: deflates δ = 0 (U itself) and lower periods
          SCd u = 0, du = 0, logd = 0;
          for (int k = 1; k <= tx.p; k++) {
            SCd un, dun;
            step(u, du, d, un, dun);
            u = un; du = dun;
            if (k < tx.p && tx.p % k == 0) logd -= du / u;
          }
          logd += du / u;   // G'/G
          const SCd st = 1.0 / logd;
          d -= st;
          if (!(std::abs(st) < 10)) break;
          if (std::abs(st) < 1e-13 * std::abs(d)) { ok = true; break; }
        }
        if (!ok) { printf("%s*%s model failed\n", name, tx.name.c_str()); continue; }
        const SCd dmodel = d;
        // homotopy from the model to the exact return map: u_{k+1} = (1 - t)(R_model - crit) + t (R_exact - crit)
        const int rU = r, nU = n;
        for (int ti = 1; ti <= 20 && ok; ti++) {
          const double t = ti / 20.0;
          bool conv = false;
          for (int it = 0; it < 40; it++) {
            SCd u = 0, du = 0, logd = 0;
            bool good = true;
            for (int k = 1; k <= tx.p && good; k++) {
              SCd um, dum;
              step(u, du, d, um, dum);
              GLJet xj;
              const SCd w = cr + u, sg = SCd(c.r, c.i) + d;
              if (!L.core.return_map(rU, nU, Cd(w.real(), w.imag()), Cd(sg.real(), sg.imag()), xj)) { good = false; break; }
              const SCd ue = SCd(xj.v.r, xj.v.i) - cr, due = SCd(xj.w.r, xj.w.i) * du + SCd(xj.s.r, xj.s.i);
              u = (1 - t) * um + t * ue;
              du = (1 - t) * dum + t * due;
              if (k < tx.p && tx.p % k == 0) logd -= du / u;
            }
            if (!good) break;
            logd += du / u;
            const SCd st = 1.0 / logd;
            d -= st;
            if (!(std::abs(st) < 1)) break;
            if (std::abs(st) < 1e-8 * std::abs(d) + 1e-14) { conv = true; break; }
          }
          ok = conv;
        }
        if (!ok) { printf("%s*%s homotopy failed\n", name, tx.name.c_str()); continue; }
        (void)dmodel;
        const SCd sm = SCd(c.r, c.i) + d;
        Cd cw(sm.real(), sm.imag());
        const int rw = r * tx.p, nw = tx.p * (n + 1) - 1;
        if (L.center(rw, nw, cw))
          printf("%s*%s %d %d %.17g %.17g %.17g %.17g\n", name, tx.name.c_str(), rw, nw, sm.real(), sm.imag(), cw.r, cw.i);
        else printf("%s*%s %d %d %.17g %.17g center failed\n", name, tx.name.c_str(), rw, nw, sm.real(), sm.imag());
      }
    }
    return 0;
  }
  if (mode == "hp") {
    // A component's center and normal-form constant at three precisions, each from the previous one's center (the
    // Fatou series with N terms beyond double): agreement to each precision's level checks the high-precision core
    const int r = atoi(argv[4]), n = atoi(argv[5]), Nh = argc > 8 ? atoi(argv[8]) : 40, Nb = argc > 9 ? atoi(argv[9]) : 0;
    Cd c(atof(argv[6]), atof(argv[7]));
    if (!L.center(r, n, c)) { printf("double: no center\n"); return 1; }
    const auto run = [&](auto core, auto tol, const char* name, auto& center) {
      typedef decltype(core.lam) C;
      typedef decltype(center.r) S;
      C w = core.crit, sg = center;
      const double last = core.newton(r, n, w, sg, C(0), tol, 60);
      GLJetT<S> x;
      core.return_map(r, n, core.crit, sg, x);
      const C ad = C(S(0.5)) * x.ww * x.s;
      // K = 4π² sin²(πp/q)/q⁴ at precision S
      S sn, cs;
      gl_sincos(gl_pi<S>() * S(double(p)) / S(double(q)), sn, cs);
      const S pi = gl_pi<S>(), KS = S(4.0) * pi * pi * sn * sn / (S(double(q)) * S(double(q)) * S(double(q)) * S(double(q)));
      const S cnf = KS * S(3.0) * pi / S(8.0) / (ad.r * ad.r + ad.i * ad.i);
      std::cout << name << ": last step " << last << "\n  center " << safe(sg.r) << " " << safe(sg.i) << "\n  C_nf " << safe(cnf) << "\n";
      center = sg;
      if (Nb > 0) {
        C cen;
        S ar;
        double conv, cusp;
        if (gl_area(core, r, n, sg, cen, ar, conv, cusp, Nb, tol, 1e3 * tol))
          std::cout << "  area_σ " << safe(ar) << "\n  C = K area_σ " << safe(KS * ar) << "  (conv " << conv << ", cusp " << cusp << ")\n";
        else std::cout << "  area failed\n";
      }
    };
    printf("double (N = %d): center %.17g %.17g\n", L.N, c.r, c.i);
    Complex<Expansion<2>> c2(Expansion<2>(c.r), Expansion<2>(c.i));
    run(gl_make_core<Expansion<2>>(p, q, L.side, Nh), 1e-30, "Expansion<2>", c2);
    Complex<Expansion<3>> c3(c2.r.x[0] ? Expansion<3>(c2.r.x[0]) + Expansion<3>(c2.r.x[1]) : Expansion<3>(0.0),
                             c2.i.x[0] ? Expansion<3>(c2.i.x[0]) + Expansion<3>(c2.i.x[1]) : Expansion<3>(0.0));
    run(gl_make_core<Expansion<3>>(p, q, L.side, Nh), 1e-44, "Expansion<3>", c3);
    return 0;
  }
  if (mode == "geq") {
    // f-equivariance of the transit g_σ(w) = Ψ_e(Φ_a(w) + σ + shift(k)): |g(f(w)) - f(g(w))| along the critical orbit
    const Cd sigma(atof(argv[4]), atof(argv[5]));
    const auto g = [&](const Cd w, int& k) {
      Cd s0, d0, dd0, x, d, dd;
      if (!L.phi_a(w, s0, d0, dd0, k)) return Cd(NAN, NAN);
      if (!L.psi(s0 + sigma + L.transit_shift(k), L.exit_petal(k), x, d, dd)) return Cd(NAN, NAN);
      return x;
    };
    Cd w = L.v;
    for (int i = 0; i < 8; i++) {
      int k1, k2;
      const Cd a = g(L.lam * w + w * w, k2), gb = g(w, k1), b = L.lam * gb + gb * gb;
      printf("w%d: entering petals %d -> %d, |g(f w) - f(g w)| = %.3e  (|g(f w)| %.3g)\n", i, k1, k2, std::hypot((a - b).r, (a - b).i), std::hypot(a.r, a.i));
      w = L.lam * w + w * w;
    }
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
