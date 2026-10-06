// Areas of filled Julia sets K(c) from the area transfer operator

#include "julia.h"
#include "arb_cc.h"
#include "debug.h"
#include "expansion_arith.h"
#include "nearest.h"
#include "print.h"
#include <flint/acb.h>
#include <flint/arb_hypgeom.h>
#include <atomic>
#include <chrono>
#include <cmath>
#include <functional>
#include <thread>
namespace mandelbrot {

using std::function;
using std::max;
using std::sqrt;
using std::abs;

static const int prec = 320;  // Bits for arb geometry, well past Expansion<3>

// Run f(i) for i < n on all cores, in chunks
static void parallel_for(const int64_t n, const function<void(int64_t)>& f) {
  const int threads = max(1u, std::thread::hardware_concurrency());
  std::atomic<int64_t> next(0);
  const int64_t chunk = max(int64_t(1), n / (16 * threads));
  vector<std::thread> pool;
  for (int t = 0; t < threads; t++)
    pool.emplace_back([&]() {
      for (;;) {
        const int64_t i0 = next.fetch_add(chunk);
        if (i0 >= n) return;
        for (int64_t i = i0; i < std::min(n, i0 + chunk); i++) f(i);
      }
    });
  for (auto& t : pool) t.join();
}

template<class S> static S rnd(const arb_t x) { return round_near<S>(arb_midref(x)); }

// sin and cos of π p / q
static void sin_cos_pi(arb_t s, arb_t c, const int64_t p, const int64_t q) {
  Arb x;
  arb_set_si(x, p);
  arb_div_si(x, x, q, prec);
  arb_sin_cos_pi(s, c, x, prec);
}

template<class S> void cheb_row(const vector<S>& nodes, const vector<S>& weights, const S s, S* row) {
  const int n = int(nodes.size());
  S sum(0);
  for (int k = 0; k < n; k++) {
    const S d = s - nodes[k];
    if (d == S(0)) {
      for (int l = 0; l < n; l++) row[l] = l == k ? S(1) : S(0);
      return;
    }
    row[k] = weights[k] / d;
    sum += row[k];
  }
  const S inv_sum = inv(sum);
  for (int k = 0; k < n; k++) row[k] = row[k] * inv_sum;
}

template<class S> void trig_row(const vector<S>& rho_x, const vector<S>& rho_y, const S ex, const S ey, S* row) {
  // sin((θ - θ_j)/2) = Im(η conj(ρ_j)), with ρ_j = e^{iθ_j/2}
  const int n = int(rho_x.size());
  S sum(0);
  for (int j = 0; j < n; j++) {
    const S sn = ey * rho_x[j] - ex * rho_y[j];
    if (sn == S(0)) {
      for (int l = 0; l < n; l++) row[l] = l == j ? S(1) : S(0);
      return;
    }
    row[j] = j & 1 ? -inv(sn) : inv(sn);
    sum += row[j];
  }
  const S inv_sum = inv(sum);
  for (int j = 0; j < n; j++) row[j] = row[j] * inv_sum;
}

// For double, plain division
template<> void cheb_row(const vector<double>& nodes, const vector<double>& weights, const double s, double* row) {
  const int n = int(nodes.size());
  double sum = 0;
  for (int k = 0; k < n; k++) {
    const double d = s - nodes[k];
    if (d == 0) { for (int l = 0; l < n; l++) row[l] = l == k; return; }
    sum += row[k] = weights[k] / d;
  }
  for (int k = 0; k < n; k++) row[k] /= sum;
}
template<> void trig_row(const vector<double>& rho_x, const vector<double>& rho_y, const double ex, const double ey,
                         double* row) {
  const int n = int(rho_x.size());
  double sum = 0;
  for (int j = 0; j < n; j++) {
    const double sn = ey * rho_x[j] - ex * rho_y[j];
    if (sn == 0) { for (int l = 0; l < n; l++) row[l] = l == j; return; }
    sum += row[j] = (j & 1 ? -1 : 1) / sn;
  }
  for (int j = 0; j < n; j++) row[j] /= sum;
}

// y = L h for the factored operator: y_i = W_i Σ_k R[i,k] Σ_j T[i,j] H[k,j]
template<class S> static void apply(const int nr, const int nt, const vector<S>& W, const vector<S>& R,
                                    const vector<S>& T, const vector<S>& h, vector<S>& y) {
  const int64_t N = int64_t(nr) * nt;
  parallel_for(N, [&](const int64_t i) {
    const S* r = &R[i * nr];
    const S* t = &T[i * nt];
    S sum(0);
    for (int k = 0; k < nr; k++) {
      const S* hk = &h[int64_t(k) * nt];
      S dot(0);
      for (int j = 0; j < nt; j++) dot += t[j] * hk[j];
      sum += r[k] * dot;
    }
    y[i] = W[i] * sum;
  });
}

// Restarted GMRES for (1 - L) x = b in double; returns iterations
static int gmres(const int nr, const int nt, const vector<double>& W, const vector<double>& R,
                 const vector<double>& T, const vector<double>& b, vector<double>& x, const double tol) {
  const int64_t N = int64_t(nr) * nt;
  const int m = 60;
  const auto A = [&](const vector<double>& v, vector<double>& out) {
    apply(nr, nt, W, R, T, v, out);
    for (int64_t i = 0; i < N; i++) out[i] = v[i] - out[i];
  };
  const auto norm = [](const vector<double>& v) { double s = 0; for (const double a : v) s += a * a; return sqrt(s); };
  const double bnorm = norm(b);
  x.assign(N, 0);
  if (bnorm == 0) return 0;
  vector<double> r(N), w(N);
  int iters = 0;
  for (int restart = 0; restart < 100; restart++) {
    A(x, w);
    for (int64_t i = 0; i < N; i++) r[i] = b[i] - w[i];
    const double beta = norm(r);
    if (beta <= tol * bnorm) return iters;
    vector<vector<double>> V(1, r);
    for (auto& a : V[0]) a /= beta;
    vector<vector<double>> H(m + 1, vector<double>(m, 0));
    vector<double> cs(m), sn(m), g(m + 1, 0);
    g[0] = beta;
    int k = 0;
    for (; k < m; k++) {
      iters++;
      A(V[k], w);
      for (int l = 0; l <= k; l++) {  // Modified Gram-Schmidt, twice
        double d = 0;
        for (int64_t i = 0; i < N; i++) d += w[i] * V[l][i];
        H[l][k] = d;
        for (int64_t i = 0; i < N; i++) w[i] -= d * V[l][i];
      }
      for (int l = 0; l <= k; l++) {
        double d = 0;
        for (int64_t i = 0; i < N; i++) d += w[i] * V[l][i];
        H[l][k] += d;
        for (int64_t i = 0; i < N; i++) w[i] -= d * V[l][i];
      }
      H[k + 1][k] = norm(w);
      V.push_back(w);
      if (H[k + 1][k]) for (auto& a : V[k + 1]) a /= H[k + 1][k];
      for (int l = 0; l < k; l++) {  // Apply previous rotations, then a new one
        const double t = cs[l] * H[l][k] + sn[l] * H[l + 1][k];
        H[l + 1][k] = -sn[l] * H[l][k] + cs[l] * H[l + 1][k];
        H[l][k] = t;
      }
      const double d = std::hypot(H[k][k], H[k + 1][k]);
      cs[k] = H[k][k] / d;
      sn[k] = H[k + 1][k] / d;
      H[k][k] = d;
      H[k + 1][k] = 0;
      g[k + 1] = -sn[k] * g[k];
      g[k] = cs[k] * g[k];
      if (abs(g[k + 1]) <= tol * bnorm) { k++; break; }
    }
    vector<double> y(k);
    for (int l = k - 1; l >= 0; l--) {
      double s = g[l];
      for (int q = l + 1; q < k; q++) s -= H[l][q] * y[q];
      y[l] = s / H[l][l];
    }
    for (int l = 0; l < k; l++)
      for (int64_t i = 0; i < N; i++) x[i] += y[l] * V[l][i];
  }
  die("gmres did not converge");
}

template<class S> JuliaResult<S> julia_area(const JuliaParams& p) {
  const auto t0 = std::chrono::steady_clock::now();
  const int nr = p.nr, nt = p.nt;
  const int64_t N = int64_t(nr) * nt;
  slow_assert(nt % 2 == 1 && nr >= 2, "need odd nt and nr ≥ 2");
  const double ac = std::hypot(p.cx, p.cy);
  const double r1 = p.r1 ? p.r1 : max(2 * ac, 0.25), r2 = p.r2;
  slow_assert(r1 - ac > r1 * r1 && r2 * r2 - r2 > ac, "annulus %g < |z| < %g does not work for |c| = %g", r1, r2, ac);
  // Preimages of A have radii in [√(r1 - |c|), √(r2 + |c|)], strictly inside (r1, r2)
  slow_assert(sqrt(r1 - ac) > r1 && sqrt(r2 + ac) < r2);

  // Chebyshev nodes in s = log r, barycentric weights (-1)^k sin((2k+1)π/(2nr)), and half angles ρ_j = e^{iπj/nt}
  Arb s1, s2, mid, half_len, a, b;
  arb_set_d(a, r1); arb_log(s1, a, prec);
  arb_set_d(a, r2); arb_log(s2, a, prec);
  arb_add(mid, s1, s2, prec); arb_mul_2exp_si(mid, mid, -1);
  arb_sub(half_len, s2, s1, prec); arb_mul_2exp_si(half_len, half_len, -1);
  vector<Arb> sk(nr);
  vector<S> nodes(nr), weights(nr), rho_x(nt), rho_y(nt);
  vector<double> nodes_d(nr), weights_d(nr), rho_xd(nt), rho_yd(nt);
  for (int k = 0; k < nr; k++) {
    sin_cos_pi(a, b, 2 * k + 1, 2 * nr);
    arb_mul(sk[k], half_len, b, prec);
    arb_add(sk[k], sk[k], mid, prec);
    if (k & 1) arb_neg(a, a);
    nodes[k] = rnd<S>(sk[k]); weights[k] = rnd<S>(a);
    nodes_d[k] = double(nodes[k]); weights_d[k] = double(weights[k]);
  }
  for (int j = 0; j < nt; j++) {
    sin_cos_pi(a, b, j, nt);
    rho_x[j] = rnd<S>(b); rho_y[j] = rnd<S>(a);
    rho_xd[j] = double(rho_x[j]); rho_yd[j] = double(rho_y[j]);
  }

  // Per collocation point z = e^{s_k} e^{iθ_j}: W = 1/(4|z - c|), the preimage log radius ½ log|z - c|, and
  // η = e^{i arg(z - c)/4}, the half angle phase of √(z - c) (the other preimage has iη).  Rows in S and double.
  vector<S> W(N), R(N * nr), T(N * nt);
  vector<double> Wd(N), Rd(N * nr), Td(N * nt);
  parallel_for(N, [&](const int64_t i) {
    const int k = int(i / nt), j = int(i % nt);
    Arb r, sn, cs, m, phi, t;
    acb_t u;
    acb_init(u);
    arb_exp(r, sk[k], prec);
    sin_cos_pi(sn, cs, 2 * j, nt);
    arb_mul(acb_realref(u), r, cs, prec);
    arb_mul(acb_imagref(u), r, sn, prec);
    arb_sub(acb_realref(u), acb_realref(u), exact_arb(p.cx), prec);
    arb_sub(acb_imagref(u), acb_imagref(u), exact_arb(p.cy), prec);
    acb_abs(m, u, prec);
    acb_arg(phi, u, prec);
    acb_clear(u);
    arb_mul_2exp_si(t, m, 2);
    arb_inv(t, t, prec);
    W[i] = rnd<S>(t);
    arb_log(t, m, prec);
    arb_mul_2exp_si(t, t, -1);
    const S sw = rnd<S>(t);
    arb_mul_2exp_si(phi, phi, -2);
    arb_sin_cos(sn, cs, phi, prec);
    const S ex = rnd<S>(cs), ey = rnd<S>(sn);
    cheb_row(nodes, weights, sw, &R[i * nr]);
    vector<S> t2(nt);
    trig_row(rho_x, rho_y, ex, ey, &T[i * nt]);
    trig_row(rho_x, rho_y, -ey, ex, t2.data());
    for (int l = 0; l < nt; l++) T[i * nt + l] += t2[l];
    Wd[i] = double(W[i]);
    for (int l = 0; l < nr; l++) Rd[i * nr + l] = double(R[i * nr + l]);
    for (int l = 0; l < nt; l++) Td[i * nt + l] = double(T[i * nt + l]);
  });

  const auto elapsed = [&]() { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
  const double t_setup = elapsed();
  double t_apply = 0, t_gmres = 0;

  // Iterative refinement: residuals in S, corrections by GMRES in double
  JuliaResult<S> res;
  vector<S> h(N, S(0)), Lh(N), rs(N);
  vector<double> rd(N), dx;
  const double eps = is_same_v<S, double> ? 1e-15 : is_same_v<S, Expansion<2>> ? 1e-31 : 1e-46;
  for (;;) {
    const double ta = elapsed();
    apply(nr, nt, W, R, T, h, Lh);
    t_apply += elapsed() - ta;
    double rmax = 0;
    for (int64_t i = 0; i < N; i++) {
      rs[i] = (S(1) - h[i]) + Lh[i];
      rd[i] = double(rs[i]);
      rmax = max(rmax, abs(rd[i]));
    }
    // Stop at the precision's floor, or once a refinement stops helping (rounding in S limits the residual)
    const bool stalled = res.refinements && rmax > 0.1 * res.residual;
    res.residual = rmax;
    if (p.verbose) print("  refinement %d: residual %.3g", res.refinements, rmax);
    if (rmax <= eps || stalled || res.refinements >= p.max_refine) break;
    const double tg = elapsed();
    res.gmres_iters += gmres(nr, nt, Wd, Rd, Td, rd, dx, 1e-14);
    t_gmres += elapsed() - tg;
    for (int64_t i = 0; i < N; i++) h[i] += S(dx[i]);
    res.refinements++;
  }

  // ∫_X h, X = {ρ(θ) ≤ |z| ≤ r2} with |ρ^2 e^{2iθ} + c| = r2: trapezoid in θ, Gauss-Legendre in r
  const int nq = p.nq ? p.nq : 2 * nt + 1, ng = p.ng ? p.ng : nr + 16;
  vector<Arb> gx(ng), gw(ng);
  for (int g = 0; g < ng; g++) arb_hypgeom_legendre_p_ui_root(gx[g], gw[g], ng, g, prec);
  vector<S> part(nq);
  parallel_for(nq, [&](const int64_t q) {
    Arb sn, cs, bb, rho, rg, wg, t, u;
    // b = Re(conj(c) e^{2iθ}), θ = 2πq/nq;  ρ^2 = -b + √(b^2 - |c|^2 + r2^2)
    sin_cos_pi(sn, cs, 4 * q, nq);
    arb_mul(bb, cs, exact_arb(p.cx), prec);
    arb_mul(t, sn, exact_arb(p.cy), prec);
    arb_add(bb, bb, t, prec);
    arb_sqr(t, bb, prec);
    arb_set_d(u, r2); arb_sqr(u, u, prec);
    arb_add(t, t, u, prec);
    arb_set_d(u, p.cx); arb_sqr(u, u, prec); arb_sub(t, t, u, prec);
    arb_set_d(u, p.cy); arb_sqr(u, u, prec); arb_sub(t, t, u, prec);
    arb_sqrt(t, t, prec);
    arb_sub(rho, t, bb, prec);
    arb_sqrt(rho, rho, prec);
    // h along this angle: v[k] = Σ_j T(θ)[j] H[k,j], with η = e^{iθ/2} = e^{iπq/nq}
    sin_cos_pi(sn, cs, q, nq);
    vector<S> trow(nt), v(nr), crow(nr);
    trig_row(rho_x, rho_y, rnd<S>(cs), rnd<S>(sn), trow.data());
    for (int k = 0; k < nr; k++) {
      S dot(0);
      for (int j = 0; j < nt; j++) dot += trow[j] * h[int64_t(k) * nt + j];
      v[k] = dot;
    }
    S sum(0);
    Arb half_w, center;
    arb_set_d(t, r2);
    arb_sub(half_w, t, rho, prec); arb_mul_2exp_si(half_w, half_w, -1);
    arb_add(center, t, rho, prec); arb_mul_2exp_si(center, center, -1);
    for (int g = 0; g < ng; g++) {
      arb_mul(rg, half_w, gx[g], prec);
      arb_add(rg, rg, center, prec);
      arb_mul(wg, half_w, gw[g], prec);
      arb_mul(wg, wg, rg, prec);  // r dr
      arb_log(t, rg, prec);
      cheb_row(nodes, weights, rnd<S>(t), crow.data());
      S hv(0);
      for (int k = 0; k < nr; k++) hv += crow[k] * v[k];
      sum += rnd<S>(wg) * hv;
    }
    part[q] = sum;
  });
  S integral(0);
  for (const auto& s : part) integral += s;
  Arb pi_r2;
  arb_const_pi(pi_r2, prec);
  arb_mul_2exp_si(pi_r2, pi_r2, 1);
  arb_div_si(pi_r2, pi_r2, nq, prec);  // 2π / nq, the trapezoid weight
  const S tw = rnd<S>(pi_r2);
  arb_const_pi(pi_r2, prec);
  arb_set_d(a, r2);
  arb_sqr(a, a, prec);
  arb_mul(pi_r2, pi_r2, a, prec);
  res.area = rnd<S>(pi_r2) - tw * integral;
  res.secs = elapsed();
  if (p.verbose)
    print("  time: setup %.2f s, residuals %.2f s, GMRES %.2f s, area %.2f s", t_setup, t_apply, t_gmres,
          res.secs - t_setup - t_apply - t_gmres);
  return res;
}

template JuliaResult<double> julia_area(const JuliaParams&);
#define JULIA(S) \
  template JuliaResult<S> julia_area(const JuliaParams&); \
  template void cheb_row(const vector<S>&, const vector<S>&, S, S*); \
  template void trig_row(const vector<S>&, const vector<S>&, S, S, S*);
JULIA(Expansion<2>)
JULIA(Expansion<3>)

}  // namespace mandelbrot
