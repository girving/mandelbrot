// Areas of filled Julia sets K(c) from the area transfer operator

#include "julia.h"
#include "arb_cc.h"
#include "debug.h"
#include "expansion_arith.h"
#include "nearest.h"
#include "print.h"
#include <flint/acb.h>
#include <flint/arb_hypgeom.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <complex>
#include <functional>
#include <random>
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

// Angular grading: e^{iθ} = (e^{iφ} + a)/(1 + a e^{iφ}) with interpolation uniform in φ, so 0 ≤ a < 1 packs
// points near θ = 0 by (1 + a)/(1 - a).  mobius(out, e^{iφ}, a) gives e^{iθ}; mobius(out, e^{iθ}, -a) inverts it.
static void mobius(acb_t out, const acb_t e, const double a) {
  acb_t num, den;
  acb_init(num); acb_init(den);
  acb_set_d(num, a);
  acb_add(num, num, e, prec);
  acb_mul_arb(den, e, exact_arb(a), prec);
  acb_add_ui(den, den, 1, prec);
  acb_div(out, num, den, prec);
  acb_clear(num); acb_clear(den);
}

// e^{i arg(e)/2} as a pair, for a unit complex e
static void half_phase(arb_t ex, arb_t ey, const acb_t e) {
  Arb t;
  acb_arg(t, e, prec);
  arb_mul_2exp_si(t, t, -1);
  arb_sin_cos(ey, ex, t, prec);
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

// Interpolation grid: nr Chebyshev points (first kind) in s = log r on [s1, s2], with barycentric weights
// (-1)^k sin((2k+1)π/(2nr)), and nt angles θ_j = 2πj/nt through their half angle phases ρ_j = e^{iπj/nt}
template<class S> struct Grid {
  int nr, nt;
  double grade;  // Angular grading a: θ_j = θ(φ_j) with φ_j = 2πj/nt
  vector<Arb> s;  // The nodes in arb, for collocation points
  vector<S> nodes, weights, rho_x, rho_y;
  int64_t size() const { return int64_t(nr) * nt; }
};
template<class S> static Grid<S> make_grid(const arb_t s1, const arb_t s2, const int nr, const int nt,
                                           const double grade) {
  slow_assert(nt % 2 == 1 && nr >= 2, "need odd nt and nr ≥ 2");
  slow_assert(0 <= grade && grade < 1, "need 0 ≤ grade < 1");
  Grid<S> g;
  g.nr = nr; g.nt = nt; g.grade = grade;
  g.s.resize(nr);
  g.nodes.resize(nr); g.weights.resize(nr); g.rho_x.resize(nt); g.rho_y.resize(nt);
  Arb mid, half_len, a, b;
  arb_add(mid, s1, s2, prec); arb_mul_2exp_si(mid, mid, -1);
  arb_sub(half_len, s2, s1, prec); arb_mul_2exp_si(half_len, half_len, -1);
  for (int k = 0; k < nr; k++) {
    sin_cos_pi(a, b, 2 * k + 1, 2 * nr);
    arb_mul(g.s[k], half_len, b, prec);
    arb_add(g.s[k], g.s[k], mid, prec);
    if (k & 1) arb_neg(a, a);
    g.nodes[k] = rnd<S>(g.s[k]); g.weights[k] = rnd<S>(a);
  }
  for (int j = 0; j < nt; j++) {
    sin_cos_pi(a, b, j, nt);
    g.rho_x[j] = rnd<S>(b); g.rho_y[j] = rnd<S>(a);
  }
  return g;
}

// Preimage data of a collocation point z = z0 + e^{s_k} e^{iθ_j} (z0 the annulus center): W = 1/(4|z - c|) =
// 1/|f'(w)|^2 and, for each preimage w = ±√(z - c), its log radius log|w - z0| and the half phase e^{iφ/2} of its
// interpolation angle φ.  With z0 = 0 both preimages share a radius.
template<class S> struct Pre { S W, sw[2], ex[2], ey[2]; };
template<class S> static vector<Pre<S>> preimages(const Grid<S>& g, const double cx, const double cy,
                                                   const acb_t z0) {
  vector<Pre<S>> pre(g.size());
  parallel_for(g.size(), [&](const int64_t i) {
    const int k = int(i / g.nt), j = int(i % g.nt);
    Arb r, sn, cs, m, phi, t;
    acb_t u, e, w, d;
    acb_init(u); acb_init(e); acb_init(w); acb_init(d);
    arb_exp(r, g.s[k], prec);
    sin_cos_pi(acb_imagref(e), acb_realref(e), 2 * j, g.nt);
    mobius(u, e, g.grade);
    acb_mul_arb(u, u, r, prec);
    acb_add(u, u, z0, prec);
    arb_sub(acb_realref(u), acb_realref(u), exact_arb(cx), prec);
    arb_sub(acb_imagref(u), acb_imagref(u), exact_arb(cy), prec);
    acb_abs(m, u, prec);
    acb_arg(phi, u, prec);
    auto& q = pre[i];
    arb_mul_2exp_si(t, m, 2);
    arb_inv(t, t, prec);
    q.W = rnd<S>(t);
    // The preimages ±w, w = √|z - c| e^{i arg(z - c)/2}, relative to the center; map each angle back to φ
    arb_sqrt(r, m, prec);
    arb_mul_2exp_si(phi, phi, -1);
    arb_sin_cos(acb_imagref(w), acb_realref(w), phi, prec);
    acb_mul_arb(w, w, r, prec);
    for (int b = 0; b < 2; b++) {
      if (b) acb_neg(w, w);
      acb_sub(d, w, z0, prec);
      acb_abs(t, d, prec);
      acb_div_arb(e, d, t, prec);
      arb_log(t, t, prec);
      q.sw[b] = rnd<S>(t);
      mobius(u, e, -g.grade);
      half_phase(cs, sn, u);
      q.ex[b] = rnd<S>(cs); q.ey[b] = rnd<S>(sn);
    }
    acb_clear(u); acb_clear(e); acb_clear(w); acb_clear(d);
  });
  return pre;
}

// L restricted to some points, acting on functions interpolated from a grid, with nb rows (R, T) per point:
//   (L h)_i = W_i Σ_b Σ_k R[i nb + b, k] Σ_j T[i nb + b, j] H[k, j]
// nb = 1 when both preimages share a radius (annulus centered at 0), summing their angular rows; else nb = 2.
template<class S> struct Factors {
  int nr, nt, nb;
  vector<S> W, R, T;
  int64_t rows() const { return int64_t(W.size()); }
};
template<class S> static Factors<S> factors(const vector<Pre<S>>& pre, const Grid<S>& g, const int nb) {
  Factors<S> F;
  const int nr = F.nr = g.nr, nt = F.nt = g.nt;
  F.nb = nb;
  const int64_t n = int64_t(pre.size());
  F.W.resize(n); F.R.resize(n * nb * nr); F.T.resize(n * nb * nt);
  parallel_for(n, [&](const int64_t i) {
    const auto& q = pre[i];
    F.W[i] = q.W;
    if (nb == 1) {
      cheb_row(g.nodes, g.weights, q.sw[0], &F.R[i * nr]);
      vector<S> t2(nt);
      trig_row(g.rho_x, g.rho_y, q.ex[0], q.ey[0], &F.T[i * nt]);
      trig_row(g.rho_x, g.rho_y, q.ex[1], q.ey[1], t2.data());
      for (int l = 0; l < nt; l++) F.T[i * nt + l] += t2[l];
    } else {
      for (int b = 0; b < 2; b++) {
        cheb_row(g.nodes, g.weights, q.sw[b], &F.R[(2 * i + b) * nr]);
        trig_row(g.rho_x, g.rho_y, q.ex[b], q.ey[b], &F.T[(2 * i + b) * nt]);
      }
    }
  });
  return F;
}
template<class S> static Factors<double> to_double(const Factors<S>& F) {
  Factors<double> D;
  D.nr = F.nr; D.nt = F.nt; D.nb = F.nb;
  const auto conv = [](const vector<S>& x, vector<double>& y) { y.resize(x.size()); for (size_t i = 0; i < x.size(); i++) y[i] = double(x[i]); };
  conv(F.W, D.W); conv(F.R, D.R); conv(F.T, D.T);
  return D;
}
template<class S> static vector<Pre<double>> to_double(const vector<Pre<S>>& pre) {
  vector<Pre<double>> d(pre.size());
  for (size_t i = 0; i < pre.size(); i++)
    d[i] = {double(pre[i].W), {double(pre[i].sw[0]), double(pre[i].sw[1])},
            {double(pre[i].ex[0]), double(pre[i].ex[1])}, {double(pre[i].ey[0]), double(pre[i].ey[1])}};
  return d;
}

template<class S> static void apply(const Factors<S>& F, const vector<S>& h, vector<S>& y) {
  const int nr = F.nr, nt = F.nt, nb = F.nb;
  y.resize(F.rows());
  parallel_for(F.rows(), [&](const int64_t i) {
    S sum(0);
    for (int b = 0; b < nb; b++) {
      const S* r = &F.R[(i * nb + b) * nr];
      const S* t = &F.T[(i * nb + b) * nt];
      for (int k = 0; k < nr; k++) {
        const S* hk = &h[int64_t(k) * nt];
        S dot(0);
        for (int j = 0; j < nt; j++) dot += t[j] * hk[j];
        sum += r[k] * dot;
      }
    }
    y[i] = F.W[i] * sum;
  });
}

typedef function<void(const vector<double>&, vector<double>&)> Op;

// Double: the generic loop's inner dot product is a reduction, which strict floating point keeps scalar.  Instead
// transpose H and accumulate axpys acc[k] += T[i,j] H[k,j] over j, which vectorize, four rows at a time so each
// load of H feeds four rows.
template<> void apply(const Factors<double>& F, const vector<double>& h, vector<double>& y) {
  const int nr = F.nr, nt = F.nt, nb = F.nb, B = 4;
  const int64_t n = F.rows() * nb;  // Rows (point, branch)
  vector<double> v(n);
  vector<double> ht(int64_t(nt) * nr);
  for (int k = 0; k < nr; k++)
    for (int j = 0; j < nt; j++) ht[int64_t(j) * nr + k] = h[int64_t(k) * nt + j];
  parallel_for((n + B - 1) / B, [&](const int64_t blk) {
    const int64_t i0 = blk * B, rows = std::min(int64_t(B), n - i0);
    double acc[B][512];
    slow_assert(nr <= 512);
    for (int b = 0; b < B; b++) for (int k = 0; k < nr; k++) acc[b][k] = 0;
    const double* t[B];
    for (int b = 0; b < B; b++) t[b] = &F.T[(i0 + std::min(int64_t(b), rows - 1)) * nt];
    for (int j = 0; j < nt; j++) {
      const double* hj = &ht[int64_t(j) * nr];
      const double t0 = t[0][j], t1 = t[1][j], t2 = t[2][j], t3 = t[3][j];
      for (int k = 0; k < nr; k++) {
        const double v = hj[k];
        acc[0][k] = fma(t0, v, acc[0][k]);
        acc[1][k] = fma(t1, v, acc[1][k]);
        acc[2][k] = fma(t2, v, acc[2][k]);
        acc[3][k] = fma(t3, v, acc[3][k]);
      }
    }
    for (int b = 0; b < rows; b++) {
      const double* r = &F.R[(i0 + b) * nr];
      double sum = 0;
      for (int k = 0; k < nr; k++) sum = fma(r[k], acc[b][k], sum);
      v[i0 + b] = sum;
    }
  });
  y.resize(F.rows());
  for (int64_t i = 0; i < F.rows(); i++) y[i] = F.W[i] * (nb == 1 ? v[i] : v[2 * i] + v[2 * i + 1]);
}

// Right preconditioned restarted GMRES for A x = b: solves A M u = b, x = M u (M may be null).  Stops at tol, or when
// a restart gains less than a factor of 2 (rounding floors vary with the grid; refinement absorbs them).  Returns
// iterations.
static int gmres(const int64_t N, const Op& A, const Op& M, const vector<double>& b, vector<double>& x,
                 const double tol, const bool verbose = false) {
  const int m = 60;
  const auto norm = [](const vector<double>& v) { double s = 0; for (const double a : v) s += a * a; return sqrt(s); };
  const double bnorm = norm(b);
  x.assign(N, 0);
  if (bnorm == 0) return 0;
  vector<double> r(N), w(N), z(N), mz(N);
  int iters = 0;
  double last = INFINITY;
  for (int restart = 0; restart < 100; restart++) {
    A(x, w);
    for (int64_t i = 0; i < N; i++) r[i] = b[i] - w[i];
    const double beta = norm(r);
    if (verbose) print("    gmres restart %d: iterations %d, residual %.3g", restart, iters, beta / bnorm);
    if (beta <= tol * bnorm || beta > 0.5 * last) return iters;
    last = beta;
    vector<vector<double>> V(1, r);
    for (auto& a : V[0]) a /= beta;
    vector<vector<double>> H(m + 1, vector<double>(m, 0));
    vector<double> cs(m), sn(m), g(m + 1, 0);
    g[0] = beta;
    int k = 0;
    for (; k < m; k++) {
      iters++;
      if (M) { M(V[k], z); A(z, w); }
      else A(V[k], w);
      for (int pass = 0; pass < 2; pass++)  // Modified Gram-Schmidt, twice
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
    z.assign(N, 0);
    for (int l = 0; l < k; l++)
      for (int64_t i = 0; i < N; i++) z[i] += y[l] * V[l][i];
    if (M) M(z, mz);
    else mz = z;
    for (int64_t i = 0; i < N; i++) x[i] += mz[i];
  }
  die("gmres did not converge");
}

// Fast evaluation of grid functions' spectral interpolants at many points, NUFFT style.  The interpolant is a
// cosine series in α = arccos x (x ∈ [-1, 1] the Chebyshev variable: T_n = cos nα, so even and 2π periodic in α)
// times a Fourier series in φ.  Put it on a grid oversampled by os in both, by a tensor product Ps H Ptᵀ costing
// N (os nr + os^2 nt), then interpolate locally with separable sw-point stencils.  Two kernels:
//   es: the exponential of semicircle ψ(x) = exp(β (√(1 - (2x/(sw h))^2) - 1)), with the grid values deconvolved
//       by 1/ψ̂ (Poisson summation: Σ_u ψ(α - α_u) e^{inα_u} ≈ ψ̂(n) e^{inα} / h), accurate to ~1e-13 at os = 2,
//       sw = 13;
//   Lagrange interpolation of the plain oversampled values, cheaper but only ~1e-3 at os = 2, sw = 6, which is
//       plenty for the preconditioner.
struct Oversampled {
  int nr, nt, os, sw, na, nf;
  bool es;
  double sa, shalf, grade, beta;
  std::complex<double> z0;  // Annulus center
  vector<double> Ps, PtT;  // na x nr and nt x nf
  struct Stencil { int32_t u[16], v[16]; double wu[16], wv[16]; };

  // The kernel in units of the grid spacing, and its Fourier transform at frequency ω (radians per spacing)
  double psi(const double x) const {
    const double t = 2 * x / sw;
    return std::abs(t) < 1 ? std::exp(beta * (std::sqrt(1 - t * t) - 1)) : 0;
  }
  double psi_hat(const double omega) const {
    double sum = 0;  // ∫ ψ(x) cos(ωx) dx over [-sw/2, sw/2], midpoint rule (ψ is smooth and vanishes to all orders)
    const int n = 2000;
    for (int i = 0; i < n; i++) {
      const double x = (i + 0.5) / n * sw - sw / 2.0;
      sum += psi(x) * std::cos(omega * x);
    }
    return sum * sw / n;
  }

  Oversampled(const Grid<double>& g, const double sa, const double sb, const double grade, const int os, const int sw,
              const bool es, const std::complex<double> z0)
      : nr(g.nr), nt(g.nt), os(os), sw(sw), na(os * g.nr), nf(os * g.nt), es(es), sa(sa), shalf((sb - sa) / 2),
        grade(grade), beta(2.30 * sw), z0(z0) {
    slow_assert(2 <= sw && sw <= 16 && sw <= na && os >= 1, "bad stencil width %d or oversampling %d", sw, os);
    Ps.resize(int64_t(na) * nr); PtT.resize(int64_t(nt) * nf);
    if (!es) {
      vector<double> row(nt);
      for (int u = 0; u < na; u++)
        cheb_row(g.nodes, g.weights, sa + shalf * (1 + std::cos((u + 0.5) * M_PI / na)), &Ps[int64_t(u) * nr]);
      for (int v = 0; v < nf; v++) {
        trig_row(g.rho_x, g.rho_y, std::cos(M_PI * v / nf), std::sin(M_PI * v / nf), row.data());
        for (int j = 0; j < nt; j++) PtT[int64_t(j) * nf + v] = row[j];
      }
      return;
    }
    // α: nodal values at α_k = (2k+1)π/(2nr) → coefficients a_n = (2 - δ_n0)/nr Σ_k f_k cos(n α_k) → grid values
    // g_u = Σ_n a_n cos(n α_u) / ψ̂(n h), with α_u = (u + 1/2) h, h = π/na, the kernel in units of h
    vector<double> da(nr), df(nt / 2 + 1);
    for (int n = 0; n < nr; n++) da[n] = 1 / psi_hat(n * M_PI / na);
    for (int m = 0; m <= nt / 2; m++) df[m] = 1 / psi_hat(m * 2 * M_PI / nf);
    for (int u = 0; u < na; u++)
      for (int k = 0; k < nr; k++) {
        double sum = 0;
        for (int n = 0; n < nr; n++)
          sum += (n ? 2.0 : 1.0) / nr * std::cos(n * (2 * k + 1) * M_PI / (2 * nr)) * da[n] * std::cos(n * (u + 0.5) * M_PI / na);
        Ps[int64_t(u) * nr + k] = sum;
      }
    // φ: nodal values at φ_j = 2πj/nt → c_m = Σ_j f_j e^{-imφ_j} / nt → g_v = Σ_|m|≤nt/2 c_m e^{imφ_v} / ψ̂(m h_φ)
    for (int v = 0; v < nf; v++)
      for (int j = 0; j < nt; j++) {
        double sum = df[0];
        for (int m = 1; m <= nt / 2; m++) sum += 2 * df[m] * std::cos(m * (2 * M_PI * v / nf - 2 * M_PI * j / nt));
        PtT[int64_t(j) * nf + v] = sum / nt;
      }
  }

  // Stencil for weight times the interpolant at z
  void stencil(const std::complex<double> zabs, const double weight, Stencil& st) const {
    const std::complex<double> z = zabs - z0;
    // α on the grid α_u = (u + 1/2) π / na, reflected evenly across 0 and π
    const double x = std::clamp((std::log(std::abs(z)) - sa) / shalf - 1, -1.0, 1.0);
    const double ua = std::acos(x) / M_PI * na - 0.5;
    const int u0 = int(std::floor(ua)) - (sw - 1) / 2;
    // φ on the grid φ_v = 2π v / nf, from e^{iθ} by the inverse grading map
    const std::complex<double> a = grade, e = z / std::abs(z), ef = (e - a) / (1.0 - a * e);
    double vf = std::arg(ef) / (2 * M_PI) * nf;
    vf -= nf * std::floor(vf / nf);
    const int v0 = int(std::floor(vf)) - (sw - 1) / 2;
    for (int q = 0; q < sw; q++) {
      if (es) {
        st.wu[q] = weight * psi(ua - u0 - q);
        st.wv[q] = psi(vf - v0 - q);
      } else {
        st.wu[q] = weight;
        st.wv[q] = 1;
        for (int r = 0; r < sw; r++) if (r != q) {
          st.wu[q] *= (ua - u0 - r) / (q - r);
          st.wv[q] *= (vf - v0 - r) / (q - r);
        }
      }
      int uu = ((u0 + q) % (2 * na) + 2 * na) % (2 * na);  // Period 2 na in u, even about u = -1/2 and na - 1/2
      if (uu >= na) uu = 2 * na - 1 - uu;
      st.u[q] = uu;
      st.v[q] = ((v0 + q) % nf + nf) % nf;
    }
  }

  // G = Ps H Ptᵀ on rows [ulo, uhi] of the oversampled grid, H = h as nr x nt, as vectorizable axpys
  void grid(const vector<double>& h, vector<double>& G, const int ulo, const int uhi) const {
    G.resize(int64_t(na) * nf);
    parallel_for(uhi - ulo + 1, [&](const int64_t du) {
      const int64_t u = ulo + du;
      vector<double> g(nt, 0.0);
      for (int k = 0; k < nr; k++) {
        const double c = Ps[u * nr + k];
        const double* hk = &h[int64_t(k) * nt];
        for (int j = 0; j < nt; j++) g[j] = fma(c, hk[j], g[j]);
      }
      double* Gu = &G[u * nf];
      for (int v = 0; v < nf; v++) Gu[v] = 0;
      for (int j = 0; j < nt; j++) {
        const double c = g[j];
        const double* pj = &PtT[int64_t(j) * nf];
        for (int v = 0; v < nf; v++) Gu[v] = fma(c, pj[v], Gu[v]);
      }
    });
  }

  double eval(const vector<double>& G, const Stencil& st) const {
    double sum = 0;
    for (int q = 0; q < sw; q++) {
      const double* Gu = &G[int64_t(st.u[q]) * nf];
      double t = 0;
      for (int r = 0; r < sw; r++) t = fma(st.wv[r], Gu[st.v[r]], t);
      sum = fma(st.wu[q], t, sum);
    }
    return sum;
  }
};

template<class S> JuliaResult<S> julia_area(const JuliaParams& p) {
  const auto t0 = std::chrono::steady_clock::now();
  const auto elapsed = [&]() { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
  const int nr = p.nr, nt = p.nt;
  const int64_t N = int64_t(nr) * nt;
  // The annulus A = {r1 < |z - z0| < r2}.  Centered at the attracting fixed point α = (1 - √(1 - 4c))/2, with
  // multiplier λ = 2α, f(z) - α = (z - α)(z + α), so |f(z) - α| ≤ |z - α| (|z - α| + |λ|): the disk |z - α| ≤ r1
  // maps strictly into itself if r1 < 1 - |λ|, and |z - α| ≥ r2 grows to escape if r2 > 1 + |λ|.  This covers
  // the whole main cardioid.  Centered at 0 (|c| < 1/4): r1 - |c| > r1^2 and r2^2 - r2 > |c|, and both preimages
  // of a point share a radius, halving the operator's rows.
  //
  // Either way V must contain the critical value c, and so the whole postcritical orbit: L's weight 1/(4|z - c|)
  // makes h singular like 1/|z - p| at each postcritical point p in A.  Around α that needs |c - α| =
  // |α||1 - α| < r1 < 1 - |λ|, which on the real axis reaches only c ∈ (-0.394, 0.236): beyond, near the
  // cardioid's satellite roots, the forward-invariant region containing the postcritical orbit is far from round.
  const std::complex<double> c(p.cx, p.cy), alpha = (1.0 - std::sqrt(1.0 - 4.0 * c)) / 2.0, lambda = 2.0 * alpha;
  const double ac = std::abs(c), al = std::abs(lambda), cd = std::abs(c - alpha);
  const bool centered = p.center == 1 || (p.center < 0 && ac >= 0.25);  // Auto: the origin when it works
  slow_assert(!centered || (al < 1 && cd < 1 - al), "c = %g + %gi has no round interior disk at 0 or at α "
              "(|c - α| = %g, 1 - |λ| = %g)", p.cx, p.cy, cd, 1 - al);
  const double r1 = p.r1 ? p.r1 : centered ? (cd + 1 - al) / 2 : max(2 * ac, 0.25);
  const double r2 = p.r2 ? p.r2 : centered ? 1.5 + al : 1.5;
  if (centered)
    slow_assert(cd < r1 && r1 + al < 1 && r2 > 1 + al, "annulus %g < |z - α| < %g does not work for |λ| = %g, "
                "|c - α| = %g", r1, r2, al, cd);
  else {
    slow_assert(r1 - ac > r1 * r1 && r2 * r2 - r2 > ac, "annulus %g < |z| < %g does not work for |c| = %g", r1, r2, ac);
    slow_assert(sqrt(r1 - ac) > r1 && sqrt(r2 + ac) < r2);  // Preimage radii [√(r1 - |c|), √(r2 + |c|)] inside
  }
  const std::complex<double> z0 = centered ? alpha : 0.0;
  acb_t z0a;  // The center in arb
  acb_init(z0a);
  if (centered) {
    acb_t t;
    acb_init(t);
    acb_set_d_d(t, p.cx, p.cy);
    acb_mul_2exp_si(t, t, 2);
    acb_neg(t, t);
    acb_add_ui(t, t, 1, prec);
    acb_sqrt(t, t, prec);
    acb_neg(t, t);
    acb_add_ui(t, t, 1, prec);
    acb_mul_2exp_si(z0a, t, -1);
    acb_clear(t);
  }
  Arb s1, s2, a;
  arb_set_d(a, r1); arb_log(s1, a, prec);
  arb_set_d(a, r2); arb_log(s2, a, prec);

  // The operator on the fine grid, in S and double
  const auto fine = make_grid<S>(s1, s2, nr, nt, p.grade);
  const auto& nodes = fine.nodes;
  const auto& weights = fine.weights;
  const auto& rho_x = fine.rho_x;
  const auto& rho_y = fine.rho_y;
  const auto pre = preimages(fine, p.cx, p.cy, z0a);
  const auto F = factors(pre, fine, centered ? 2 : 1);
  const auto Fd = to_double(F);

  const Op A = [&](const vector<double>& v, vector<double>& out) {
    apply(Fd, v, out);
    for (int64_t i = 0; i < N; i++) out[i] = v[i] - out[i];
  };
  // Preconditioner for the slow modes as α nears the boundary of the cardioid (|λ| → 1).  Let w+ = σ√(z - c) be
  // the inverse branch fixing α (σ = ±1).  Near a parabolic parameter the cycle born from α (the fixed point q near
  // the cusp, the 2-cycle near -3/4, ...) is attracting for w+ with rate ~|λ|, and w+'s iterates converge to it
  // from the whole annulus.  Split L = L+ + L- by branch: L+ carries a cluster of eigenvalues → 1 (dilations in the
  // cycle's Koenigs coordinate), while L- adds little to the leading eigenvalue (O(ε^2) at the cusp), so
  // L- (1 - L+)^-1 stays well conditioned.  And L+^m g(z) = |(w+^m)'(z)|^2 g(w+^m(z)) has L's factored form for
  // any m: one interpolation row at the point w+^m(z).  So
  //   M = Π_{j<J} (1 + L+^{2^j}) = Σ_{n<2^J} L+^n ≈ (1 - L+)^-1
  // costs J single-row matvecs, with 2^J steps enough that |λ|^{2^{J+1}} ≤ 1/100.
  const double sigma = std::abs(std::sqrt(alpha - c) - alpha) < std::abs(std::sqrt(alpha - c) + alpha) ? 1 : -1;
  const double mu = 1 / al;
  const int J = p.pre == 0 || al == 0 ? 0 : max(0, int(std::ceil(std::log2(std::log(100.0) / (2 * std::log(mu))))));
  const bool use_pre = p.pre == 1 || (p.pre < 0 && J >= 4);
  // Each level L+^{2^l} evaluates a grid function's spectral interpolant at the points w+^{2^l}(z_i).  Full spectral
  // rows cost a matvec per level, and local stencils on the collocation grid are not enough: near q, L+ carries even
  // grid-scale oscillations with weight ≈ μ^-2 ≈ 1, so the slow cluster includes Nyquist-scale modes, which local
  // interpolation damps.  Oversampled interpolation handles them, and each level computes only the oversampled rows
  // its points touch, a narrow band at deep levels where every point is near q.
  const auto fine_d = make_grid<double>(s1, s2, nr, nt, p.grade);
  const Oversampled over(fine_d, rnd<double>(s1), rnd<double>(s2), p.grade, p.oversample, p.stencil, false, z0);
  const Oversampled over_fast(fine_d, rnd<double>(s1), rnd<double>(s2), p.grade, 2, p.fast_width, true, z0);
  vector<std::complex<double>> zs(N);  // Collocation points in double
  for (int64_t i = 0; i < N; i++) {
    const std::complex<double> e = std::polar(1.0, 2 * M_PI * (i % nt) / nt), a = p.grade;
    zs[i] = z0 + std::exp(fine_d.nodes[i / nt]) * (e + a) / (1.0 + a * e);
  }
  vector<vector<Oversampled::Stencil>> levels(use_pre ? J : 0, vector<Oversampled::Stencil>(N));
  vector<int> ulo(J, over.na), uhi(J, -1);
  if (use_pre) {
    parallel_for(N, [&](const int64_t i) {
      std::complex<double> z = zs[i];
      double jac = 1;
      int64_t m = 0;
      for (int l = 0; l < J; l++) {
        for (; m < (int64_t(1) << l); m++) {  // Advance to w+^{2^l}(z), accumulating |w+'|^2 = 1/(4|z - c|)
          jac /= 4 * std::abs(z - c);
          z = sigma * std::sqrt(z - c);
        }
        over.stencil(z, jac, levels[l][i]);
      }
    });
    for (int l = 0; l < J; l++)
      for (const auto& st : levels[l])
        for (int q = 0; q < over.sw; q++) {
          ulo[l] = std::min(ulo[l], int(st.u[q]));
          uhi[l] = std::max(uhi[l], int(st.u[q]));
        }
  }
  const Op M = !use_pre ? Op() : Op([&](const vector<double>& r, vector<double>& out) {
    out = r;
    vector<double> G, t(N);
    for (int l = 0; l < J; l++) {
      over.grid(out, G, ulo[l], uhi[l]);
      parallel_for(N, [&](const int64_t i) { t[i] = over.eval(G, levels[l][i]); });
      for (int64_t i = 0; i < N; i++) out[i] += t[i];
    }
  });

  // Fast approximate L for the double precision corrections: both preimages of each collocation point through the
  // oversampled grid, O(N (nr + nt)) rather than O(N^2).  Refinement residuals still use the exact L in S.
  vector<Oversampled::Stencil> fast(p.fast ? 2 * N : 0);
  if (p.fast)
    parallel_for(N, [&](const int64_t i) {
      const std::complex<double> w = std::sqrt(zs[i] - c);
      const double W = 1 / (4 * std::abs(zs[i] - c));
      over_fast.stencil(w, W, fast[2 * i]);
      over_fast.stencil(-w, W, fast[2 * i + 1]);
    });
  const Op Afast = [&](const vector<double>& v, vector<double>& out) {
    vector<double> G;
    over_fast.grid(v, G, 0, over_fast.na - 1);
    out.resize(N);
    parallel_for(N, [&](const int64_t i) {
      out[i] = v[i] - (over_fast.eval(G, fast[2 * i]) + over_fast.eval(G, fast[2 * i + 1]));
    });
  };
  if (p.fast && p.verbose) {  // Accuracy of the fast L on a random vector
    std::mt19937_64 rng(5);
    std::uniform_real_distribution<double> u(-1, 1);
    vector<double> v(N), e, f;
    for (auto& x : v) x = u(rng);
    A(v, e);
    Afast(v, f);
    double d = 0, m = 0;
    for (int64_t i = 0; i < N; i++) { d = max(d, abs(e[i] - f[i])); m = max(m, abs(v[i])); }
    print("  fast L: max |(L - L_fast) v| / max |v| = %.3g on random v", d / m);
  }
  if (p.verbose)
    print("  annulus %.6g < |z - (%.6g + %.6gi)| < %.6g, |λ| %.6f%s", r1, z0.real(), z0.imag(), r2, al,
          use_pre ? tfm::format(", parabolic preconditioner: branch %+g, J = %d", sigma, J) : "");
  const double t_setup = elapsed();
  double t_apply = 0, t_gmres = 0;

  JuliaResult<S> res;
  // L is positive, so its leading eigenvalue is real and positive: power iteration from 1, normalized in max norm
  if (p.eig) {
    vector<double> v(N, 1.0), w(N);
    double last = 0;
    for (int it = 0; it < 100000; it++) {
      apply(Fd, v, w);
      double m = 0;
      for (const double a : w) m = max(m, abs(a));
      for (int64_t i = 0; i < N; i++) v[i] = w[i] / m;
      res.rho = m;
      if (p.verbose && it % 1000 == 0) print("  power iteration %d: %.12f", it, m);
      if (it > 10 && abs(m - last) <= 1e-13 * m) break;
      last = m;
    }
  }

  // Iterative refinement: residuals in S, corrections by GMRES in double
  vector<S> h(N, S(0)), Lh(N), rs(N);
  vector<double> rd(N), dx;
  const double eps = is_same_v<S, double> ? 1e-15 : is_same_v<S, Expansion<2>> ? 1e-31 : 1e-46;
  for (;;) {
    const double ta = elapsed();
    apply(F, h, Lh);
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
    if (p.verbose) print("  refinement %d: residual %.3g%s", res.refinements, rmax, p.eig ? tfm::format(" (rho %.12f)", res.rho) : "");
    if (rmax <= eps || stalled || res.refinements >= p.max_refine) break;
    const double tg = elapsed();
    res.gmres_iters += gmres(N, p.fast ? Afast : A, M, rd, dx, 1e-14, p.verbose);
    t_gmres += elapsed() - tg;
    for (int64_t i = 0; i < N; i++) h[i] += S(dx[i]);
    res.refinements++;
  }

  // ∫_X h, X = {ρ(θ) ≤ |z - z0| ≤ r2} = {z ∈ A : |f(z) - z0| ≥ r2}: trapezoid in φ (weighted by dθ/dφ),
  // Gauss-Legendre in r.  Centered at 0, |ρ^2 e^{2iθ} + c| = r2; at α, ρ |ρ e^{iθ} + λ| = r2.
  const int nq = p.nq ? p.nq : 2 * nt + 1, ng = p.ng ? p.ng : nr + 16;
  vector<Arb> gx(ng), gw(ng);
  for (int g = 0; g < ng; g++) arb_hypgeom_legendre_p_ui_root(gx[g], gw[g], ng, g, prec);
  vector<S> part(nq);
  parallel_for(nq, [&](const int64_t q) {
    Arb sn, cs, bb, rho, rg, wg, t, u, jac;
    // θ = θ(φ) for φ = 2πq/nq, dθ/dφ = (1 - a^2)/|1 + a e^{iφ}|^2, b = Re(conj(c) e^{2iθ}), ρ^2 = -b + √(b^2 - |c|^2 + r2^2)
    acb_t e, et;
    acb_init(e); acb_init(et);
    sin_cos_pi(acb_imagref(e), acb_realref(e), 2 * q, nq);
    mobius(et, e, p.grade);
    acb_mul_arb(e, e, exact_arb(p.grade), prec);
    acb_add_ui(e, e, 1, prec);
    acb_abs(jac, e, prec);
    arb_sqr(jac, jac, prec);
    arb_set_d(u, p.grade); arb_sqr(u, u, prec); arb_set_ui(t, 1); arb_sub(t, t, u, prec);
    arb_div(jac, t, jac, prec);
    if (!centered) {
      acb_sqr(e, et, prec);
      arb_mul(bb, acb_realref(e), exact_arb(p.cx), prec);
      arb_mul(t, acb_imagref(e), exact_arb(p.cy), prec);
      arb_add(bb, bb, t, prec);
      arb_sqr(t, bb, prec);
      arb_set_d(u, r2); arb_sqr(u, u, prec);
      arb_add(t, t, u, prec);
      arb_set_d(u, p.cx); arb_sqr(u, u, prec); arb_sub(t, t, u, prec);
      arb_set_d(u, p.cy); arb_sqr(u, u, prec); arb_sub(t, t, u, prec);
      arb_sqrt(t, t, prec);
      arb_sub(rho, t, bb, prec);
      arb_sqrt(rho, rho, prec);
    } else {
      // g(ρ) = ρ^2 (ρ^2 + 2 b ρ + |λ|^2) - r2^2 with b = Re(λ e^{-iθ}), increasing on the root's bracket: a double
      // bisection, then Newton in arb
      acb_t lam;
      acb_init(lam);
      acb_mul_2exp_si(lam, z0a, 1);
      acb_conj(e, et);
      acb_mul(e, e, lam, prec);
      arb_set(bb, acb_realref(e));
      Arb l2;
      acb_abs(l2, lam, prec);
      arb_sqr(l2, l2, prec);
      acb_clear(lam);
      const double bd = rnd<double>(bb), l2d = rnd<double>(l2);
      const auto g = [&](const double x) { return x * x * (x * x + 2 * bd * x + l2d) - r2 * r2; };
      double lo = r1, hi = r2;
      slow_assert(g(lo) < 0 && g(hi) > 0);
      for (int it = 0; it < 60; it++) { const double md = (lo + hi) / 2; (g(md) < 0 ? lo : hi) = md; }
      arb_set_d(rho, (lo + hi) / 2);
      Arb gv, dg, x2, r22;
      arb_set_d(r22, r2); arb_sqr(r22, r22, prec);
      for (int it = 0; it < 5; it++) {  // g = x^4 + 2b x^3 + l2 x^2 - r2^2, g' = 4x^3 + 6b x^2 + 2 l2 x
        arb_sqr(x2, rho, prec);
        arb_mul(gv, bb, rho, prec); arb_mul_2exp_si(gv, gv, 1); arb_add(gv, gv, x2, prec); arb_add(gv, gv, l2, prec);
        arb_mul(gv, gv, x2, prec); arb_sub(gv, gv, r22, prec);
        arb_mul_2exp_si(dg, x2, 2); arb_mul(t, bb, rho, prec); arb_mul_ui(t, t, 6, prec); arb_add(dg, dg, t, prec);
        arb_mul_2exp_si(t, l2, 1); arb_add(dg, dg, t, prec); arb_mul(dg, dg, rho, prec);
        arb_div(gv, gv, dg, prec);
        arb_sub(rho, rho, gv, prec);
        mag_zero(arb_radref(rho.x));  // Newton's iterate is just a point; keep the ball from growing
      }
    }
    acb_clear(e); acb_clear(et);
    // h along this angle: v[k] = Σ_j T(φ)[j] H[k,j], with η = e^{iφ/2} = e^{iπq/nq}
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
    part[q] = rnd<S>(jac) * sum;
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
  acb_clear(z0a);
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
