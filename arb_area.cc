// Mandelbrot area via Arb

#include "debug.h"
#include "known.h"
#include "poly.h"
#include "print.h"
#include "wall_time.h"
namespace mandelbrot {

using std::min;
using std::max;
using std::runtime_error;
using std::span;

// C = a*A - 2^(log2_b)*B
void poly_linear_sub_si_2exp(Poly& C, slong a, const Poly& A, const slong log2_b, const Poly& B, slong prec) {
  const auto n = max(A.length(), B.length());
  arb_poly_zero(C);
  arb_poly_fit_length(C, n);
  _arb_poly_set_length(C, n);
  for (int i = 0; i < n; i++) {
    auto c = C.mutable_at(i);
    arb_mul_2exp_si(c, B[i], log2_b);
    arb_neg(c, c);
    arb_addmul_si(c, A[i], a, prec);
  }
}

// h = log escape(k, z*e^-g)^(2^-k) in mandelbrot-area-cupy.ipynb notation
void escape(Poly& h, Poly& dh, const int k, const Poly& g, const Poly& dg,
            const int n, const int dn, const int prec) {
  const bool verbose = false;

  // Base case
  arb_poly_zero(h);
  arb_poly_zero(dh);

  Poly t, dt, s, ds, u;
  for (int i = 1; i <= k; i++) {
    const int p = 1 << i;

    // t = (1-p)g - p h
    poly_linear_sub_si_2exp(t, 1-p, g, i, h, prec);
    poly_linear_sub_si_2exp(dt, 1-p, dg, i, dh, prec);
    if (verbose) stats(t, "linear %d", i);

    // h += (1/p) log(1 + z^(p-1)exp(t))
    poly_log1p_exp_shift_refine(s, t, p-1, n-(p-1), prec);  // s = z^(1-p) log(1 + z^(p-1) exp(t))
    poly_sub_shift_series(u, t, s, p-1, dn-(p-1), prec);    // u = t - z^(p-1) s
    poly_exp_refine(t, u, dn-(p-1), prec);                  // t = exp(t - z^(p-1) s)
    safe_poly_mullow(ds, t, dt, dn-(p-1), prec);            // ds = exp(t - z^(p-1) s) dt
    arb_poly_scalar_mul_2exp_si(s, s, -i);                  // s /= p
    arb_poly_scalar_mul_2exp_si(ds, ds, -i);                // ds /= p
    poly_add_shift_series(h, h, s, p-1, n, prec);           // h += z^(p-1) s
    poly_add_shift_series(dh, dh, ds, p-1, dn, prec);       // dh += z^(p-1) ds
  }
}

// escape(k, g) + g
void implicit(Poly& F, Poly& dF, const int k, const Poly& g, const Poly& dg,
              const int n, const int dn, const int prec) {
  escape(F, dF, k, g, dg, n, dn, prec);
  arb_poly_add_series(F, F, g, n, prec);
  arb_poly_add_series(dF, dF, dg, dn, prec);
}

void area(arb_t mu, const Poly& f, const int prec) {
  arb_zero(mu);
  const auto n = f.length();
  Arb t;
  for (slong i = 0; i < n; i++) {
    // mu += (1-i) f[i]^2
    arb_poly_get_coeff_arb(t, f, i);
    arb_sqr(t, t, prec);
    arb_addmul_si(mu, t, 1-i, prec);
  }
  arb_const_pi(t, prec);
  arb_mul(mu, mu, t, prec);
}

void arb_areas(const int max_k, const int prec) {
  print("prec = %d (%d digits)\n\n", prec, int(prec*log10(2)));

  // f = 1, so g = log f = 0
  Poly g;
  print("k 0:\n  f = 1\n  g = 0");

  // dg = 1
  const Poly dg(1);

  Poly g0, F, dF, ignore;
  for (int k = 1; k <= max_k; k++) {
    print("\nk %d:", k);
    for (int refine = 0; refine < 2; refine++) {
      print("  refine %d:", refine);
      const auto start = wall_time();
      const int p = 1 << k;
      const int dp = refine ? p : p / 2;

      if (refine) {
        // Newton update all terms
        poly_mid(g0, g);
        implicit(F, ignore, k, g0, dg, p, 0, prec);
        implicit(ignore, dF, k, g, dg, p, p, prec);
        poly_div_refine(F, F, dF, p, prec);
        poly_intersect_sub(g, g0, F, p, prec);
      } else {
        // Newton update only the high terms
        implicit(F, dF, k, g, dg, p, dp, prec);
        F.assert_low_zero(dp);
        F >>= dp;
        poly_div_refine(F, F, dF, dp, prec);
        F <<= dp;
        arb_poly_sub_series(g, g, F, p, prec);
      }
      const auto elapsed = wall_time() - start;

      // Report results
      Poly f;
      poly_exp_refine(f, g, p, prec);
      Arb mu;
      area(mu, f, prec);
      print("    mu = %.10g %s", mu, mu.safe());
      if (k < 4) {
        print("    f = %.3g", f);
        print("    g = %.3g", g);
      }
      if (0) stats(f, "f");
      print("    time = %.3g s", elapsed.seconds());

      // Check against known results
      const span<const Known> knowns(known_areas);
      if (k < int(knowns.size())) {
        Arb known;
        arb_set_str(known, knowns[k].value, prec);
        slow_assert(arb_overlaps(mu, known), "No overlap with known area %.20g", known);
      }
    }
  }
}

}  // namespace mandelbrot
