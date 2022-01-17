// Mandelbrot area via Arb

#include "known.h"
#include "poly.h"
#include "print.h"
#include "wall_time.h"
namespace mandelbrot {

using std::min;
using std::max;
using std::runtime_error;

// C = a*A - 2^(log2_b)*B
void poly_linear_sub_si_2exp(Poly& C, const slong a, const Poly& A, const slong log2_b, const Poly& B, slong prec) {
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

void escape(Poly& h, Poly& dh, const int k, const Poly& f, const Poly& df,
            const int n, const int dn, const int prec) {
  // h = escape(k, z*f)^(2^-k) in mandelbrot-area-cupy.ipynb notation
  const bool verbose = false;

  // Handle h = 1 base case specially to avoid log(0)
  if (k == 0) {
    arb_poly_one(h);
    arb_poly_zero(dh);
    return;
  }

  // log_f = log(f)
  Poly log_f, dlog_f;
  poly_log_refine(log_f, f, n, prec);
  poly_div_refine(dlog_f, df, f, dn, prec);

  Poly log_h, dlog_h, t, dt, s, ds, u;
  for (int i = 1; i <= k; i++) {
    const int p = 1 << i;

    // t = (p-1)*log_f - p*log_h
    poly_linear_sub_si_2exp(t, p-1, log_f, i, log_h, prec);
    poly_linear_sub_si_2exp(dt, p-1, dlog_f, i, dlog_h, prec);
    if (verbose) stats(t, "linear %d", i);

    // log_h = log_h + (1/p) log(1 + z^(p-1)exp(t))
    poly_log1p_exp_shift_refine(s, t, p-1, n-(p-1), prec);     // s = z^(1-p) log(1 + z^(p-1) exp(t))
    poly_sub_shift_series(u, t, s, p-1, dn-(p-1), prec);       // u = t - z^(p-1) s
    poly_exp_refine(t, u, dn-(p-1), prec);                     // u = exp(t - z^(p-1) s)
    safe_poly_mullow(ds, t, dt, dn-(p-1), prec);               // ds = exp(t - z^(p-1) s) dt
    arb_poly_scalar_mul_2exp_si(s, s, -i);                     // s /= p
    arb_poly_scalar_mul_2exp_si(ds, ds, -i);                   // ds /= p
    poly_add_shift_series(log_h, log_h, s, p-1, n, prec);      // log_h += z^(p-1) s
    poly_add_shift_series(dlog_h, dlog_h, ds, p-1, dn, prec);  // dlog_h += z^(p-1) ds
  }

  // h = exp(log_h)
  poly_exp_refine(h, log_h, n, prec);
  arb_poly_mullow(dh, h, dlog_h, dn, prec);
}

void implicit(Poly& F, Poly& dF, const int k, const Poly& f, const Poly& df,
              const int n, const int dn, const int prec) {
  // escape(k, z/f) - 1/f

  // fi = 1 / f
  Poly fi, dfi;
  poly_inv_refine(fi, f, n, prec);
  arb_poly_mullow(dfi, fi, fi, dn, prec);
  arb_poly_mullow(dfi, dfi, df, dn, prec);
  arb_poly_neg(dfi, dfi);

  // escape(k, fi) - fi
  escape(F, dF, k, fi, dfi, n, dn, prec);
  arb_poly_sub_series(F, F, fi, n, prec);
  arb_poly_sub_series(dF, dF, dfi, dn, prec);
}

void area(arb_t mu, const Poly& f, const int prec) {
  arb_zero(mu);
  const int n = f.length();
  Arb t;
  for (int i = 0; i < n; i++) {
    // mu += (1-i) f[i]^2
    arb_poly_get_coeff_arb(t, f, i);
    arb_sqr(t, t, prec);
    arb_addmul_si(mu, t, 1-i, prec);
  }
  arb_const_pi(t, prec);
  arb_mul(mu, mu, t, prec);
}

void toplevel() {
  const int prec = 200;
  const int repeats = 2;
  print("prec = %d (%d digits)\n\n", prec, int(prec*log10(2)));

  // f = 1
  Poly f(1);
  print("k 0:\n  f = %.3g", f);

  // df = 1
  Poly df(1);

  Poly f0, F, dF, ignore;
  for (int k = 1; k <= 14; k++) {
    print("\nk %d:", k);
    for (int refine = 0; refine < repeats; refine++) {
      print("  refine %d:", refine);
      const auto start = wall_time();
      const int p = 1 << k;
      const int dp = refine ? p : p / 2;

      if (refine) {
        // Newton update all terms
        poly_mid(f0, f);
        implicit(F, ignore, k, f0, df, p, 0, prec);
        implicit(ignore, dF, k, f, df, p, p, prec);
        poly_div_refine(F, F, dF, p, prec);
        poly_intersect_sub(f, f0, F, p, prec);
      } else {
        // Newton update only the high terms
        implicit(F, dF, k, f, df, p, dp, prec);
        F.assert_low_zero(dp);
        F >>= dp;
        poly_div_refine(F, F, dF, dp, prec);
        F <<= dp;
        arb_poly_sub_series(f, f, F, p, prec);
      }
      const auto elapsed = wall_time() - start;

      // Report results
      Arb mu;
      area(mu, f, prec);
      print("    mu = %.10g %s", mu, mu.safe(20));
      if (k < 4)
        print("    f = %.3g", f);
      if (0) stats(f, "f");
      print("    time = %.3g s", elapsed.seconds());

      // Check against known results
      if (k < known_ks) {
        Arb known;
        arb_set_str(known, known_areas[k], prec);
        if (!arb_overlaps(mu, known))
          throw runtime_error(format("No overlap with known area %.20g", known));
      }
    }
  }
}

}  // namespace mandelbrot

int main() {
  try {
    mandelbrot::toplevel();
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
