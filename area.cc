// Mandelbrot area via Arb

#include "known.h"
#include "poly.h"
#include "print.h"
#include "wall_time.h"
namespace mandelbrot {

using std::min;
using std::max;
using std::runtime_error;

void poly_mid(Poly& mid, const Poly& f) {
  const slong n = f.length();
  arb_poly_zero(mid);
  arb_poly_fit_length(mid, n);
  _arb_poly_set_length(mid, n);
  Arb m;
  for (int i = 0; i < n; i++) {
    arb_get_mid_arb(m, f[i]);
    arb_poly_set_coeff_arb(mid, i, m);
  }
}

void intersect_sub(arb_t z, const arb_t x, const arb_t y, const slong prec) {
  Arb xy;
  arb_sub(xy, x, y, prec);
  if (!arb_overlaps(z, xy))
    throw runtime_error("intersect_sub failure");
  arb_intersection(z, z, xy, prec);
}

void poly_intersect_sub(Poly& C, const Poly& A, const Poly& B, const slong n, const slong prec) {
  arb_poly_zero(C);
  arb_poly_fit_length(C, n);
  _arb_poly_set_length(C, n);
  for (int i = 0; i < n; i++)
    intersect_sub(C.mutable_at(i), A[i], B[i], prec);
}

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

void escape(arb_poly_t h, arb_poly_t dh, const int k, const arb_poly_t f, const arb_poly_t df,
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
  arb_poly_log_series(log_f, f, n, prec);
  arb_poly_div_series(dlog_f, df, f, dn, prec);

  Poly log_h, dlog_h, t, dt, s, ds;
  for (int i = 1; i <= k; i++) {
    const int p = 1 << i;

    // t = (p-1)*log_f - p*log_h
    poly_linear_sub_si_2exp(t, p-1, log_f, i, log_h, prec);
    poly_linear_sub_si_2exp(dt, p-1, dlog_f, i, dlog_h, prec);
    if (verbose) stats(t, "linear %d", i);

    // log_h = log_h + (1/p) log(1 + z^(p-1)exp(t))
    arb_poly_exp_series(s, t, n-(p-1), prec);           // s = exp(t)
    if (verbose) stats(s, "exp(t) %d", i);
    safe_poly_mullow(ds, s, dt, dn-(p-1), prec);        // ds = exp(t) dt
    arb_poly_shift_left(s, s, p-1);                     // s = z^(p-1) exp(t)
    arb_poly_shift_left(ds, ds, p-1);                   // ds = z^(p-1) exp(t) dt
    arb_poly_log1p_series(t, s, n, prec);               // t = log(1 + z^(p-1) exp(t))
    if (verbose) stats(t, "log1p %d", i);
    arb_poly_add_si(s, s, 1, prec);                     // s = 1 + z^(p-1) exp(t)
    arb_poly_div_series(dt, ds, s, dn, prec);           // dt = d/df log(1 + z^(p-1) exp(t))
    arb_poly_scalar_mul_2exp_si(t, t, -i);              // t /= p
    arb_poly_scalar_mul_2exp_si(dt, dt, -i);            // dt /= p
    if (verbose) stats(t, "delta %d", i);
    arb_poly_add_series(log_h, log_h, t, n, prec);      // log_h += t
    arb_poly_add_series(dlog_h, dlog_h, dt, dn, prec);  // dlog_h += dt
    if (verbose) stats(log_h, "log_h %d", i);
  }

  // h = exp(log_h)
  arb_poly_exp_series(h, log_h, n, prec);
  arb_poly_mullow(dh, h, dlog_h, dn, prec);
}

void implicit(arb_poly_t F, arb_poly_t dF, const int k, const arb_poly_t f, const arb_poly_t df,
              const int n, const int dn, const int prec) {
  // escape(k, z/f) - 1/f

  // fi = 1 / f
  Poly fi, dfi;
  arb_poly_inv_series(fi, f, n, prec);
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
  printf("prec %d\n\n", prec);

  // f = 1
  Poly f(1);
  print("k 0:\n  f = %.3g", f);

  // df = 1
  Poly df(1);

  Poly f0, F, dF, ignore;
  for (int k = 1; k < 13; k++) {
    print("\nk %d:", k);
    for (int refine = 0; refine < repeats; refine++) {
      print("  refine %d:", refine);
      const auto start = wall_time();
      const int p = 1 << k;
      const int dp = refine ? p : p / 2;

      if (refine) {
        // Newton update all terms
        poly_mid(f0, f);
        implicit(ignore, dF, k, f0, df, p, p, prec);
        implicit(F, ignore, k, f0, df, p, 0, prec);
        arb_poly_div_series(F, F, dF, p, prec);
        poly_intersect_sub(f, f0, F, p, prec);
      } else {
        // Newton update only the high terms
        implicit(F, dF, k, f, df, p, dp, prec);
        F.assert_low_zero(dp);
        F >>= dp;
        arb_poly_div_series(F, F, dF, dp, prec);
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
    std::cerr << e.what() << std::endl;
    return 1;
  }
  return 0;
}
