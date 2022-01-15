/* Mandelbrot area via Arb */

#include <arb.h>
#include <arb_poly.h>
#include <stdarg.h>
#include <stdbool.h>
#include "known.h"

void stats(const arb_poly_t f, const char* format, ...) {
  printf("    ");
  va_list args;
  va_start(args, format);
  vfprintf(stdout, format, args);
  va_end(args);

  arb_t c, min, max, radius;
  arb_init(c);
  arb_init(min);
  arb_init(max);
  arb_init(radius);

  const int prec = 1000;
  const slong degree = arb_poly_degree(f);
  arb_pos_inf(min);
  arb_neg_inf(max);
  for (int i = 0; i <= degree; i++) {
    arb_poly_get_coeff_arb(c, f, i);
    arb_min(min, min, c, prec);
    arb_max(max, max, c, prec);
    arb_get_rad_arb(c, c);
    arb_max(radius, radius, c, prec);
  }

  printf(": degree %ld, coeffs [", degree);
  arb_printd(min, 3);
  printf(", ");
  arb_printd(max, 3);
  printf("], radius ");
  arb_printd(radius, 3);
  printf("\n");

  arb_clear(c);
  arb_clear(min);
  arb_clear(max);
  arb_clear(radius);
}

void safe_poly_mullow(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, slong n, slong prec) {
  if (n <= 0)
    arb_poly_zero(C);
  else
    arb_poly_mullow(C, A, B, n, prec);
}

void poly_mid(arb_poly_t mid, const arb_poly_t f) {
  const slong n = arb_poly_length(f);
  arb_poly_fit_length(mid, n);
  arb_t m;
  arb_init(m);
  for (int i = 0; i < n; i++) {
    arb_get_mid_arb(m, arb_poly_get_coeff_ptr(f, i));
    arb_poly_set_coeff_arb(mid, i, m);
  }
  arb_clear(m);
}

void intersect_sub(arb_t z, const arb_t x, const arb_t y, const slong prec) {
  arb_t xy;
  arb_init(xy);
  arb_sub(xy, x, y, prec);
  if (!arb_overlaps(z, xy)) {
    printf("intersect_sub failure\n");
    exit(1);
  }
  arb_intersection(z, z, xy, prec);
  arb_clear(xy);
}

void poly_intersect_sub(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, const slong n, const slong prec) {
  arb_t c, a, b;
  arb_init(c);
  arb_init(a);
  arb_init(b);
  for (int i = 0; i < n; i++) {
    arb_poly_get_coeff_arb(c, C, i);
    arb_poly_get_coeff_arb(a, A, i);
    arb_poly_get_coeff_arb(b, B, i);
    intersect_sub(c, a, b, prec);
    arb_poly_set_coeff_arb(C, i, c);
  }
  arb_clear(c);
  arb_clear(a);
  arb_clear(b);
}

void escape(arb_poly_t h, arb_poly_t dh, const int k, const arb_poly_t f, const arb_poly_t df, const int prec) {
  /* h = escape(k, z*f)^(2^-k) in mandelbrot-area-cupy.ipynb notation */
  const int n = 1 << k;
  const int dn = n / 2;
  const bool verbose = false;

  /* Handle h = 1 base case specially to avoid log(0) */
  if (k == 0) {
    arb_poly_one(h);
    arb_poly_zero(dh);
    return;
  }

  /* Prepare temporaries */
  arb_t a;
  arb_init(a);
  arb_poly_t log_h, dlog_h, log_f, dlog_f, t, dt, s, ds;
  arb_poly_init(log_h);
  arb_poly_init(dlog_h);
  arb_poly_init(log_f);
  arb_poly_init(dlog_f);
  arb_poly_init(t);
  arb_poly_init(dt);
  arb_poly_init(s);
  arb_poly_init(ds);

  /* log_f = log(f) */
  arb_poly_log_series(log_f, f, n, prec);
  arb_poly_div_series(dlog_f, df, f, dn, prec);

  for (int i = 0; i < k; i++) {
    const int p = 1 << (i+1);

    /* t = (p-1)*log_f - p*log_h */
    arb_set_ui(a, p-1);
    arb_poly_scalar_mul(t, log_f, a, prec);
    arb_poly_scalar_mul(dt, dlog_f, a, prec);
    arb_set_ui(a, p);
    arb_poly_scalar_mul(s, log_h, a, prec);
    arb_poly_scalar_mul(ds, dlog_h, a, prec);
    arb_poly_sub_series(t, t, s, n, prec);
    arb_poly_sub_series(dt, dt, ds, dn, prec);
    if (verbose) stats(t, "linear %d", i);

    /* log_h = log_h + (1/p) log(1 + z^(p-1)exp(t)) */
    arb_poly_exp_series(s, t, n-(p-1), prec);           /* s = exp(t) */
    if (verbose) stats(s, "exp(t) %d", i);
    if (0 && k == 11 && i == 2) {
      printf("\nt = ");
      arb_poly_printd(t, 3);
      printf("\n\nexp(t) = ");
      arb_poly_printd(s, 3);
      printf("\n\n");
    }
    safe_poly_mullow(ds, s, dt, dn-(p-1), prec);        /* ds = exp(t) dt */
    arb_poly_shift_left(s, s, p-1);                     /* s = z^(p-1) exp(t) */
    arb_poly_shift_left(ds, ds, p-1);                   /* ds = z^(p-1) exp(t) dt */
    arb_poly_log1p_series(t, s, n, prec);               /* t = log(1 + z^(p-1) exp(t)) */
    if (verbose) stats(t, "log1p %d", i);
    arb_poly_add_si(s, s, 1, prec);                     /* s = 1 + z^(p-1) exp(t) */
    arb_poly_div_series(dt, ds, s, dn, prec);           /* dt = d/df log(1 + z^(p-1) exp(t)) */
    arb_poly_scalar_mul_2exp_si(t, t, -(i+1));          /* t /= p */
    arb_poly_scalar_mul_2exp_si(dt, dt, -(i+1));        /* dt /= p */
    if (verbose) stats(t, "delta %d", i);
    arb_poly_add_series(log_h, log_h, t, n, prec);      /* log_h += t */
    arb_poly_add_series(dlog_h, dlog_h, dt, dn, prec);  /* dlog_h += dt */
    if (verbose) stats(log_h, "log_h %d", i);
  }

  /* h = exp(log_h) */
  arb_poly_exp_series(h, log_h, n, prec);
  arb_poly_mullow(dh, h, dlog_h, dn, prec);

  /* Clear temporaries */
  arb_clear(a);
  arb_poly_clear(log_h);
  arb_poly_clear(dlog_h);
  arb_poly_clear(log_f);
  arb_poly_clear(dlog_f);
  arb_poly_clear(t);
  arb_poly_clear(dt);
  arb_poly_clear(s);
  arb_poly_clear(ds);
}

void implicit(arb_poly_t F, arb_poly_t dF, const int k, const arb_poly_t f, const arb_poly_t df, const int prec) {
  /* escape(k, z/f) - 1/f */
  const int n = 1 << k;
  const int dn = n / 2;

  /* Prepare temporaries */
  arb_poly_t fi, dfi;
  arb_poly_init(fi);
  arb_poly_init(dfi);

  /* fi = 1 / f */
  arb_poly_inv_series(fi, f, n, prec);
  arb_poly_mullow(dfi, fi, fi, dn, prec);
  arb_poly_mullow(dfi, dfi, df, dn, prec);
  arb_poly_neg(dfi, dfi);

  /* escape(k, fi) - fi */
  escape(F, dF, k, fi, dfi, prec);
  arb_poly_sub_series(F, F, fi, n, prec);
  arb_poly_sub_series(dF, dF, dfi, dn, prec);

  /* Clean temporaries */
  arb_poly_clear(fi);
  arb_poly_clear(dfi);
}

void area(arb_t mu, const arb_poly_t f, const int prec) {
  /* Prepare temporary */
  arb_t t;
  arb_init(t);

  arb_zero(mu);
  const int n = arb_poly_length(f);
  for (int i = 0; i < n; i++) {
    /* mu += (1-i) f[i]^2 */
    arb_poly_get_coeff_arb(t, f, i);
    arb_sqr(t, t, prec);
    arb_addmul_si(mu, t, 1-i, prec);
  }
  arb_const_pi(t, prec);
  arb_mul(mu, mu, t, prec);

  /* Clear temporary */
  arb_clear(t);
}

void checked_shift_right(arb_poly_t f, const arb_poly_t g, const slong n) {
  /* Verify low terms are zero */
  arb_t a;
  arb_init(a);
  for (int i = 0; i < n; i++) {
    arb_poly_get_coeff_arb(a, g, i);
    if (!arb_contains_si(a, 0)) {
      printf("checked_shift_right failure\n");
      exit(1);
    }
  }
  arb_clear(a);

  /* Shift the zero terms away */
  arb_poly_shift_right(f, g, n);
}

int main() {
  const int prec = 640;
  const int repeats = 2;
  printf("prec %d\n\n", prec);

  /* Prepare temporaries */
  arb_t a, mu, known, F0;
  arb_init(a);
  arb_init(mu);
  arb_init(known);
  arb_init(F0);
  arb_poly_t f, f0, df, F, dF;
  arb_poly_init(f);
  arb_poly_init(f0);
  arb_poly_init(df);
  arb_poly_init(F);
  arb_poly_init(dF);

  /* f = 1 */
  arb_poly_one(f);
  printf("k 0:\n  f = ");
  arb_poly_printd(f, 3);
  printf("\n");

  /* df = 1 */
  arb_poly_one(df);

  for (int k = 1; k < 32; k++) {
    printf("\nk %d:\n", k);
    for (int refine = 0; refine < repeats; refine++) {
      printf("  refine %d:\n", refine);
      const int p = 1 << k;
      const int dp = refine ? p : p / 2;

      /* Evaluate implicit equation and its derivative */
      if (refine)
        poly_mid(f0, f);
      implicit(F, dF, k, refine ? f0 : f, df, prec);
      arb_poly_get_coeff_arb(F0, F, dp);
      if (0) {
        stats(F, "F");
        stats(dF, "dF");
      }

      if (refine) {  /* Newton update all terms */
        arb_poly_div_series(F, F, dF, p, prec);
        poly_intersect_sub(f, f0, F, p, prec);
      } else {  /* Newton update only the high terms */
        checked_shift_right(F, F, dp);
        arb_poly_div_series(F, F, dF, dp, prec);
        arb_poly_shift_left(F, F, dp);
        arb_poly_sub_series(f, f, F, p, prec);
      }

      /* Report results */
      area(mu, f, prec);
      printf("    mu = ");
      arb_printd(mu, 10);
      printf(" ");
      arb_printn(mu, 20, 0);
      printf("\n    F = [");
      arb_printd(F0, 3);
      printf("]z^%d + O(z^%d)\n", dp, dp+1);
      if (k < 4) {
        printf("    f = ");
        arb_poly_printd(f, 3);
        printf("\n");
      }
      if (0) stats(f, "f");

      /* Check against known results */
      if (k < known_ks) {
        arb_set_str(known, known_areas[k], prec);
        if (!arb_overlaps(mu, known)) {
          printf("  ERROR: No overlap with known area ");
          arb_printd(known, 20);
          printf("\n");
          exit(1);
        }
      }
    }
  }

  /* Clean temporaries */
  arb_clear(a);
  arb_clear(mu);
  arb_clear(F0);
  arb_clear(known);
  arb_poly_clear(f);
  arb_poly_clear(f0);
  arb_poly_clear(df);
  arb_poly_clear(F);
  arb_poly_clear(dF);

  return 0;
}
