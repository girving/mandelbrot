/* Mandelbrot area via Arb */

#include <arb.h>
#include <arb_poly.h>

void safe_poly_mullow(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, slong n, slong prec) {
  if (n <= 0)
    arb_poly_zero(C);
  else
    arb_poly_mullow(C, A, B, n, prec);
}

void escape(arb_poly_t h, arb_poly_t dh, const int k, const arb_poly_t f, const arb_poly_t df, const int prec) {
  /* h = escape(k, z*f)^(2^-k) in mandelbrot-area-cupy.ipynb notation */
  const int n = 1 << k;
  const int dn = n / 2;

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

    /* log_h = log_h + (1/p) log(1 + z^(p-1)exp(t)) */
    arb_poly_exp_series(s, t, n-(p-1), prec);           /* s = exp(t) */
    safe_poly_mullow(ds, s, dt, dn-(p-1), prec);        /* ds = exp(t) dt */
    arb_poly_shift_left(s, s, p-1);                     /* s = z^(p-1) exp(t) */
    arb_poly_shift_left(ds, ds, p-1);                   /* ds = z^(p-1) exp(t) dt */
    arb_poly_log1p_series(t, s, n, prec);               /* t = log(1 + z^(p-1) exp(t)) */
    arb_poly_add_si(s, s, 1, prec);                     /* s = 1 + z^(p-1) exp(t) */
    arb_poly_div_series(dt, ds, s, dn, prec);           /* dt = d/df log(1 + z^(p-1) exp(t)) */
    arb_poly_scalar_mul_2exp_si(t, t, -(i+1));          /* t /= p */
    arb_poly_scalar_mul_2exp_si(dt, dt, -(i+1));        /* dt /= p */
    arb_poly_add_series(log_h, log_h, t, n, prec);      /* log_h += t */
    arb_poly_add_series(dlog_h, dlog_h, dt, dn, prec);  /* dlog_h += dt */
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

/* Known areas, computed with precision 2000 */
const int known_ks = 15+1;
const char* known_areas[known_ks] = {
  "3.1415926535897932385 +/- 3.74e-20",  /* k 0 */
  "3.1415926535897932385 +/- 3.74e-20",  /* k 1 */
  "2.6998061866787285643 +/- 3.84e-21",  /* k 2 */
  "2.4636540388624509718 +/- 1.11e-20",  /* k 3 */
  "2.3089221574735649437 +/- 1.86e-20",  /* k 4 */
  "2.1863969605631909502 +/- 3.77e-20",  /* k 5 */
  "2.1117843557698974143 +/- 3.66e-20",  /* k 6 */
  "2.0290979917037465762 +/- 2.62e-20",  /* k 7 */
  "1.9793859246155479841 +/- 4.66e-21",  /* k 8 */
  "1.9277162229003024017 +/- 3.17e-20",  /* k 9 */
  "1.8959430758033194112 +/- 3.71e-20",  /* k 10 */
  "1.8546567767819607804 +/- 1.97e-20",  /* k 11 */
  "1.8346557326261556582 +/- 3.03e-20",  /* k 12 */
  "1.8061788858565225232 +/- 4.46e-20",  /* k 13 */
  "1.7863897165519852571 +/- 2.07e-20",  /* k 14 */
  "1.7668376741940830366 +/- 2.45e-20",  /* k 15 */
};

int main() {
  const int prec = 2000;
  printf("prec %d\n\n", prec);

  /* Prepare temporaries */
  arb_t a, mu, known, F0, r;
  arb_init(a);
  arb_init(mu);
  arb_init(known);
  arb_init(F0);
  arb_init(r);
  arb_poly_t f, df, F, dF;
  arb_poly_init(f);
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

  for (int k = 1;; k++) {
    printf("\nk %d:\n", k);
    const int p = 1 << k;
    const int dp = p / 2;

    /* Evaluate implicit equation and its derivative */
    implicit(F, dF, k, f, df, prec);

    /* Verify low terms are zero, and shift away */
    for (int i = 0; i < dp; i++) {
      arb_poly_get_coeff_arb(a, F, i);
      if (!arb_contains_si(a, 0)) {
        printf("Newton failure\n");
        exit(1);
      }
    }
    arb_poly_shift_right(F, F, dp);
    arb_poly_get_coeff_arb(F0, F, 0);

    /* Newton step */
    arb_poly_div_series(F, F, dF, dp, prec);
    arb_poly_shift_left(F, F, dp);
    arb_poly_sub_series(f, f, F, p, prec);

    /* Report results */
    area(mu, f, prec);
    arb_get_rad_arb(r, mu);
    printf("  mu = ");
    arb_printn(mu, 20, 0);
    printf(" (radius ");
    arf_printd(arb_midref(r), 3);
    printf(")\n  F = [");
    arb_printd(F0, 3);
    printf("]z^%d + O(z^%d)\n", dp, dp+1);
    if (k < 4) {
      printf("  f = ");
      arb_poly_printd(f, 3);
      printf("\n");
    }

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

  /* Clean temporaries */
  arb_clear(a);
  arb_clear(mu);
  arb_clear(F0);
  arb_clear(known);
  arb_clear(r);
  arb_poly_clear(f);
  arb_poly_clear(df);
  arb_poly_clear(F);
  arb_poly_clear(dF);

  return 0;
}
