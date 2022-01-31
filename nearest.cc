// Accurate constants via arb

#include "nearest.h"
#include "acb_cc.h"
#include "arb_cc.h"
#include "expansion.h"
#include "fmpq_cc.h"
#include <vector>
namespace mandelbrot {

using std::vector;

template<class S> optional<S> round_nearest(const arb_t c, const int prec) {
  Arf a, b;
  arb_get_interval_arf(a, b, c, prec);
  const auto a_ = round_near<S>(a);
  const auto b_ = round_near<S>(b);
  optional<S> c_;
  if (a_ == b_)
    c_ = a_;
  return c_;
}
template<class S> optional<Complex<S>> round_nearest(const acb_t z, const int prec) {
  optional<Complex<S>> z_;
  const auto r = round_nearest<S>(acb_realref(z), prec);
  if (r) {
    const auto i = round_nearest<S>(acb_imagref(z), prec);
    if (i)
      z_ = Complex<S>(*r, *i);
  }
  return z_;
}

template<int n> Expansion<n> RoundNear<Expansion<n>>::round(const arf_t c) {
  const int max_prec = 1600;
  Expansion<n> e;
  Arb dc, t;
  for (int prec = 60*n; prec <= max_prec; prec <<= 1) {
    arb_set_arf(dc, c);
    for (int i = 0; i < n; i++) {
      const auto x = round_nearest<double>(dc, prec);
      if (!x) goto fail;
      e.x[i] = *x;
      arb_set_d(t, *x);
      arb_sub(dc, dc, t, prec);
    }
    return e;
    fail:;
  }
  die("ran out of precision rounding arf_t to Expansion<%d> (max prec = %d)", n, max_prec);
}

// 𝜋
template<class S> S nearest_pi() {
  static const S pi = nearest<S>([](const int prec) {
    Arb c; arb_const_pi(c, prec); return c;
  });
  return pi;
}

template<class S> S nearest_sqrt(const int64_t a, const int64_t b) {
  Arb r;
  return nearest<S>([a,b,&r](const int prec) {
    arb_set_si(r, a);
    arb_div_ui(r, r, b, prec);
    Arb s;
    arb_sqrt(s, r, prec);
    return s;
  });
}

static inline void cis_pi(Acb& z, const Fmpq& t, const int prec) {
  arb_sin_cos_pi_fmpq(acb_imagref(z.x), acb_realref(z.x), t, prec);
}

// exp(2𝜋i a/b)
template<class S> Complex<S> nearest_twiddle(const int64_t a, const int64_t b) {
  Fmpq t;
  fmpq_set_si(t, 2*a, b);
  return nearest<S>([t=move(t)](const int prec) {
    Acb z; cis_pi(z, t, prec); return z;
  });
}

// exp(2𝜋i a/b) for a ∈ [0,zs.size())
template<class S> void nearest_twiddles(span<Complex<S>> zs, const int64_t b) {
  const int max_prec = 1600;
  const int64_t n = zs.size();
  // Write n <= n0 * n1, so that j = j0*n1 + j1
  const auto n0 = int64_t(ceil(sqrt(double(n))));
  const auto n1 = (n + n0 - 1) / n0;
  Fmpq t;
  Acb z, zs0;
  vector<Acb> zs1(n1);
  for (int prec = 200; prec <= max_prec; prec <<= 1) {
    // Compute low twiddles
    for (int64_t j1 = 0; j1 < n1; j1++) {
      fmpq_set_si(t, 2*j1, b);
      cis_pi(zs1[j1], t, prec);
    }
    // Compute all twiddles
    for (int64_t j0 = 0; j0 < n0; j0++) {
      fmpq_set_si(t, 2*j0*n1, b);
      cis_pi(zs0, t, prec);
      for (int64_t j1 = 0; j1 < n1; j1++) {
        const auto j = j0*n1 + j1;
        if (j >= n) break;
        acb_mul(z, zs0, zs1[j1], prec);
        const auto z_ = round_nearest<S>(z, prec);
        if (!z_) goto fail;
        zs[j0*n1 + j1] = *z_;
      }
    }
    return;  // Success!
    fail:;  // Not enough precision
  }
  die("nearest_twiddles ran out of precision (max prec = %g)", max_prec);
}

#define NEAREST(S) \
  template S nearest_pi(); \
  template S nearest_sqrt(const int64_t, const int64_t); \
  template Complex<S> nearest_twiddle(const int64_t, const int64_t); \
  template void nearest_twiddles(span<Complex<S>>, const int64_t); 
NEAREST(double)

#define EXPANSION(n) \
  template struct RoundNear<Expansion<n>>; \
  NEAREST(Expansion<n>)
EXPANSION(2)
EXPANSION(3)

}  // namespace mandelbrot
