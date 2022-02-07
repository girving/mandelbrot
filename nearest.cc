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

// exp(2𝜋i a/b) for a ∈ [0,zs.size()).
template<class S> int64_t nearest_twiddles(span<Complex<S>> zs, const int64_t b, const int fast_prec) {
  // Write n <= n0 * n1, so that j = j0*n1 + j1
  const int64_t n = zs.size();
  const auto n0 = int64_t(ceil(sqrt(double(n))));
  const auto n1 = (n + n0 - 1) / n0;

  int64_t fallbacks = 0;
  vector<Acb> factors(n0 + n1);  // Each twiddle(j0*n1, b), then each twiddle(j1, b)
  #pragma omp parallel
  {
    // Compute j0*n1 and j1 twiddles
    Fmpq t;
    #pragma omp for
    for (int64_t j = 0; j < n0+n1; j++) {
      const auto j0 = j, j1 = j - n0;
      fmpq_set_si(t, j < n0 ? 2*j0*n1 : 2*j1, b);
      cis_pi(factors[j], t, fast_prec);
    }

    // First compute factored using low precision, filling in holes using nonfactored evaluation
    Acb z;
    int64_t thread_fallbacks = 0;
    #pragma omp for
    for (int64_t j = 0; j < n; j++) {
      const auto j0 = j / n1, j1 = j - j0*n1;
      acb_mul(z, factors[j0], factors[n0 + j1], fast_prec);
      if (auto r = round_nearest<S>(z, fast_prec))
        zs[j] = *r;
      else {
        zs[j] = nearest_twiddle<S>(j, b);
        thread_fallbacks++;
      }
    }
    #pragma omp atomic
    fallbacks += thread_fallbacks;
  }
  return fallbacks;
}

#define NEAREST(S) \
  template S nearest_pi(); \
  template S nearest_sqrt(const int64_t, const int64_t); \
  template Complex<S> nearest_twiddle(const int64_t, const int64_t); \
  template int64_t nearest_twiddles(span<Complex<S>>, const int64_t, const int fast_prec);
NEAREST(double)

#define EXPANSION(n) \
  template struct RoundNear<Expansion<n>>; \
  NEAREST(Expansion<n>)
EXPANSION(2)
EXPANSION(3)

}  // namespace mandelbrot
