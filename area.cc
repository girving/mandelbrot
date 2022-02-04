// Mandelbrot area via custom power series

#include "area.h"
#include "arb_cc.h"
#include "device.h"
#include "expansion.h"
#include "nearest.h"
#include "debug.h"
#include "known.h"
#include "print.h"
#include "series.h"
#include "wall_time.h"
namespace mandelbrot {

using std::min;
using std::max;
using std::runtime_error;
using std::swap;

// z = a*x - 2^(log2_b)*y
DEF_LOOP(sub_si_2exp_loop, nz, i, (S* z, const int a, const S* x, const int xnz, const int b, const S* y, const int ynz),
  const auto xi = i < xnz ? x[i] : S(0);
  const auto yi = i < ynz ? y[i] : S(0);
  z[i] = a*xi - ldexp(yi, b);)
SERIES_EXP(sub_si_2exp, z, (class SA,class SB), (a,x,b,y), (a,x=x.view(),b,y=y.view()),
          (const int64_t a, const SA& x, const int b, const SB& y)) {
  const auto nk = min(x.known(), y.known());
  const auto nz = min(nk, max(x.nonzero(), y.nonzero()));
  z.set_counts(nk, nz);
  sub_si_2exp_loop(nz, z.data(), a, x.data(), x.nonzero(), b, y.data(), y.nonzero());
}

// h = log escape(k, z*e^-g)^(2^-k) in mandelbrot-area-cupy.ipynb notation
template<class T> void escape(Series<T>& h, Series<T>& dh, const int k,
                              SeriesView<const T> g, SeriesView<const T> dg,
                              const int64_t n, const int64_t dn) {
  // Base case
  h.set_scalar(n, 0);
  dh.set_scalar(dn, 0);

  Series<T> t(n), dt(dn), s(n), u(dn);
  for (int i = 1; i <= k; i++) {
    const auto p = int64_t(1) << i;

    // t = (1-p)g - p h
    t = sub_si_2exp(1-p, g, i, h);
    dt = sub_si_2exp(1-p, dg, i, dh.low(dn-(p-1)));

    // h += (1/p) log(1 + z^(p-1)exp(t))
    s = log1p_exp(t.low(n-(p-1)), p-1);  // s = z^(1-p) log(1 + z^(p-1) exp(t))
    u = t.low(dn-(p-1));
    u.high_sub(p-1, s);                  // u = t - z^(p-1) s
    t = exp(u);                          // t = exp(t - z^(p-1) s)
    dt = mul(t, dt.low(dn-(p-1)));       // ds = exp(t - z^(p-1) s) dt
    s = ldexp(s, -i);                    // s /= p
    dt = ldexp(dt, -i);                  // ds /= p
    h.high_add(p-1, s);                  // h += z^(p-1) s
    dh.high_add(p-1, dt);                // h += z^(p-1) ds
  }
}

// escape(k, g) + g
template<class T> void implicit(Series<T>& F, Series<T>& dF, const int k,
                                SeriesView<const T> g, SeriesView<const T> dg,
                                const int64_t n, const int64_t dn) {
  escape(F, dF, k, g, dg, n, dn);
  F += g;
  dF += dg;
}

template<class A> typename Series<A>::Scalar area(const Series<A>& f) {
  typedef typename Series<A>::Scalar S;
  S mu = 0;
  const auto nz = f.nonzero();
  const auto& hf = host_copy(f);
  for (int64_t i = 0; i < nz; i++)
    mu += (1-i) * sqr(hf[i]);
  return nearest_pi<S>() * mu;
}

template<class T> void areas(const int max_k, const double tol) {
  typedef Undevice<T> S;

  // f = 1, so g = log f = 0
  Series<T> g(1);
  g.set_scalar(1, 0);
  print("k 0:\n  f = 1\n  g = 0");

  for (int k = 1; k <= max_k; k++) {
    print("\nk %d:", k);
    for (int refine = 0; refine < 2; refine++) {
      print("  refine %d:", refine);
      const auto start = wall_time();
      const int p = 1 << k;
      const int dp = refine ? p : p / 2;

      // Reallocate and extend
      g.copy(p).swap(g);
      g.set_known(p);

      // dg = 1
      Series<T> dg(dp);
      dg.set_scalar(dp, 1);

      if (refine) {
        // Newton update all terms
        static_assert(!is_interval<S>);
        const auto& g0 = g;  // Valid until S is an interval type
        Series<T> F(p), dF(p), ignore(0);
        implicit<T>(F, dF, k, g, dg, p, p);
        implicit<T>(F, ignore, k, g0, dg, p, 0);
        dg = div(F, dF);
        g.high_sub(p/2, dg.high(p/2));
      } else {
        // Newton update only the high terms
        Series<T> F(p), dF(dp);
        implicit<T>(F, dF, k, g, dg, p, dp);
        dg = div(F.high(dp), dF);
        g.high_sub(dp, dg);
      }
      const auto elapsed = wall_time() - start;

      // Report results
      Series<T> f(p);
      f = exp(g);
      const S mu = area(f);
      print("    mu = %s", safe(mu));
      if (k < 4) {
        print("    f = %.3g", host_copy(f));
        print("    g = %.3g", host_copy(g));
      }
      print("    time = %.3g s", elapsed.seconds());

      // Check against known results
      const span<const Known> knowns(known_areas);
      if (k < int(knowns.size())) {
        const int prec = 1000;
        Arb known, error;
        arb_set_str(known, knowns[k].value, prec);
        const Arb ours = exact_arb(mu);
        arb_sub(error, known, ours, prec);
        const double e = bound(error);
        print("    error = %.3g", e);
        const auto t = refine ? tol : 1e-6;
        slow_assert(e <= t, "error %g > %g", e, t);
      }
    }
  }
}

#define AREAS(T) \
  template void areas<T>(const int max_k, const double);
AREAS(double)
AREAS(Expansion<2>)
IF_CUDA(
  AREAS(Device<double>)
  AREAS(Device<Expansion<2>>)
)

}  // namespace mandelbrot
