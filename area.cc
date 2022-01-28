// Mandelbrot area via custom power series

#include "area.h"
#include "arb_cc.h"
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

// z = a*x - 2^(log2_b)*y
SERIES_EXP(sub_si_2exp, z, (class A,class B), (a,x,b,y),
          (const int64_t a, const Series<A>& x, const int b, const Series<B>& y)) {
  const auto nk = min(x.known(), y.known());
  const auto nz = min(nk, max(x.nonzero(), y.nonzero()));
  const auto both = min(x.nonzero(), y.nonzero());
  const auto x_only = min(nz, x.nonzero());
  const auto y_only = min(nz, y.nonzero());
  const auto zp = z.data();
  for (int64_t i = 0; i < both; i++)
    zp[i] = a*x[i] - ldexp(y[i], b);
  for (int64_t i = both; i < x_only; i++)
    zp[i] = a*x[i];
  for (int64_t i = both; i < y_only; i++)
    zp[i] = -ldexp(y[i], b);
  z.set_counts(nk, nz);
}

// h = log escape(k, z*e^-g)^(2^-k) in mandelbrot-area-cupy.ipynb notation
template<class S> void escape(Series<S>& h, Series<S>& dh, const int k,
                              const Series<const S>& g, const Series<const S>& dg,
                              const int64_t n, const int64_t dn) {
  // Base case
  h.set_scalar(n, 0);
  dh.set_scalar(dn, 0);

  Series<S> t(n), dt(dn), s(n), u(dn);
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
template<class S> void implicit(Series<S>& F, Series<S>& dF, const int k,
                                const Series<const S>& g, const Series<const S>& dg,
                                const int64_t n, const int64_t dn) {
  escape(F, dF, k, g, dg, n, dn);
  F += g;
  dF += dg;
}

template<class A> remove_const_t<A> area(const Series<A>& f) {
  typedef remove_const_t<A> S;
  S mu = 0;
  const auto nz = f.nonzero();
  for (int64_t i = 0; i < nz; i++)
    mu += (1-i) * sqr(f[i]);
  return nearest_pi<S>() * mu;
}

template<class S> void areas(const int max_k) {
  // f = 1, so g = log f = 0
  Series<S> g(1);
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
      Series<S> dg(dp);
      dg.set_scalar(dp, 1);

      if (refine) {
        // Newton update all terms
        static_assert(!is_interval<S>);
        const auto& g0 = g;  // Valid until S is an interval type
        Series<S> F(p), dF(p), ignore(0);
        implicit<S>(F, dF, k, g, dg, p, p);
        implicit<S>(F, ignore, k, g0, dg, p, 0);
        dg = div(F, dF);
        g.high_sub(p/2, dg.high(p/2));
      } else {
        // Newton update only the high terms
        Series<S> F(p), dF(dp);
        implicit<S>(F, dF, k, g, dg, p, dp);
        dg = div(F.high(dp), dF);
        g.high_sub(dp, dg);
      }
      const auto elapsed = wall_time() - start;

      // Report results
      Series<S> f(p);
      f = exp(g);
      const S mu = area(f);
      print("    mu = %.10g", mu);
      if (k < 4) {
        print("    f = %.3g", f);
        print("    g = %.3g", g);
      }
      print("    time = %.3g s", elapsed.seconds());

      // Check against known results
      if (k < known_ks) {
        Arb known;
        arb_set_str(known, known_areas[k], 200);
        const S k = round_near<S>(arb_midref(known.x));
        const S e = abs(mu - k);
        print("    error = %.3g", e);
        const S tol = 1e-6;
        slow_assert(e <= tol, "error %g > %g", e, tol);
      }
    }
  }
}

#define AREAS(S) \
  template void areas<S>(const int max_k);
AREAS(double)

}  // namespace mandelbrot
