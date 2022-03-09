// Mandelbrot area via custom power series

#include "area.h"
#include "arb_cc.h"
#include "device.h"
#include "expansion_arith.h"
#include "nearest.h"
#include "debug.h"
#include "known.h"
#include "print.h"
#include "series.h"
#include "wall_time.h"
namespace mandelbrot {

using std::min;
using std::max;
using std::make_tuple;
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

template<class S> S area(SeriesView<const S> f) {
  // Sum via arb out of paranoia
  return nearest<S>([f](const int prec) {
    // mu = sum_i (1-i) f[i]^2
    Arb mu;
    const auto reduce = [prec](Arb& y, const Arb& x) {
      arb_add(y, y, x, prec);
    };
    const auto map = [f,prec](const int64_t i) {
      auto t = exact_arb(f[i]);
      arb_sqr(t, t, prec);
      arb_mul_si(t, t, 1-i, prec);
      return t;
    };
    map_reduce(mu, reduce, map, f.nonzero());
    Arb pi;
    arb_const_pi(pi, prec);
    arb_mul(mu, mu, pi, prec);
    return mu;
  });
}

// f = 1, so g = log f = 0
template<class T> Series<T> bottcher_base() {
  Series<T> g(1);
  g.set_scalar(1, 0);
  print("k 0:\n  f = 1\n  g = 0");
  return g;
}

// Determine k.  We assume it's a power of two.
int known_to_k(const int64_t known) {
  const int k = int(countr_zero(uint64_t(known)));
  slow_assert(known == int64_t(1) << k, "known = %d is not a power of two", known);
  return k;
}

template<class T> tuple<Series<T>,Undevice<T>> bottcher_step(Series<T>& g, const double tol) {
  typedef Undevice<T> S;

  // Determine k
  const int k0 = known_to_k(g.known());
  const int k = k0 + 1;
  const int64_t p = 1 << k;
  print("\nk %d:", k);

  // One step of Newton extension, then one step of Newton refinement
  // The loop runs for 2 iterations, but returns within for scoping reasons
  for (int refine = 0; /* refine < 2 */; refine++) {
    const auto start = wall_time();
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

    // Estimate area!
    Series<T> f(p);
    f = exp(g);
    const auto& hf = host_copy(f);
    const S mu = area<S>(hf);
    const auto elapsed = wall_time() - start;

    // Check against known results
    double error = 0;
    string error_s;
    const span<const Known> knowns(known_areas);
    if (k < int(knowns.size())) {
      const int prec = 1000;
      Arb known, error_a;
      arb_set_str(known, knowns[k].value, prec);
      const Arb ours = exact_arb(mu);
      arb_sub(error_a, known, ours, prec);
      error = bound(error_a);
      error_s = format(", error = %.3g", error);
    }

    // Report results
    print("  k %d, %.3g s: mu = %s%s", k, elapsed.seconds(), safe(mu), error_s);
    if (k < 4) {
      print("    f = %.3g", hf);
      print("    g = %.3g", host_copy(g));
    }

    // Report Bittner comparisons
    if (refine) {
      for (const auto& b : bittner_areas) {
        if (p/2 <= b.terms && b.terms <= p) {
          const auto ours = area<S>(hf.low(b.terms));
          const int prec = 1000;
          auto diff = exact_arb(ours);
          arb_sub(diff, diff, exact_arb(b.value), prec);
          print("  bittner %d = %.5g, ours = %s, error = %.3g", b.terms, b.value, safe(ours), bound(diff));
        }
      }
    }

    // Bail if we're inaccurate
    const auto goal = refine ? tol : 1e-6;
    slow_assert(error <= goal, "error %g > %g", error, goal);

    // Return both f and mu if we're done
    if (refine) return make_tuple(move(f), mu);
  }
}

template<class T> void areas(const int max_k, const double tol) {
  auto g = bottcher_base<T>();
  for (int k = 1; k <= max_k; k++)
    bottcher_step(g, tol);
}

template<class T> void write_bottcher(const string& output, const string& mode,
                                      const Undevice<T> mu, SeriesView<const T> f, SeriesView<const T> g) {
  if constexpr (is_device<T>)
    return write_bottcher<Undevice<T>>(output, mode, mu, host_copy(f), host_copy(g));
  else {
    const int k = known_to_k(f.known());
    slow_assert(f.known() == g.known());
    slow_assert(f[0] == 1 && g[0] == 0, "f[0] = %g, g[0] = %g", f[0], g[0]);  // Make sure f and g aren't flipped
    const auto write = [&output,&mode,k,mu=mu](const string& n, const string& name, const auto& x) {
      write_series(
          format("%s/%c-k%d", output, n, k),
          {name, format("mode = %s", mode), format("k = %d", k), format("mu = %s", safe(mu))},
          x);
    };
    write("g", "g = log(f)", g);
    write("f", "f = f(z) = 1/phi(1/z)", f);
  }
}

template<class T> Series<T> read_bottcher(const string& input) {
  auto [comments, g] = read_series<T>(input);
  print("reading from '%s':", input);
  for (const auto& c : comments)
    print("  %s", c);
  slow_assert(comments.size() && comments[0] == "g = log(f)", "bad comments: %s", comments);
  return move(g);
}

#define SERIES(T) \
  template Series<T> bottcher_base(); \
  template tuple<Series<T>,Undevice<T>> bottcher_step(Series<T>& g, const double tol); \
  template void areas<T>(const int max_k, const double); \
  template void write_bottcher(const string&, const string&, const Undevice<T>, \
                               SeriesView<const T>, SeriesView<const T> g); \
  template Series<T> read_bottcher(const string&);
#define AREAS(S) \
  template S area(SeriesView<const S> f); \
  SERIES(S)
AREAS(double)
AREAS(Expansion<2>)
IF_CUDA(
  SERIES(Device<double>)
  SERIES(Device<Expansion<2>>)
)

}  // namespace mandelbrot
