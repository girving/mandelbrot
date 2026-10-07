// Area of a filled Julia set K(c), |c| < 1/4, from the area transfer operator (see julia.h)
//   julia_area --c -0.2 0 --nr 40 --nt 121 --prec td

#include "julia.h"
#include "arb_cc.h"
#include "expansion_arith.h"
#include "print.h"
#include <argparse.hpp>
#include <flint/arb.h>
using namespace mandelbrot;

template<class S> static string digits(const S x) {
  Arb a;
  if constexpr (is_same_v<S, double>) arb_set_d(a, x);
  else a = exact_arb(x);
  char* s = arb_get_str(a, is_same_v<S, double> ? 17 : is_same_v<S, Expansion<2>> ? 33 : 49, ARB_STR_NO_RADIUS);
  const string r(s);
  flint_free(s);
  return r;
}

template<class S> static void run(const JuliaParams& p) {
  const auto r = julia_area<S>(p);
  print("c = %.17g + %.17gi, nr %d, nt %d: area %s", p.cx, p.cy, p.nr, p.nt, digits(r.area));
  print("  residual %.3g, %d refinements, %d GMRES iterations, %.2f s", r.residual, r.refinements, r.gmres_iters,
        r.secs);
  if (p.eig) print("  leading eigenvalue of L: %.12f", r.rho);
}

int main(int argc, char** argv) {
  argparse::ArgumentParser program("julia_area");
  program.add_argument("--c").help("c = x y").nargs(2).scan<'g', double>().default_value(vector<double>{0, 0});
  program.add_argument("--r1").help("inner radius (0: max(2|c|, 1/4))").scan<'g', double>().default_value(0.0);
  program.add_argument("--r2").help("outer radius").scan<'g', double>().default_value(1.5);
  program.add_argument("--nr").help("Chebyshev points in log r").scan<'i', int>().default_value(32);
  program.add_argument("--nt").help("angles (odd)").scan<'i', int>().default_value(101);
  program.add_argument("--grade").help("angular grading toward θ = 0, in [0, 1)").scan<'g', double>().default_value(0.0);
  program.add_argument("--pre").help("parabolic preconditioner: 1 on, 0 off, -1 auto").scan<'i', int>().default_value(-1);
  program.add_argument("--stencil").help("preconditioner interpolation width").scan<'i', int>().default_value(6);
  program.add_argument("--oversample").help("preconditioner oversampling").scan<'i', int>().default_value(2);
  program.add_argument("--nq").help("quadrature angles (0: 2 nt + 1)").scan<'i', int>().default_value(0);
  program.add_argument("--ng").help("Gauss-Legendre points (0: nr + 16)").scan<'i', int>().default_value(0);
  program.add_argument("--prec").help("double, dd, or td").default_value(string("double"));
  program.add_argument("--eig").help("estimate L's leading eigenvalue").default_value(false).implicit_value(true);
  program.add_argument("--verbose").default_value(false).implicit_value(true);
  program.parse_args(argc, argv);
  JuliaParams p;
  const auto c = program.get<vector<double>>("--c");
  p.cx = c[0]; p.cy = c[1];
  p.r1 = program.get<double>("--r1"); p.r2 = program.get<double>("--r2");
  p.nr = program.get<int>("--nr"); p.nt = program.get<int>("--nt");
  p.grade = program.get<double>("--grade");
  p.pre = program.get<int>("--pre");
  p.stencil = program.get<int>("--stencil");
  p.oversample = program.get<int>("--oversample");
  p.nq = program.get<int>("--nq"); p.ng = program.get<int>("--ng");
  p.verbose = program.get<bool>("--verbose");
  p.eig = program.get<bool>("--eig");
  const auto prec = program.get<string>("--prec");
  if (prec == "double") run<double>(p);
  else if (prec == "dd") run<Expansion<2>>(p);
  else if (prec == "td") run<Expansion<3>>(p);
  else die("unknown --prec %s", prec);
}
