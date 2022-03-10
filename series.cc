// Series arithmetic

#include "series.h"
#include "expansion_arith.h"
#include "poly.h"
#include <fstream>
namespace mandelbrot {

using std::endl;
using std::exception;
using std::ifstream;
using std::is_same_v;
using std::ofstream;
using std::runtime_error;
using std::string_view;

Series<double> approx(const Poly& x, const int64_t n) {
  slow_assert(n >= x.length());
  Series<double> y(n);
  y.set_counts(n, n);
  for (int64_t i = 0; i < n; i++)
    y[i] = arf_get_d(arb_midref(x[i]), ARF_RND_NEAR);
  return y;
}

double error(SeriesView<const double> x, SeriesView<const double> y, const bool relative) {
  if (x.known() != y.known())
    return numeric_limits<double>::infinity();
  double e = 0;
  const auto both = min(x.nonzero(), y.nonzero());
  for (int64_t i = 0; i < both; i++) {
    double d = abs(x[i] - y[i]);
    if (relative)
      d /= max(1., abs(y[i]));
    e = max(e, d);
  }
  for (int64_t i = both; i < x.nonzero(); i++)
    e = max(e, abs(x[i]));
  for (int64_t i = both; i < y.nonzero(); i++) {
    double d = abs(y[i]);
    if (relative)
      d = min(d, 1.);
    e = max(e, d);
  }
  return e;
}

double error(SeriesView<const double> x, initializer_list<double>&& ys, const bool relative) {
  return error(x, Series<double>(move(ys)), relative);
}

double error(SeriesView<const double> x, const Poly& y, const bool relative) {
  if (x.known() < y.length())
    return numeric_limits<double>::infinity();
  return error(x, approx(y, x.known()), relative);
}

DEF_SERIAL(add_scalar_kernel, (S* ys, const S a), ys[0] += a;)

template<class T> void add_scalar(Series<T>& x, const typename Series<T>::Scalar a) {
  slow_assert(x.nonzero());
  add_scalar_kernel(x.data(), a);
}

DEF_LOOP(high_addsub_loop, n, i, (S* y, const S* x, const int ynz, const int xnz, const int sign, const int s),
  const auto yi = i < ynz ? y[i] : S(0);
  auto xi = uint32_t(i-s) < uint32_t(xnz) ? x[i-s] : S(0);
  if (sign < 0) xi = -xi;
  y[i] = yi + xi;)

template<class T> void high_addsub(Series<T>& y, const int sign, const int64_t s, SeriesView<add_const_t<T>> x) {
  const auto ynz = y.nonzero(), xnz = x.nonzero();
  const auto nk = min(y.known(), x.known() + s);
  const auto nz = min(nk, max(ynz, xnz ? xnz + s : 0));
  slow_assert(abs(sign) == 1 && !y.alias(x) && nz <= y.limit());
  high_addsub_loop(nz, y.data(), x.data(), ynz, xnz, sign, s);
  y.set_counts(nk, nz);
}

DEF_LOOP(mul1p_post_loop, post, i, (S* z, const S* x, const int s, const int xnz),
  z[i] = (i < s ? S(0) : z[i]) + (i < xnz ? x[i] : S(0));)

template<class T> void mul1p_post(Series<T>& z, SeriesView<add_const_t<T>> x,
                                  const int64_t post, const int64_t s, const int64_t xnz) {
  mul1p_post_loop(post, z.data(), x.data(), s, xnz);
}

static string time_str() {
  // From https://stackoverflow.com/questions/16357999/current-date-and-time-as-string/16358264
  const auto t = std::time(nullptr);
  const auto tm = *std::localtime(&t);
  return format("time = %s", std::put_time(&tm, "%F %T"));
}

template<class T> void write_series(const string& path, const vector<string>& comments, SeriesView<const T> x) {
  ofstream out(path);

  // Add comments for known(), nonzero(), and time
  auto cs = comments;
  cs.push_back(format("known = %d", x.known()));
  cs.push_back(format("nonzero = %d", x.nonzero()));
  cs.push_back(time_str());

  // Write comments at top of file
  for (const auto& c : cs)
    out << "# " << c << endl;

  // Write series terms in plain text
  string s;
  const auto& hx = host_copy(x);
  const auto reduce = [](string& y, const string& x) { y += x; };
  const auto map = [&hx](const int64_t i) { auto s = safe(hx[i]); s += '\n'; return s; };
  map_reduce(s, reduce, map, x.nonzero());
  out << s;
}

template<class T> tuple<vector<string>,Series<T>> read_series(const string& path) {
  if constexpr (is_device<T>) {
    const auto [cs, x] = read_series<Undevice<T>>(path);
    Series<T> dx(x.nonzero());
    host_to_device(dx, x);
    return make_tuple(cs, move(dx));
  } else {
    const auto trim = [](string_view& v) {
      while (v.size() && isspace(v[0]))
        v.remove_prefix(1);
    };
    const auto number = [](const string_view c, const string_view n) {
      try {
        const int64_t i = stol(string(n));
        if (i >= 0) return i;
      } catch (const exception&) {}
      throw runtime_error(format("failed to parse comment '%s' as nonnegative integer", c));
    };
    ifstream in(path);
    int64_t known = -1, nonzero = -1, terms = -1;
    vector<string> comments;
    Series<T> x;
    string line;
    while (getline(in, line)) {
      string_view v(line);
      const auto c = v;
      trim(v);
      if (!v.size()) continue;  // Skip blank lines
      if (v[0] == '#') {  // Comment!
        v.remove_prefix(1);
        trim(v);
        comments.emplace_back(v);
        if (v.starts_with("known =")) {
          slow_assert(known < 0);
          v.remove_prefix(7);
          known = number(c, v);
        } else if (v.starts_with("nonzero =")) {
          slow_assert(nonzero < 0);
          v.remove_prefix(9);
          nonzero = number(c, v);
        } else if (v.starts_with("terms =")) {
          // Older files had 'terms' instead of 'nonzero' and 'known'
          slow_assert(known < 0 && nonzero < 0);
          v.remove_prefix(7);
          known = nonzero = number(c, v);
        }
        if (known >= 0 && nonzero >= 0) {
          slow_assert(known >= nonzero);
          Series<T>(nonzero).swap(x);
          x.set_counts(known, nonzero);
          terms = 0;
        }
      } else {  // Number, hopefully
        slow_assert(terms >= 0);
        slow_assert(terms < nonzero);
        T a;
        if constexpr (is_same_v<T,double>) a = stod(string(v));
        else a = T(v);
        x[terms++] = a;
      }
    }
    slow_assert(terms == nonzero);
    return make_tuple(comments, move(x));
  }
}

#define Ss(S) \
  template void add_scalar(Series<S>&, const typename Series<S>::Scalar); \
  template void high_addsub(Series<S>&, const int, const int64_t, SeriesView<const S>); \
  template void mul1p_post(Series<S>&, SeriesView<const S>, const int64_t, const int64_t, const int64_t); \
  template void write_series(const string& path, const vector<string>& comments, SeriesView<const S> x); \
  template tuple<vector<string>,Series<S>> read_series(const string& path);
Ss(double)
Ss(Expansion<2>)
IF_CUDA(
  Ss(Device<double>)
  Ss(Device<Expansion<2>>)
)

}  // namespace mandelbrot
