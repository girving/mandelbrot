// Series arithmetic

#include "series.h"
#include "expansion_arith.h"
#include "numpy.h"
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

DEF_LOOP(high_addsub_ldexp_loop, n, i,
         (S* y, const S* x, const int ynz, const int xnz, const int sign, const int s, const int b),
  const auto yi = i < ynz ? y[i] : S(0);
  auto xi = ldexp(uint32_t(i-s) < uint32_t(xnz) ? x[i-s] : S(0), b);
  if (sign < 0) xi = -xi;
  y[i] = yi + xi;)

template<class T> void high_addsub_ldexp(Series<T>& y, const int sign, const int b, const int64_t s,
                                         SeriesView<add_const_t<T>> x) {
  const auto ynz = y.nonzero(), xnz = x.nonzero();
  const auto nk = min(y.known(), x.known() + s);
  const auto nz = min(nk, max(ynz, xnz ? xnz + s : 0));
  slow_assert(abs(sign) == 1 && !y.alias(x) && nz <= y.limit());
  if (!b)
    high_addsub_loop(nz, y.data(), x.data(), ynz, xnz, sign, s);
  else
    high_addsub_ldexp_loop(nz, y.data(), x.data(), ynz, xnz, sign, s, b);
  y.set_counts(nk, nz);
}

DEF_LOOP(mul1p_middle_loop, n, i, (S* z, const S* x, const int nx),
  z[i] = i < nx ? x[i] : 0;)

template<class T> void mul1p_middle(Series<T>& z, const T* x, const int64_t nx) {
  mul1p_middle_loop(z.nonzero(), z.data(), x, nx);
}

static string time_str() {
  // From https://stackoverflow.com/questions/16357999/current-date-and-time-as-string/16358264
  const auto t = std::time(nullptr);
  const auto tm = *std::localtime(&t);
  return tfm::format("time = %s", std::put_time(&tm, "%F %T"));
}

template<class T> void write_series(const string& path, const vector<string>& comments, SeriesView<const T> x,
                                    const int64_t batch_size) {
  ofstream out(path);

  // Add comments for known(), nonzero(), and time
  auto cs = comments;
  cs.push_back(tfm::format("known = %d", x.known()));
  cs.push_back(tfm::format("nonzero = %d", x.nonzero()));
  cs.push_back(time_str());

  // Write comments at top of file
  for (const auto& c : cs)
    out << "# " << c << endl;

  // Write series terms in plain text, in batches
  const auto& hx = host_copy(x);
  for (int64_t start = 0; start < x.nonzero(); start += batch_size) {
    string s;
    const auto reduce = [](string& y, const string& x) { y += x; };
    const auto map = [&hx,start](const int64_t i) { auto s = safe(hx[start + i]); s += '\n'; return s; };
    map_reduce(s, reduce, map, min(batch_size, relu(x.nonzero() - start)));
    out << s;
  }
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
      throw runtime_error(tfm::format("failed to parse comment '%s' as nonnegative integer", c));
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
        slow_assert(terms >= 0 && terms < nonzero);
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

template<class T> void write_series_npy(const string& path, SeriesView<const T> x) {
  slow_assert(path.ends_with(".npy"));

  // Prepare shape
  vector<int> shape = {x.nonzero()};
  typedef Undevice<T> S;
  if constexpr (is_same_v<S,Expansion<2>>) shape.push_back(2);
  else { static_assert(is_same_v<S,double>); }
  const auto header = numpy_header(shape);

  // Copy to host
  const auto hx = host_copy(x);

  // Write file
  FILE* f = fopen(path.c_str(), "wb");
  slow_assert(f);
  fwrite(header.c_str(), header.size(), 1, f);
  fwrite(hx.data(), sizeof(S)*x.nonzero(), 1, f);
  fclose(f);
}

#define Ss(S) \
  template void add_scalar(Series<S>&, const typename Series<S>::Scalar); \
  template void high_addsub_ldexp(Series<S>&, const int, const int, const int64_t, SeriesView<const S>); \
  template void mul1p_middle(Series<S>&, const S*, const int64_t); \
  template void write_series(const string& path, const vector<string>& comments, SeriesView<const S> x, const int64_t); \
  template tuple<vector<string>,Series<S>> read_series(const string& path); \
  template void write_series_npy(const string&, SeriesView<const S>);
Ss(double)
Ss(Expansion<2>)
IF_CUDA(
  Ss(Device<double>)
  Ss(Device<Expansion<2>>)
)

}  // namespace mandelbrot
