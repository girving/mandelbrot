// Renormalization structure in the Böttcher energy map over external angle
//
// Input: angle_map.npy of shape [octaves, bins], where row j is the energy of Böttcher octave j
// (sum over 2^j <= m < 2^(j+1) of (m-1) f_m^2) as a density over half-circle angle bins, folded
// under θ -> 1-θ (b_n is real, so the full-circle density is symmetric).
//
// A. Sub-wake tracking: does the p/q sub-wake of a tuned component (e.g. the period 2 disk) carry energy
//    proportional to the cardioid's p/q wake at a fixed lag in octaves?
// B. Pullback self-similarity: map a wake with root period r onto the circle by θ -> 2^r θ, and find the
//    lag L for which its angular distribution at octave j best matches the full distribution at j - L.

#include "angles.h"
#include "debug.h"
#include "numpy.h"
#include "print.h"
#include "wall_time.h"
#include <bit>
#include <cmath>
#include <functional>
#include <numeric>
#include <vector>
namespace mandelbrot {
namespace {

using std::function;
using std::gcd;
using std::vector;

struct Map {
  int octaves, bins;         // bins covers the full circle
  vector<vector<double>> e;  // e[j][i] = energy of octave j in θ ∈ [i/bins, (i+1)/bins)
};

Map read_map(const string& path) {
  const auto a = read_numpy(path);
  slow_assert(a.shape.size() == 2, "expected 2D angle map");
  Map m{int(a.shape[0]), 2 * int(a.shape[1]), {}};
  const int half = int(a.shape[1]);
  for (int j = 0; j < m.octaves; j++) {
    vector<double> row(m.bins);
    for (int i = 0; i < m.bins; i++)
      row[i] = 0.5 * a.data[int64_t(j) * half + (i < half ? i : m.bins - 1 - i)];
    m.e.push_back(move(row));
  }
  return m;
}

double center(const Map& m, const int i) { return (i + 0.5) / m.bins; }

// Energy per octave inside any of the given wakes
vector<double> energy(const Map& m, const vector<Wake>& ws) {
  vector<double> E(m.octaves);
  for (int i = 0; i < m.bins; i++) {
    const double t = center(m, i);
    for (const auto& w : ws)
      if (w.contains(t)) {
        for (int j = 0; j < m.octaves; j++) E[j] += m.e[j][i];
        break;
      }
  }
  return E;
}

// Fraction of variance of detrended log2 E over [j0,j1] explained by a period-P pattern
double periodic_fit(const vector<double>& E, const int j0, const int j1, const int P) {
  const int n = j1 - j0 + 1;
  vector<double> y(n);
  double sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (int k = 0; k < n; k++) {
    y[k] = std::log2(E[j0 + k]);
    const double x = j0 + k;
    sx += x; sy += y[k]; sxx += x*x; sxy += x*y[k];
  }
  const double b = (n*sxy - sx*sy) / (n*sxx - sx*sx), a = (sy - b*sx) / n;
  vector<double> r(n), mean(P), count(P);
  double var = 0;
  for (int k = 0; k < n; k++) {
    r[k] = y[k] - a - b*(j0 + k);
    var += r[k]*r[k];
    mean[(j0 + k) % P] += r[k];
    count[(j0 + k) % P]++;
  }
  double left = 0;
  for (int k = 0; k < n; k++) {
    const int c = (j0 + k) % P;
    const double d = r[k] - mean[c] / count[c];
    left += d*d;
  }
  return 1 - left / var;
}

int best_period(const vector<double>& E, const int j0, const int j1, double& explained) {
  int best = 0;
  explained = -1;
  for (int P = 2; P <= 10; P++) {
    const double f = periodic_fit(E, j0, j1, P);
    if (f > explained + 0.02) { explained = f; best = P; }  // Prefer the smallest period that does as well
  }
  return best;
}

// Std of log2(X(j) / Y(j - s)) over j in [j0, j1]
double log_ratio_std(const vector<double>& X, const vector<double>& Y, const int s, const int j0, const int j1) {
  double s1 = 0, s2 = 0;
  int n = 0;
  for (int j = j0; j <= j1; j++) {
    if (j - s < 1) continue;
    const double r = std::log2(X[j] / Y[j - s]);
    s1 += r; s2 += r*r; n++;
  }
  const double mean = s1 / n;
  return std::sqrt(std::max(0.0, s2 / n - mean*mean));
}

// Angular distribution of octave j pulled back from wakes ws by θ -> 2^r θ, on K full-circle bins
vector<double> pullback(const Map& m, const vector<Wake>& ws, const int r, const int j, const int K) {
  vector<double> P(K);
  double total = 0;
  for (int i = 0; i < m.bins; i++) {
    const double t = center(m, i);
    for (const auto& w : ws)
      if (w.contains(t)) {
        const double u = std::ldexp(t, r);
        P[int((u - std::floor(u)) * K) % K] += m.e[j][i];
        total += m.e[j][i];
        break;
      }
  }
  for (auto& p : P) p /= total;
  return P;
}

vector<double> full(const Map& m, const int j, const int K) {
  return pullback(m, {Wake{{0, 1}, {1, 1}}}, 0, j, K);
}

double tv(const vector<double>& P, const vector<double>& Q) {
  double d = 0;
  for (size_t i = 0; i < P.size(); i++) d += std::abs(P[i] - Q[i]);
  return d / 2;
}

// Angular distribution (2^D bins) of octave j's energy on angles that untune to at least D digits,
// and the fraction of the wake's energy that does
vector<double> untuned(const Map& m, const Wake& w, const int j, const int D, double& captured) {
  vector<double> P(1 << D);
  double total = 0, wake = 0;
  for (int i = 0; i < m.bins; i++) {
    if (!w.contains(center(m, i))) continue;
    wake += m.e[j][i];
    uint64_t prefix;
    if (untune(w, uint64_t(i), std::countr_zero(unsigned(m.bins)), D, prefix) == D) {
      P[prefix] += m.e[j][i];
      total += m.e[j][i];
    }
  }
  captured = total / wake;
  for (auto& p : P) p /= total;
  return P;
}

vector<Wake> cardioid_wakes(const int q) {
  vector<Wake> ws;
  for (int p = 1; p < q; p++)
    if (gcd(p, q) == 1) ws.push_back(cardioid_wake(p, q));
  return ws;
}

vector<Wake> tuned(const Wake& w, const vector<Wake>& vs) {
  vector<Wake> r;
  for (const auto& v : vs) r.push_back(tune(w, v));
  return r;
}

void run(const string& path) {
  const auto t0 = wall_time();
  const auto m = read_map(path);
  print("read %s: %d octaves, %d full-circle bins (%.3f s)", path, m.octaves, m.bins, (wall_time() - t0).seconds());
  const int j0 = 12, j1 = m.octaves - 1;
  const auto half = cardioid_wake(1, 2);

  // A. Periodicity of each wake class, and sub-wake tracking
  print("\nA. Octave periodicity of wake energy over j = %d..%d (detrended log2 energy)", j0, j1);
  print("   class                               period  explained");
  vector<vector<double>> C(7), D(6);
  for (int q = 2; q <= 6; q++) {
    C[q] = energy(m, cardioid_wakes(q));
    double f;
    const int P = best_period(C[q], j0, j1, f);
    print("   cardioid p/q wakes, q = %d            %6d  %8.0f%%", q, P, 100*f);
  }
  for (int q = 2; q <= 5; q++) {
    D[q] = energy(m, tuned(half, cardioid_wakes(q)));
    double f;
    const int P = best_period(D[q], j0, j1, f);
    print("   period-2 disk p/q sub-wakes, q = %d   %6d  %8.0f%%", q, P, 100*f);
  }
  print("\n   Tracking: std of log2(sub-wake_q(j) / cardioid-wake_q(j - s)) over j = %d..%d; best lags marked", j0 + 2, j1);
  string header = "   q   ";
  for (int s = 0; s <= 10; s++) header += tfm::format("  s=%-3d", s);
  print(header);
  for (int q = 2; q <= 5; q++) {
    string row = tfm::format("   %d   ", q);
    double best = 1e9;
    int bs = 0;
    vector<double> sd;
    for (int s = 0; s <= 10; s++) {
      sd.push_back(log_ratio_std(D[q], C[q], s, j0 + 2, j1));
      if (sd.back() < best) { best = sd.back(); bs = s; }
    }
    for (int s = 0; s <= 10; s++) row += tfm::format(" %5.2f%s", sd[s], s == bs ? "*" : " ");
    print(row);
  }
  print("   (std 0 = exact proportionality; for comparison, the raw log2 std of each class is ~0.3-0.6)");

  // B. Pullback self-similarity
  const int K = 256;
  print("\nB. Pullback θ -> 2^r θ from wakes, vs the full distribution at octave j - L: mean TV distance (%%), j = 16..%d, K = %d", j1, K);
  struct Case { string name; vector<Wake> ws; int r; };
  vector<Case> cases;
  for (int q = 2; q <= 5; q++) {
    for (int p = 1; p < q; p++)
      if (gcd(p, q) == 1 && 2*p <= q)  // By symmetry p and q-p give mirrored results
        cases.push_back({tfm::format("cardioid %d/%d wake (r=%d)", p, q, q), {cardioid_wake(p, q)}, q});
  }
  for (int q = 2; q <= 3; q++)
    cases.push_back({tfm::format("period-2 disk 1/%d sub-wake (r=%d)", q, 2*q), {tune(half, cardioid_wake(1, q))}, 2*q});
  header = tfm::format("   %-34s", "wake");
  for (int L = 0; L <= 10; L++) header += tfm::format(" L=%-3d", L);
  print(header);
  for (const auto& c : cases) {
    vector<double> mean(11);
    for (int L = 0; L <= 10; L++) {
      int n = 0;
      for (int j = 16; j <= j1; j++) {
        if (j - L < 1) continue;
        mean[L] += tv(pullback(m, c.ws, c.r, j, K), full(m, j - L, K));
        n++;
      }
      mean[L] /= n;
    }
    int bL = 0;
    for (int L = 1; L <= 10; L++) if (mean[L] < mean[bL]) bL = L;
    string row = tfm::format("   %-34s", c.name);
    for (int L = 0; L <= 10; L++) row += tfm::format(" %4.0f%s", 100*mean[L], L == bL ? "*" : " ");
    print(row);
  }
  // Baseline: how different are full distributions at different octaves at all?
  {
    string row = tfm::format("   %-34s", "baseline: full(j) vs full(j-L)");
    for (int L = 0; L <= 10; L++) {
      double s = 0; int n = 0;
      for (int j = 16; j <= j1; j++) if (j - L >= 1) { s += tv(full(m, j, K), full(m, j - L, K)); n++; }
      row += tfm::format(" %4.0f ", 100*s/n);
    }
    print(row);
  }
  // C. Inverse-tuning pullback: angles near the tuned copy of M (words of the tuning wake) decoded to M's angles
  print("\nC. Inverse tuning: copy's octave j vs the full distribution at octave j' (TV distance %%, 2^D bins)");
  struct Tuning { string name; Wake w; int D; };
  for (const auto& T : {Tuning{"period-2 copy (0->01, 1->10)", half, 7}, Tuning{"period-3 copy (1/3 bulb)", cardioid_wake(1, 3), 5}}) {
    print("   %s, D = %d: best j' for each j, TV there, TV at j' = j, and energy captured", T.name, T.D);
    print("     j   best j'   TV(best)   TV(j'=j)   j/j'   captured");
    for (int j = 12; j <= j1; j++) {
      double captured;
      const auto P = untuned(m, T.w, j, T.D, captured);
      int bj = 1;
      double best = 2;
      for (int jp = 1; jp <= j; jp++) {
        const double d = tv(P, full(m, jp, 1 << T.D));
        if (d < best) { best = d; bj = jp; }
      }
      print("    %2d   %6d   %7.0f%%   %7.0f%%   %5.2f   %6.1f%%", j, bj, 100*best, 100*tv(P, full(m, j, 1 << T.D)),
            double(j) / bj, 100*captured);
    }
  }
  print("\ntotal %.3f s", (wall_time() - t0).seconds());
}

}  // namespace
}  // namespace mandelbrot

int main(const int argc, const char** argv) {
  using namespace mandelbrot;
  try {
    slow_assert(argc == 2, "usage: %s <angle_map.npy>", argv[0]);
    run(argv[1]);
    return 0;
  } catch (const std::exception& e) {
    die(e.what());
  }
}
