// External angles

#include "angles.h"
#include "debug.h"
#include <algorithm>
#include <map>
#include <tuple>
#include <numeric>
#include <unordered_map>
#include <vector>
namespace mandelbrot {

using std::gcd;
using std::vector;

Wake cardioid_wake(const int p, const int q) {
  slow_assert(1 <= p && p < q && q <= 62 && gcd(p, q) == 1, "bad rotation number %d/%d", p, q);
  // The unique doubling q-cycle x_0 < ... < x_{q-1} with rotation number p/q maps x_i -> x_{i+p mod q}.
  // Exactly p of its points lie in [1/2,1), so the binary digit of x_i is [i >= q-p], and x_i's digits
  // are those of x_i, x_{i+p}, x_{i+2p}, ...
  vector<uint64_t> k(q);
  for (int i = 0; i < q; i++)
    for (int m = 0; m < q; m++)
      k[i] = k[i] << 1 | ((i + int64_t(m) * p) % q >= q - p);
  // The wake is the arc between cyclically adjacent points at distance 1/(2^q-1)
  for (int i = 0; i + 1 < q; i++)
    if (k[i+1] - k[i] == 1)
      return Wake{{k[i], q}, {k[i+1], q}};
  die("no wake found for %d/%d", p, q);
}

Periodic tune(const Wake& w, const Periodic& a) {
  const int r = w.lo.q;
  slow_assert(w.hi.q == r && int64_t(r) * a.q <= 62, "tune: periods %d, %d too large", r, a.q);
  uint64_t k = 0;
  for (int i = a.q - 1; i >= 0; i--)
    k = k << r | (a.k >> i & 1 ? w.hi.k : w.lo.k);
  return Periodic{k, r * a.q};
}

Wake tune(const Wake& w, const Wake& v) {
  return Wake{tune(w, v.lo), tune(w, v.hi)};
}

int untune(const Wake& w, const uint64_t bits, const int nbits, const int D, uint64_t& prefix) {
  const int r = w.lo.q;
  slow_assert(w.hi.q == r && nbits <= 64, "untune: bad arguments");
  prefix = 0;
  int d = 0;
  for (; d < D && (d + 1) * r <= nbits; d++) {
    const uint64_t word = bits >> (nbits - (d + 1) * r) & ((uint64_t(1) << r) - 1);
    if (word == w.lo.k) prefix = prefix << 1;
    else if (word == w.hi.k) prefix = prefix << 1 | 1;
    else break;
  }
  return d;
}

namespace {

// Point update, range max (or min, with negated values) over [lo, hi)
struct MaxTree {
  int64_t n;
  vector<int64_t> t;
  explicit MaxTree(const int64_t size) : n(1) { while (n < size) n *= 2; t.assign(2*n, -1); }
  void set(int64_t i, const int64_t v) { for (t[i += n] = v; i >>= 1;) t[i] = std::max(t[2*i], t[2*i+1]); }
  int64_t max(int64_t lo, int64_t hi) const {
    int64_t r = -1;
    for (lo += n, hi += n; lo < hi; lo >>= 1, hi >>= 1) {
      if (lo & 1) r = std::max(r, t[lo++]);
      if (hi & 1) r = std::max(r, t[--hi]);
    }
    return r;
  }
};

}  // namespace

vector<int> maximal_tuning(const vector<Root>& roots) {
  // Roots by (period, lo word, hi word)
  std::map<std::tuple<int, uint64_t, uint64_t>, int> index;
  const auto key = [](const int q, const uint64_t lo, const uint64_t hi) { return std::make_tuple(q, lo, hi); };
  for (size_t i = 0; i < roots.size(); i++) {
    const auto& w = roots[i].w;
    slow_assert(w.lo.q == w.hi.q && w.lo.q <= 62, "maximal_tuning: period %d too large", w.lo.q);
    index[key(w.lo.q, w.lo.k, w.hi.k)] = int(i);
  }
  vector<int> parent(roots.size(), -1);
  for (size_t i = 0; i < roots.size(); i++) {
    const auto& w = roots[i].w;
    const int p = w.lo.q;
    for (int d = 2; d < p && parent[i] < 0; d++) {
      if (p % d) continue;
      const uint64_t mask = (uint64_t(1) << d) - 1;
      // Distinct d-bit words of both angles
      uint64_t words[2];
      int n = 0;
      bool ok = true;
      for (const uint64_t k : {w.lo.k, w.hi.k})
        for (int b = 0; b < p && ok; b += d) {
          const uint64_t u = k >> b & mask;
          if (n > 0 && words[0] == u) continue;
          if (n > 1 && words[1] == u) continue;
          if (n == 2) ok = false;
          else words[n++] = u;
        }
      if (!ok || n != 2) continue;
      const auto it = index.find(key(d, std::min(words[0], words[1]), std::max(words[0], words[1])));
      if (it != index.end()) parent[i] = it->second;
    }
  }
  return parent;
}

namespace {

bool exact_period(const uint64_t k, const int p) {
  const uint64_t M = (uint64_t(1) << p) - 1;
  for (int d = 1; d < p; d++)
    if (p % d == 0 && k % (M / ((uint64_t(1) << d) - 1)) == 0) return false;
  return true;
}

// Lavaurs' pairing of the given angles (all of exact period ≤ max_period in an arc whose chords pair among
// themselves), sorted here
vector<Root> lavaurs_pair(vector<Periodic> all, const int max_period) {
  std::sort(all.begin(), all.end(), [](const Periodic& a, const Periodic& b) {
    typedef unsigned __int128 U;
    return U(a.k) * ((uint64_t(1) << b.q) - 1) < U(b.k) * ((uint64_t(1) << a.q) - 1);
  });
  const int64_t N = all.size();

  // left[i] = partner of a chord's left endpoint at i; right[i] = N - 1 - partner of a right endpoint at i
  MaxTree left(N), right(N);
  vector<Root> roots;
  for (int p = 2; p <= max_period; p++) {
    vector<int64_t> pos;  // Positions of period p angles, increasing
    for (int64_t i = 0; i < N; i++)
      if (all[i].q == p) pos.push_back(i);
    const int64_t m = pos.size();
    vector<int64_t> next(m + 1);  // Union-find to the next unpaired index
    std::iota(next.begin(), next.end(), 0);
    const auto find = [&](int64_t i) {
      while (next[i] != i) i = next[i] = next[next[i]];
      return i;
    };
    const uint64_t M = (uint64_t(1) << p) - 1;
    for (int64_t i = 0; i < m; i++) {
      if (find(i) != i) continue;
      const int64_t a = pos[i];
      int64_t j = find(i + 1);
      for (;;) {
        slow_assert(j < m, "lavaurs: no partner for %s", all[a]);
        const int64_t b = pos[j];
        // A right endpoint in (a,b) whose partner is before a would mean we skipped a's region
        slow_assert(N - 1 - right.max(a + 1, b) > a || right.max(a + 1, b) < 0, "lavaurs: skipped region");
        const int64_t far = left.max(a + 1, b);
        if (far > b) {  // A chord starts in (a,b) and ends beyond b: jump past its end
          j = find(std::lower_bound(pos.begin(), pos.end(), far) - pos.begin());
          continue;
        }
        const Wake w{all[a], all[b]};
        bool same = false;
        uint64_t x = w.lo.k;
        for (int s = 0; s < p; s++, x = 2*x % M)
          if (x == w.hi.k) same = true;
        roots.push_back(Root{w, same});
        left.set(a, b);
        right.set(b, N - 1 - a);
        next[i] = i + 1;
        next[j] = j + 1;
        break;
      }
    }
  }
  return roots;
}

}  // namespace

vector<Root> lavaurs(const int max_period) {
  slow_assert(max_period <= 24, "lavaurs: max_period %d too large", max_period);
  // All angles of exact period 2..max_period
  vector<Periodic> all;
  for (int p = 2; p <= max_period; p++)
    for (uint64_t k = 1; k < (uint64_t(1) << p) - 1; k++)
      if (exact_period(k, p)) all.push_back(Periodic{k, p});
  return lavaurs_pair(std::move(all), max_period);
}

vector<Root> lavaurs(const int max_period, const Wake& within) {
  const int q = within.lo.q;
  slow_assert(within.hi.q == q && max_period <= 62 && q <= max_period, "lavaurs: bad period %d or wake", max_period);
  // Angles of exact period q..max_period in [lo, hi]: k/M ≥ lo.k/Q iff k Q ≥ lo.k M (Q = 2^q - 1), in 128 bits
  typedef unsigned __int128 U;
  const uint64_t Q = (uint64_t(1) << q) - 1;
  vector<Periodic> all;
  for (int p = q; p <= max_period; p++) {
    const uint64_t M = (uint64_t(1) << p) - 1;
    for (uint64_t k = uint64_t((U(within.lo.k) * M + Q - 1) / Q); U(k) * Q <= U(within.hi.k) * M; k++)
      if (exact_period(k, p)) all.push_back(Periodic{k, p});
  }
  return lavaurs_pair(std::move(all), max_period);
}

vector<int32_t> wake_owner(const vector<Root>& roots, const int64_t n) {
  // Wakes starting in [0,1/2), sorted by start
  vector<int32_t> order;
  for (size_t r = 0; r < roots.size(); r++)
    if (roots[r].w.lo.value() < 0.5) order.push_back(int32_t(r));
  std::sort(order.begin(), order.end(), [&](const int32_t a, const int32_t b) {
    return roots[a].w.lo.value() < roots[b].w.lo.value();
  });
  vector<int32_t> owner(n/2 + 1, -1), stack;
  size_t next = 0;
  for (int64_t i = 0; i <= n/2; i++) {
    const double t = double(i) / double(n);
    // Open every wake starting before t, first closing those that ended before it starts
    for (; next < order.size() && roots[order[next]].w.lo.value() < t; next++) {
      const auto& w = roots[order[next]].w;
      while (stack.size() && roots[stack.back()].w.hi.value() <= w.lo.value()) stack.pop_back();
      if (w.hi.value() > t) stack.push_back(order[next]);
    }
    while (stack.size() && roots[stack.back()].w.hi.value() <= t) stack.pop_back();
    if (stack.size()) owner[i] = stack.back();
  }
  return owner;
}

}  // namespace mandelbrot
