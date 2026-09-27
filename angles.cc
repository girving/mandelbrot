// External angles

#include "angles.h"
#include "debug.h"
#include <numeric>
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

vector<Root> lavaurs(const int max_period) {
  slow_assert(max_period <= 24, "lavaurs: max_period %d too large", max_period);
  // Exact comparison of a/(2^p-1) and b/(2^q-1)
  const auto less = [](const Periodic& a, const Periodic& b) {
    return a.k * ((uint64_t(1) << b.q) - 1) < b.k * ((uint64_t(1) << a.q) - 1);
  };
  // Chord (c,d) crosses chord (a,b) iff exactly one of c,d lies strictly between a and b
  const auto inside = [&](const Periodic& x, const Wake& w) { return less(w.lo, x) && less(x, w.hi); };
  vector<Root> roots;
  for (int p = 2; p <= max_period; p++) {
    const uint64_t M = (uint64_t(1) << p) - 1;
    // Angles of exact period p, in increasing order
    vector<Periodic> angles;
    for (uint64_t k = 1; k < M; k++) {
      bool exact = true;
      for (int d = 1; d < p && exact; d++)
        if (p % d == 0 && k % (M / ((uint64_t(1) << d) - 1)) == 0) exact = false;
      if (exact) angles.push_back(Periodic{k, p});
    }
    vector<bool> paired(angles.size());
    for (size_t i = 0; i < angles.size(); i++) {
      if (paired[i]) continue;
      bool done = false;
      for (size_t j = i + 1; j < angles.size() && !done; j++) {
        if (paired[j]) continue;
        const Wake w{angles[i], angles[j]};
        bool crosses = false;
        for (const auto& r : roots)
          if (inside(w.lo, r.w) != inside(w.hi, r.w)) { crosses = true; break; }
        if (crosses) continue;
        // Satellite iff hi is on lo's doubling cycle
        bool same = false;
        uint64_t x = w.lo.k;
        for (int s = 0; s < p; s++, x = 2*x % M)
          if (x == w.hi.k) same = true;
        roots.push_back(Root{w, same});
        paired[i] = paired[j] = true;
        done = true;
      }
      slow_assert(done, "lavaurs: failed to pair %s at period %d", angles[i], p);
    }
  }
  return roots;
}

}  // namespace mandelbrot
