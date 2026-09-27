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

}  // namespace mandelbrot
