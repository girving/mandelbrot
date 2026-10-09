// External angle tests

#include "angles.h"
#include "tests.h"
#include <algorithm>
#include <numeric>
#include <tuple>
#include <vector>
namespace mandelbrot {
namespace {

using std::gcd;
using std::sort;
using std::vector;

// Brute force: search all q-cycles of doubling on k/(2^q-1) for the one with rotation number p/q
Wake brute_wake(const int p, const int q) {
  const uint64_t M = (uint64_t(1) << q) - 1;
  for (uint64_t k = 1; k < M; k++) {
    vector<uint64_t> orbit = {k};
    for (int i = 1; i < q; i++) orbit.push_back(2 * orbit.back() % M);
    if (*std::min_element(orbit.begin(), orbit.end()) != k || 2 * orbit.back() % M != k) continue;
    vector<uint64_t> s = orbit;
    sort(s.begin(), s.end());
    if (std::adjacent_find(s.begin(), s.end()) != s.end()) continue;  // Not exact period q
    bool rotation = true;
    for (int i = 0; i < q; i++) {
      const auto j = std::find(s.begin(), s.end(), 2 * s[i] % M) - s.begin();
      if ((j - i + q) % q != p) rotation = false;
    }
    if (!rotation) continue;
    for (int i = 0; i + 1 < q; i++)
      if (s[i+1] - s[i] == 1) return Wake{{s[i], q}, {s[i+1], q}};
  }
  die("brute_wake failed");
}

TEST(cardioid_wake_known) {
  ASSERT_EQ(cardioid_wake(1, 2), (Wake{{1, 2}, {2, 2}}));  // 1/3, 2/3: root -3/4
  ASSERT_EQ(cardioid_wake(1, 3), (Wake{{1, 3}, {2, 3}}));  // 1/7, 2/7
  ASSERT_EQ(cardioid_wake(2, 3), (Wake{{5, 3}, {6, 3}}));  // 5/7, 6/7
  ASSERT_EQ(cardioid_wake(1, 4), (Wake{{1, 4}, {2, 4}}));  // 1/15, 2/15
  ASSERT_EQ(cardioid_wake(2, 5), (Wake{{9, 5}, {10, 5}}));  // 9/31, 10/31
}

TEST(cardioid_wake_brute) {
  for (int q = 2; q <= 12; q++)
    for (int p = 1; p < q; p++)
      if (gcd(p, q) == 1)
        ASSERT_EQ(cardioid_wake(p, q), brute_wake(p, q)) << tfm::format("p/q = %d/%d", p, q);
}

TEST(cardioid_wakes_disjoint) {
  // Wakes are disjoint, symmetric under θ -> 1-θ, and their lengths sum to nearly 1 (Douady)
  vector<Wake> ws;
  double total = 0;
  for (int q = 2; q <= 20; q++)
    for (int p = 1; p < q; p++)
      if (gcd(p, q) == 1) {
        const auto w = cardioid_wake(p, q), v = cardioid_wake(q - p, q);
        const uint64_t M = (uint64_t(1) << q) - 1;
        ASSERT_EQ(v, (Wake{{M - w.hi.k, q}, {M - w.lo.k, q}}));
        ASSERT_EQ(w.hi.k - w.lo.k, uint64_t(1));
        ws.push_back(w);
        total += w.length();
      }
  sort(ws.begin(), ws.end(), [](const Wake& a, const Wake& b) { return a.lo.value() < b.lo.value(); });
  for (size_t i = 0; i + 1 < ws.size(); i++)
    ASSERT_LT(ws[i].hi.value(), ws[i+1].lo.value());
  ASSERT_LT(1 - 1e-4, total);
  ASSERT_LT(total, 1);
}

TEST(tune) {
  const auto half = cardioid_wake(1, 2);
  // Period doubling: the 1/2 bulb of the period 2 disk is the period 4 disk with angles 2/5 = 6/15, 3/5 = 9/15
  ASSERT_EQ(tune(half, half), (Wake{{6, 4}, {9, 4}}));
  // The 1/3 bulb of the period 2 disk has angles 22/63 and 25/63
  ASSERT_EQ(tune(half, cardioid_wake(1, 3)), (Wake{{22, 6}, {25, 6}}));
  // Tuning by the main cardioid's wake (0 -> 0, 1 -> 1 with 1-bit words) is the identity
  const Wake identity{{0, 1}, {1, 1}};
  for (int q = 2; q <= 8; q++)
    for (int p = 1; p < q; p++)
      if (gcd(p, q) == 1)
        ASSERT_EQ(tune(identity, cardioid_wake(p, q)), cardioid_wake(p, q));
  // Tuned wakes lie inside the tuning wake and are disjoint
  vector<Wake> ws;
  for (int q = 2; q <= 10; q++)
    for (int p = 1; p < q; p++)
      if (gcd(p, q) == 1) {
        const auto w = tune(half, cardioid_wake(p, q));
        ASSERT_TRUE(half.lo.value() < w.lo.value() && w.hi.value() < half.hi.value()) << w;
        ws.push_back(w);
      }
  sort(ws.begin(), ws.end(), [](const Wake& a, const Wake& b) { return a.lo.value() < b.lo.value(); });
  for (size_t i = 0; i + 1 < ws.size(); i++)
    ASSERT_LT(ws[i].hi.value(), ws[i+1].lo.value());
}

TEST(untune) {
  const auto half = cardioid_wake(1, 2), third = cardioid_wake(1, 3);
  for (const auto& w : {half, third}) {
    const int r = w.lo.q;
    for (int q = 1; q <= 8; q++)
      for (uint64_t k = 0; k < (uint64_t(1) << q); k++) {
        // tune(w, k/(2^q-1)) repeats a word of r*q bits; its first r*q bits decode back to k's q bits
        const auto t = tune(w, Periodic{k, q});
        uint64_t prefix;
        ASSERT_EQ(untune(w, t.k, r * q, q, prefix), q);
        ASSERT_EQ(prefix, k);
        // Asking for fewer digits stops early
        ASSERT_EQ(untune(w, t.k, r * q, q - 1, prefix), q - 1);
        ASSERT_EQ(prefix, k >> 1);
      }
  }
  // Decoding stops at the first word outside {01, 10}
  uint64_t prefix;
  ASSERT_EQ(untune(half, 0b01101101, 8, 4, prefix), 2);  // 01 10 11 01
  ASSERT_EQ(prefix, uint64_t(0b01));
  ASSERT_EQ(untune(half, 0b00, 2, 4, prefix), 0);
  // Too few bits for another word
  ASSERT_EQ(untune(half, 0b011, 3, 4, prefix), 1);
}

// Reference Lavaurs: check each candidate chord against every earlier chord
vector<Root> slow_lavaurs(const int max_period) {
  const auto less = [](const Periodic& a, const Periodic& b) {
    return a.k * ((uint64_t(1) << b.q) - 1) < b.k * ((uint64_t(1) << a.q) - 1);
  };
  const auto inside = [&](const Periodic& x, const Wake& w) { return less(w.lo, x) && less(x, w.hi); };
  vector<Root> roots;
  for (int p = 2; p <= max_period; p++) {
    const uint64_t M = (uint64_t(1) << p) - 1;
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
      for (size_t j = i + 1; j < angles.size(); j++) {
        if (paired[j]) continue;
        const Wake w{angles[i], angles[j]};
        bool crosses = false;
        for (const auto& r : roots)
          if (inside(w.lo, r.w) != inside(w.hi, r.w)) { crosses = true; break; }
        if (crosses) continue;
        bool same = false;
        uint64_t x = w.lo.k;
        for (int s = 0; s < p; s++, x = 2*x % M)
          if (x == w.hi.k) same = true;
        roots.push_back(Root{w, same});
        paired[i] = paired[j] = true;
        break;
      }
    }
  }
  return roots;
}

TEST(lavaurs_fast_vs_slow) {
  const auto fast = lavaurs(11), slow = slow_lavaurs(11);
  ASSERT_EQ(fast.size(), slow.size());
  for (size_t i = 0; i < fast.size(); i++) {
    ASSERT_EQ(fast[i].w, slow[i].w) << i;
    ASSERT_EQ(fast[i].satellite, slow[i].satellite);
  }
}

TEST(lavaurs) {
  const int P = 18;
  const auto roots = lavaurs(P);
  // Component counts per period: OEIS A000740
  const int comps[] = {0, 1, 1, 3, 6, 15, 27, 63, 120, 252, 495, 1023, 2010, 4095, 8127, 16365, 32640, 65535, 130788};
  vector<int> count(P + 1), sats(P + 1);
  for (const auto& r : roots) {
    ASSERT_EQ(r.w.lo.q, r.w.hi.q);
    ASSERT_LT(r.w.lo.value(), r.w.hi.value());
    count[r.w.lo.q]++;
    sats[r.w.lo.q] += r.satellite;
  }
  for (int p = 2; p <= P; p++) {
    ASSERT_EQ(count[p], comps[p]) << tfm::format("period %d", p);
    // Satellites of period p: each component of period k | p, k < p, has phi(p/k) of them
    int want = 0;
    for (int k = 1; k < p; k++)
      if (p % k == 0) {
        int phi = 0;
        for (int a = 1; a <= p / k; a++) phi += gcd(a, p / k) == 1;
        want += comps[k] * phi;
      }
    ASSERT_EQ(sats[p], want) << tfm::format("period %d", p);
  }
  const auto has = [&](const Wake& w, const bool satellite) {
    for (const auto& r : roots)
      if (r.w == w) return r.satellite == satellite;
    return false;
  };
  // Known small roots
  ASSERT_TRUE(has(Wake{{1, 2}, {2, 2}}, true));    // -3/4
  ASSERT_TRUE(has(Wake{{3, 3}, {4, 3}}, false));   // -1.75, primitive period 3
  ASSERT_TRUE(has(Wake{{3, 4}, {4, 4}}, false));   // primitive period 4 in the 1/3 limb
  ASSERT_TRUE(has(Wake{{7, 4}, {8, 4}}, false));   // -1.94, primitive period 4
  ASSERT_TRUE(has(Wake{{6, 4}, {9, 4}}, true));    // -5/4, satellite of the period 2 disk
  // Every cardioid wake and every tuned wake of the period 2 disk is a satellite root
  for (int q = 2; q <= P; q++)
    for (int p = 1; p < q; p++)
      if (gcd(p, q) == 1) {
        ASSERT_TRUE(has(cardioid_wake(p, q), true)) << tfm::format("%d/%d", p, q);
        if (2*q <= P) ASSERT_TRUE(has(tune(cardioid_wake(1, 2), cardioid_wake(p, q)), true));
      }
  // Wakes are nested or disjoint
  for (size_t i = 0; i < roots.size(); i += 701)
    for (size_t j = 0; j < roots.size(); j += 503) {
      const auto &a = roots[i].w, &b = roots[j].w;
      const bool lo_in = b.contains(a.lo.value()), hi_in = b.contains(a.hi.value());
      ASSERT_EQ(lo_in, hi_in) << a << " vs " << b;
    }
}

TEST(lavaurs_within) {
  const int P = 16;
  const auto all = lavaurs(P);
  for (const auto& [p, q] : {std::pair(1, 3), std::pair(2, 5), std::pair(3, 7), std::pair(1, 4)}) {
    const auto w = cardioid_wake(p, q);
    vector<Root> inside;
    for (const auto& r : all)
      if (w.lo.value() <= r.w.lo.value() && r.w.hi.value() <= w.hi.value()) inside.push_back(r);
    auto sub = lavaurs(P, w);
    const auto key = [](const Root& r) { return std::tuple(r.w.lo.q, r.w.lo.k, r.w.hi.k); };
    std::sort(inside.begin(), inside.end(), [&](auto& a, auto& b) { return key(a) < key(b); });
    std::sort(sub.begin(), sub.end(), [&](auto& a, auto& b) { return key(a) < key(b); });
    ASSERT_EQ(sub.size(), inside.size()) << p << "/" << q;
    for (size_t i = 0; i < sub.size(); i++) {
      ASSERT_EQ(sub[i].w, inside[i].w);
      ASSERT_EQ(sub[i].satellite, inside[i].satellite);
    }
  }
}

TEST(lavaurs_within_long) {
  // Beyond 29 bits: roots inside an inner wake agree whether paired within it or within an outer wake
  const int P = 31;
  const auto outer = lavaurs(P, cardioid_wake(5, 11));
  const Root* inner = nullptr;
  for (const auto& r : outer)
    if (r.w.lo.q == 14 && !r.satellite) { inner = &r; break; }
  ASSERT_TRUE(inner != nullptr);
  const auto key = [](const Root& r) { return std::tuple(r.w.lo.q, r.w.lo.k, r.w.hi.k); };
  const auto in = [&](const Root& r) {  // r's wake inside inner's (exact, 128 bits)
    typedef unsigned __int128 U;
    const auto le = [](const Periodic& a, const Periodic& b) {
      return U(a.k) * ((uint64_t(1) << b.q) - 1) <= U(b.k) * ((uint64_t(1) << a.q) - 1);
    };
    return le(inner->w.lo, r.w.lo) && le(r.w.hi, inner->w.hi);
  };
  vector<Root> a, b = lavaurs(P, inner->w);
  for (const auto& r : outer) if (in(r)) a.push_back(r);
  std::sort(a.begin(), a.end(), [&](auto& x, auto& y) { return key(x) < key(y); });
  std::sort(b.begin(), b.end(), [&](auto& x, auto& y) { return key(x) < key(y); });
  ASSERT_LT(1u, b.size());
  ASSERT_EQ(a.size(), b.size());
  for (size_t i = 0; i < a.size(); i++) ASSERT_EQ(a[i].w, b[i].w);
  // maximal_tuning works past 29 bits: the inner root's own satellites are tuned by it
  const auto parent = maximal_tuning(b);
  int tuned = 0;
  for (size_t i = 0; i < b.size(); i++) tuned += parent[i] >= 0;
  ASSERT_LT(0, tuned);
}

TEST(maximal_tuning) {
  const int P = 16;
  const auto roots = lavaurs(P);
  const auto parent = maximal_tuning(roots);
  // Known cases: the period-4 disk at -1.3107 (6/15, 9/15) is tuned by the period-2 disk (1/3, 2/3); the 1/4
  // bulb (1/15, 2/15) and the period-3 bulbs are not
  int tuned4 = -1, bulb4 = -1, disk2 = -1;
  for (size_t i = 0; i < roots.size(); i++) {
    const auto& w = roots[i].w;
    if (w.lo == Periodic{6, 4} && w.hi == Periodic{9, 4}) tuned4 = int(i);
    if (w.lo == Periodic{1, 4} && w.hi == Periodic{2, 4}) bulb4 = int(i);
    if (w.lo == Periodic{1, 2} && w.hi == Periodic{2, 2}) disk2 = int(i);
    if (w.lo.q == 3) ASSERT_EQ(parent[i], -1);
  }
  ASSERT_TRUE(tuned4 >= 0 && bulb4 >= 0 && disk2 >= 0);
  ASSERT_EQ(parent[tuned4], disk2);
  ASSERT_EQ(parent[bulb4], -1);
  // Parents are non-renormalizable with smaller periods dividing p, and every root's copy decodes back
  vector<int64_t> N(P + 1), nr(P + 1);
  N[1] = 1;
  for (size_t i = 0; i < roots.size(); i++) {
    const int p = roots[i].w.lo.q;
    N[p]++;
    if (parent[i] < 0) { nr[p]++; continue; }
    const auto& w0 = roots[parent[i]].w;
    ASSERT_EQ(parent[parent[i]], -1);
    ASSERT_TRUE(p % w0.lo.q == 0 && w0.lo.q < p);
  }
  // Counting: tuned roots of period p are τ_W0(W') for non-renormalizable W0 of period d | p, 1 < d < p, and
  // any W' of period p / d, each in exactly one maximal copy
  for (int p = 2; p <= P; p++) {
    int64_t tuned = 0;
    for (int d = 2; d < p; d++)
      if (p % d == 0) tuned += nr[d] * N[p / d];
    ASSERT_EQ(N[p] - nr[p], tuned) << tfm::format("p %d: %d roots, %d non-renormalizable", p, N[p], nr[p]);
  }
}

TEST(wake_owner) {
  const auto roots = lavaurs(9);
  for (const int64_t n : {64, 1000, 4096, 12345}) {
    const auto owner = wake_owner(roots, n);
    ASSERT_EQ(int64_t(owner.size()), n/2 + 1);
    for (int64_t i = 0; i <= n/2; i++) {
      // Brute force: smallest wake containing i/n
      const double t = double(i) / n;
      int32_t best = -1;
      for (size_t r = 0; r < roots.size(); r++)
        if (roots[r].w.contains(t) && (best < 0 || roots[r].w.length() < roots[best].w.length()))
          best = int32_t(r);
      ASSERT_EQ(owner[i], best) << tfm::format("n %d, i %d", n, i);
    }
  }
}

}  // namespace
}  // namespace mandelbrot
