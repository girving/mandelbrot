// External angle tests

#include "angles.h"
#include "tests.h"
#include <algorithm>
#include <numeric>
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

TEST(lavaurs) {
  const int P = 14;
  const auto roots = lavaurs(P);
  // Component counts per period: OEIS A000740
  const int comps[] = {0, 1, 1, 3, 6, 15, 27, 63, 120, 252, 495, 1023, 2010, 4095, 8127};
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
  for (size_t i = 0; i < roots.size(); i += 7)
    for (size_t j = 0; j < roots.size(); j += 5) {
      const auto &a = roots[i].w, &b = roots[j].w;
      const bool lo_in = b.contains(a.lo.value()), hi_in = b.contains(a.hi.value());
      ASSERT_EQ(lo_in, hi_in) << a << " vs " << b;
    }
}

}  // namespace
}  // namespace mandelbrot
