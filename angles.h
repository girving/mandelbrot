// External angles: exact periodic binary angles, cardioid wakes, and Douady-Hubbard tuning
#pragma once

#include <cstdint>
#include <ostream>
#include <vector>
namespace mandelbrot {

using std::ostream;

// The angle k/(2^q - 1), whose binary expansion repeats the q-bit word k.  q need not be the minimal period.
struct Periodic {
  uint64_t k;
  int q;

  double value() const { return double(k) / double((uint64_t(1) << q) - 1); }
  bool operator==(const Periodic& a) const { return k == a.k && q == a.q; }
  friend ostream& operator<<(ostream& out, const Periodic& a) { return out << a.k << "/(2^" << a.q << "-1)"; }
};

// A wake: the arc of angles between the two external rays landing at a root, with lo < hi as numbers in [0,1)
struct Wake {
  Periodic lo, hi;

  double length() const { return hi.value() - lo.value(); }
  bool contains(const double t) const { return lo.value() < t && t < hi.value(); }
  bool operator==(const Wake& w) const { return lo == w.lo && hi == w.hi; }
  friend ostream& operator<<(ostream& out, const Wake& w) { return out << '(' << w.lo << ", " << w.hi << ')'; }
};

// Wake of the p/q satellite bulb of the main cardioid (gcd(p,q) = 1, 1 <= p < q <= 62).  Its length is 1/(2^q-1).
Wake cardioid_wake(const int p, const int q);

// Tuning by the component whose root has wake w: replace each binary digit of a with w.lo's word (for 0)
// or w.hi's word (for 1).  Requires w.lo.q == w.hi.q and w.lo.q * a.q <= 62.
Periodic tune(const Wake& w, const Periodic& a);
Wake tune(const Wake& w, const Wake& v);

// Inverse tuning on the leading digits of the binary fraction bits / 2^nbits: decode successive words of
// length w.lo.q (w.lo's word -> 0, w.hi's word -> 1), stopping at the first other word or after D digits.
// Returns the number of digits decoded, with the decoded digits in prefix.
int untune(const Wake& w, const uint64_t bits, const int nbits, const int D, uint64_t& prefix);

// Roots of all hyperbolic components of periods 2 through max_period (<= 24), as the wakes bounded by
// their two landing rays, via Lavaurs' algorithm: for each period in turn, pair each unpaired angle of
// exact period p (in increasing order) with the next unpaired one whose chord crosses no earlier chord.
struct Root {
  Wake w;          // w.lo.q = w.hi.q = period
  bool satellite;  // Satellite roots have both rays on one doubling cycle; primitive roots have two cycles
};
std::vector<Root> lavaurs(const int max_period);

}  // namespace mandelbrot
