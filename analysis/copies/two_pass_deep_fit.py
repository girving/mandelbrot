"""Deep two-pass size estimates (limb_families size; cluster_deep.sh): the k- and m-dependence of a_s = A_card |s|².

  python3 two_pass_deep_fit.py deep.txt[.gz]

Per sequence and m: g(t) = k^4 m^6 a_s (sin πt / πt)^6 with t = m/k (the gate model's shape divided out), its spread
over t ≤ 1/2 (flat if the model is exact), and its t → 0 limit C_s(m) m^6 from a polynomial fit in t over t ≤ 1/4.
Then C_s(m) m^6 against m, with local exponents: a constant means C ∝ m^-6, a drift shows corrections."""
import sys, gzip, math
from collections import defaultdict
import numpy as np
A_CARD = 3 * math.pi / 8

def load(path):
    a = defaultdict(dict)
    for l in (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)):
        r = l.split(); name, k = r[0].rsplit('_', 1); seq, m = name.split('m')
        a[(seq, int(m))][int(k)] = A_CARD * float(r[6])
    return a

if __name__ == '__main__':
    A = load(sys.argv[1])
    for seq in sorted({s for s, m in A}):
        print(seq)
        prev = None
        for m in sorted(m for s, m in A if s == seq):
            d = A[(seq, m)]
            g = {k: d[k] * k**4 * m**6 * (math.sin(math.pi * m / k) / (math.pi * m / k))**6 for k in d}
            ks = [k for k in d if m / k <= 0.25]
            if len(ks) < 6: continue
            t = np.array([m / k for k in ks]); y = np.array([g[k] for k in ks])
            c = np.polynomial.polynomial.polyfit(t, y, 5)
            c4 = np.polynomial.polynomial.polyfit(t, y, 4)
            half = [g[k] for k in d if abs(m / k - 0.5) < 0.02]
            ex = '' if prev is None else '%6.3f' % (math.log(prev[1] / c[0]) / math.log(m / prev[0]))
            print('  m=%5d  C m^6 = %.8e (fit error %.0e)  local exponent of C m^6 %s   g(1/2)/g(0) = %s   slope dg/dt / g(0) = %+.3f' %
                  (m, c[0], abs(c[0] - c4[0]) / c[0], ex or '   -- ', '%.4f' % (half[0] / c[0]) if half else '  --  ', c[1] / c[0]))
            prev = (m, c[0])
