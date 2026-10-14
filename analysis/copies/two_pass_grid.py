"""Two-pass families on the grid a(k, m), m < k ≤ KMAX (cluster_families.sh): scaling in t = m/k.

  python3 two_pass_grid.py grid.txt[.gz]

Per sequence: g(k, m) = a k^4 m^6 against t = m/k, tabulated for several m at fixed t (a collapse onto Φ(t) means
the passes couple only through m/k), the same divided by the model's (πt / sin πt)^6, and the sums
S(m) = Σ_{k>m} a(k, m) with their local exponent in m.  Model: near α = -1/2 with δ = iπ/(2k), f² is a rotation by
θ ≈ π/k with a cubic term, linear in u = 1/w² (w = z - α, normalized): u ↦ e^{2iθ}(u + 2); a partial pass of m turns
ending at the exit starts at |u_0| = 2|sin mθ / sin θ| ≈ (2k/π) sin πt, contributing |u_0|^{3/2} to Λ, so area ∝
|u_0|^-6 ∝ m^-6 (πt / sin πt)^6."""
import sys, gzip, math
from collections import defaultdict

def load(path):
    a = defaultdict(dict)
    for l in (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)):
        r = l.split()
        if r[3] == 'failed': continue
        name, k = r[0].rsplit('_', 1); seq, m = name.split('m')
        a[seq][(int(k), int(m))] = float(r[5]) + float(r[10])
    return a

if __name__ == '__main__':
    A = load(sys.argv[1])
    for seq in sorted(A):
        a = A[seq]; ms = sorted({m for k, m in a}); kmax = max(k for k, m in a)
        print('%s: %d points, m ≤ %d, k ≤ %d' % (seq, len(a), ms[-1], kmax))
        print('   g = a k^4 m^6 at t = m/k:   ' + '  '.join('t=%-9.3g' % t for t in (1 / 2, 1 / 4, 1 / 8, 1 / 16)))
        for m in (2, 4, 8, 16, 32):
            row = []
            for d in (2, 4, 8, 16):
                k = m * d
                row.append('%.5e' % (a[(k, m)] * k**4 * m**6) if (k, m) in a else '     --    ')
            print('   m=%2d                        %s' % (m, '  '.join(row)))
        print('   g (sin πt / πt)^6:           ' + '  '.join('t=%-9.3g' % t for t in (1 / 2, 1 / 4, 1 / 8, 1 / 16)))
        for m in (4, 8, 16, 32):
            row = []
            for d in (2, 4, 8, 16):
                k = m * d; t = m / k
                row.append('%.5e' % (a[(k, m)] * k**4 * m**6 * (math.sin(math.pi * t) / (math.pi * t))**6) if (k, m) in a else '     --    ')
            print('   m=%2d                        %s' % (m, '  '.join(row)))
        S = {m: sum(v for (k, mm), v in a.items() if mm == m) for m in ms}
        print('   S(m) = Σ_{m<k≤%d} a:' % kmax, '  '.join('m=%d %.4e (%.2f)' % (m, S[m], math.log(S[m // 2] / S[m]) / math.log(2) if m // 2 in S and m > 1 else 0)
                                                   for m in (1, 2, 4, 8, 16, 32, 64) if m in S))
