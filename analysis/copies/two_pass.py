"""Two-pass seahorse families: keys with a second run of the gate word, the dominant families of the census.

  python3 two_pass.py M > lines;  limb_families custom 2 < lines > jobs   (names <seq>m<m>, k from K to 4K)
  bulb_batch --out out < jobs;  python3 two_pass.py fit out

E0: 00(10)^m 1 / 00(10)^{m-1}110 (j = 2m)          E2: 0011(01)^{m-1}001 / 0011(01)^{m-1}010 (j = 2m+2)
O1: 0011(01)^m / 0011(01)^{m-1}10 (j = 2m+1)       OL: 010(01)^{m-1}001 / 010(01)^{m-1}010 (j = 2m+1)
For each m the limb index runs over k = K·{1, 5/4, 3/2, 2, 5/2, 3, 4}, K = max(12, j), inside the stable range
j ≤ 2k + 1."""
import sys, math
from collections import defaultdict
SEQ = {
    'E0': lambda m: ('00' + '10' * m + '1', '00' + '10' * (m - 1) + '110'),
    'E2': lambda m: ('0011' + '01' * (m - 1) + '001', '0011' + '01' * (m - 1) + '010'),
    'O1': lambda m: ('0011' + '01' * m, '0011' + '01' * (m - 1) + '10'),
    'OL': lambda m: ('010' + '01' * (m - 1) + '001', '010' + '01' * (m - 1) + '010'),
}
def neville(xs, ys):
    p = list(ys)
    for m in range(1, len(xs)):
        p = [(xs[i+m] * p[i] - xs[i] * p[i+1]) / (xs[i+m] - xs[i]) for i in range(len(p) - 1)]
    return p[0]

def fit(path):
    """C(m) = lim_k a k^4 (error: against the fit without the smallest k), local exponent in m, σ(m) and its steps"""
    fam = defaultdict(dict)
    for l in open(path):
        r = l.split(); name, k = r[0].rsplit('_', 1)
        fam[name][int(k)] = (complex(float(r[3]), float(r[4])) + 0.75, float(r[5]) + float(r[10]))
    res = defaultdict(dict)
    for name, F in fam.items():
        seq, m = name.split('m'); ks = sorted(F)
        C = neville([1 / k for k in ks], [F[k][1] * k**4 for k in ks])
        C2 = neville([1 / k for k in ks[1:]], [F[k][1] * k**4 for k in ks[1:]])
        s = neville([1 / k for k in ks], [1j * math.pi / F[k][0] - 2 * k for k in ks])
        res[seq][int(m)] = (C, abs(C2 / C - 1), s)
    for seq in sorted(res):
        R = res[seq]; print(seq)
        for m in sorted(R):
            C, e, s = R[m]
            ex = math.log(R[m-1][0] / C) / math.log(m / (m - 1)) if m - 1 in R else 0
            ds = s - R[m-1][2] if m - 1 in R else 0
            print('  m=%2d C=%.6e (err %.0e) local exponent %6.3f  σ=%+.5f%+.5fi  Δσ·m=%+.4f%+.4fi' %
                  (m, C, e, ex, s.real, s.imag, (ds * m).real, (ds * m).imag))

if __name__ == '__main__' and sys.argv[1] == 'fit':
    fit(sys.argv[2])
elif __name__ == '__main__':
    for name, f in SEQ.items():
        for m in range(1, int(sys.argv[1]) + 1):
            lo, hi = f(m); j = len(lo) - 3; K = max(12, j)
            ks = sorted({round(K * x) for x in (1, 1.25, 1.5, 2, 2.5, 3, 4)})
            print('%sm%d %s %s %s' % (name, m, lo, hi, ','.join(map(str, ks))))
