"""Model seeds from M-side runs: σ_M = lim (iπ/δ_k - 2k) by Neville in 1/k over a family's centers (δ = c + 3/4), and the
model phase σ = -σ_M/2 + 3πi/8 + 1/2.

  python3 seeds.py area.txt[.gz] [kmin] > seeds   (lines "name σ_re σ_im")   (bulb_batch / local output: name_k ... c)"""
import sys, gzip, math
from collections import defaultdict

def neville(xs, ys):
    p = list(ys)
    for m in range(1, len(xs)):
        p = [(xs[i+m] * p[i] - xs[i] * p[i+1]) / (xs[i+m] - xs[i]) for i in range(len(p) - 1)]
    return p[0]

def seeds(path, kmin=16, npts=10):
    c = defaultdict(dict)
    for l in (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)):
        r = l.split()
        if r[3] == 'failed': continue
        nm, k = r[0].rsplit('_', 1)
        c[nm][int(k)] = complex(float(r[3]), float(r[4])) + 0.75
    out = {}
    for nm, d in c.items():
        ks = [k for k in sorted(d) if k >= kmin]
        if len(ks) < 4: continue
        ks = ks[-npts:]
        sM = neville([1 / k for k in ks], [1j * math.pi / d[k] - 2 * k for k in ks])
        out[nm] = -sM / 2 + 3j * math.pi / 8 + 0.5
    return out

if __name__ == '__main__':
    for nm, s in seeds(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 16).items():
        print('%s %.15f %.15f' % (nm, s.real, s.imag))
