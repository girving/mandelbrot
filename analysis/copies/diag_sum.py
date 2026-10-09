"""Σ of the two-transit sector's constants D = lim a(k, label) k^4.

  python3 diag_sum.py labels area(.gz) shapes neg_area  (labels: "H<i> u|form|d|v" for the insertion-rule run; shapes and
                                                       neg_area: diag_negative.py's shapes file and areas)
Labels from the limb data use their own D; negative offsets d ≤ -(|u|+|v|)-1 not among them come from the stable-word
run for d ≥ -16, and beyond that from a c |d|^-6 model fitted on d in [-14, -8] (the 1/k extrapolation needs k ≫ |d|)."""
import sys, gzip
from collections import defaultdict
from decimal import Decimal as D, getcontext
getcontext().prec = 40

def neville(xs, ys):
    p = list(ys)
    for m in range(1, len(xs)):
        p = [(xs[i+m] * p[i] - xs[i] * p[i+1]) / (xs[i+m] - xs[i]) for i in range(len(p) - 1)]
    return p[0]

def constants(path, kmin=16):
    a = defaultdict(dict)
    for l in (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)):
        r = l.split()
        if r[3] == 'failed': continue
        nm, k = r[0].rsplit('_', 1); a[nm][int(k)] = D(float(r[5])) + D(float(r[10]))
    out = {}
    for nm, d in a.items():
        ks = [k for k in sorted(d) if k >= kmin]
        if len(ks) < 5: continue
        e1 = neville([D(1) / k for k in ks], [d[k] * k**4 for k in ks])
        e2 = neville([D(1) / k for k in ks[1:]], [d[k] * k**4 for k in ks[1:]])
        out[nm] = (float(e1), float(abs(e1 - e2)))
    return out

if __name__ == '__main__':
    lab = dict(l.split() for l in open(sys.argv[1]))
    Dl = {lab[nm]: v for nm, v in constants(sys.argv[2]).items()}
    shapes = [l.strip() for l in open(sys.argv[3]) if l.strip()]
    neg = constants(sys.argv[4])
    S_lab = sum(v for v, e in Dl.values()); E_lab = sum(e for v, e in Dl.values())
    S_neg = 0.0; tail = 0.0
    for i, sh in enumerate(shapes):
        u, form, v = sh.split('|')
        vals = {}
        for nm, (val, e) in neg.items():
            si, d = nm[1:].split('d')
            if int(si) == i: vals[int(d)] = val
        for d, val in vals.items():
            if d >= -16 and '%s|%s|%d|%s' % (u, form, d, v) not in Dl: S_neg += val
        # model beyond -16, per parity, from d in [-14, -8]
        for par in (0, 1):
            fit = [vals[d] * abs(d)**6 for d in vals if -14 <= d <= -8 and d % 2 == par]
            if fit:
                c = sum(fit) / len(fit)
                tail += sum(c * abs(d)**-6.0 for d in range(-17, -100000, -1) if d % 2 == par)
    print('labels from the limbs: %d, Σ D = %.10e (extrapolation spread Σ %.1e)' % (len(Dl), S_lab, E_lab))
    print('negative offsets (stable words, d ≥ -16, new labels): %.6e; model tail d < -16: %.2e' % (S_neg, tail))
    print('two-transit sector Σ D ≈ %.10e' % (S_lab + S_neg + tail))
