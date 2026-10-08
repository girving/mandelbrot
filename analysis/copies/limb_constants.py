"""Family constants C_w = lim a_w(k) k^4 for all seahorse-valley families (limb_families.cc), and their decay in the
extra period j.

  limb_families jobs J k1,k2,... > jobs; bulb_batch --out out < jobs; python3 limb_constants.py out

Per family: Neville in 1/k through all its k (error estimate: against the fit without the smallest k), the phase
σ_w = lim iπ/δ − 2k.  Per j: count, Σ C_w, largest C_w, the share of the top 10%, and the phase-plane area
Σ 16 C_w/π² (|dδ/dσ|² = π²/|2k+σ|^4 ≈ π²/(16 k^4))."""
import sys, math
from collections import defaultdict

def neville(xs, ys):
    p = list(ys)
    for m in range(1, len(xs)):
        p = [(xs[i+m] * p[i] - xs[i] * p[i+1]) / (xs[i+m] - xs[i]) for i in range(len(p) - 1)]
    return p[0]

fam = defaultdict(dict)
for l in open(sys.argv[1]):
    r = l.split()
    if r[3] == 'failed': continue
    j, i, k = r[0].split('_')
    fam[(int(j[1:]), int(i))][int(k)] = (complex(float(r[3]), float(r[4])) + 0.75, float(r[5]) + float(r[10]))
rows = []
for (j, i), F in fam.items():
    ks = sorted(F)
    C = neville([1 / k for k in ks], [F[k][1] * k**4 for k in ks])
    C2 = neville([1 / k for k in ks[1:]], [F[k][1] * k**4 for k in ks[1:]])
    sig = neville([1 / k for k in ks], [1j * math.pi / F[k][0] - 2 * k for k in ks])
    rows.append((j, i, C, abs(C - C2) / C, sig))
if len(sys.argv) > 2:
    with open(sys.argv[2], 'w') as f:
        for j, i, C, e, s in sorted(rows): f.write('%d %d %.12e %.1e %.10f %.10f\n' % (j, i, C, e, s.real, s.imag))
print('relative error estimate of C_w: median %.1e, max %.1e' % (sorted(r[3] for r in rows)[len(rows) // 2], max(r[3] for r in rows)))
print('  j     N   Σ C_w        ratio   max C_w      top10%%  phase area   Σ C_w/N')
prev = None
tot = 0
for j in sorted({r[0] for r in rows}):
    Cs = sorted((r[2] for r in rows if r[0] == j), reverse=True)
    S = sum(Cs); tot += S
    top = sum(Cs[:max(1, len(Cs) // 10)]) / S
    print('%3d %5d  %.4e  %s  %.4e  %5.1f%%  %.4e  %.3e' % (j, len(Cs), S, '%.4f' % (S / prev) if prev else '  --  ',
          Cs[0], 100 * top, 16 * S / math.pi**2, S / len(Cs)))
    prev = S
print('Σ over j: %.6e  (phase-plane area %.6e)' % (tot, 16 * tot / math.pi**2))
