"""Accuracy of the area-from-center formula C ≈ C_t |Θ_r' Π H'|^-2 against exact areas (boundary tracing) for the
heavy children of a gl_census level: the census's raw output (formula) vs its level file (exact for C > exact-min).

  python3 gl_formula_err.py single.txt prefix level"""
import sys, math
single, prefix, r = sys.argv[1], sys.argv[2], int(sys.argv[3])
Ct = {}; nf = ns = 0
for l in open(single):
    f = l.split()
    if f[0].startswith('sat'): Ct['B%d' % ns] = float(f[3]); ns += 1
    else: Ct['F%d' % nf] = float(f[3]); nf += 1
key = lambda c: (round(c.real % 1, 7) % 1, round(c.imag, 7))
exact = {}
for l in open('%s_r%d.txt' % (prefix, r)):
    f = l.split(); exact[key(complex(float(f[3]), float(f[4])))] = (float(f[5]), f[6] == '1')
form = {}
for l in open('%s_raw_r%d.txt' % (prefix, r)):
    f = l.split(); nm = f[0].rsplit('|', 2)[0]; t = nm.split('~')[1]
    if int(f[4]) or float(f[3]) >= 1: continue
    k = key(complex(float(f[1]), float(f[2])))
    if k not in form: form[k] = (Ct[t] * float(f[3]), Ct[t])
rows = []
for k, (Cf, Ctg) in form.items():
    if k in exact and not exact[k][1] and exact[k][0] > 1e-8 and abs(exact[k][0] - Cf) > 0:
        rows.append((exact[k][0], Cf, Ctg))
rows.sort(reverse=True)
errs = sorted(abs(math.log(Cf / Ce)) for Ce, Cf, _ in rows)
q = lambda f: errs[min(len(errs) - 1, int(f * len(errs)))]
print('level %d: %d exact primitive children > 1e-8; |log(formula/exact)| median %.2e, IQR-ish 25/75%% %.2e %.2e, 90%% %.2e' % (
    r, len(rows), q(.5), q(.25), q(.75), q(.9)))
# against the relative size of the child (C/C_t): the formula's error should scale with the target-relative size
for lo, hi in ((1e-12, 1e-6), (1e-6, 1e-4), (1e-4, 1e-2), (1e-2, 1)):
    e = sorted(abs(math.log(Cf / Ce)) for Ce, Cf, Ctg in rows if lo <= Ce / Ctg < hi)
    if e: print('   C/C_t in [%.0e, %.0e): %4d children, median |log| %.2e' % (lo, hi, len(e), e[len(e) // 2]))
