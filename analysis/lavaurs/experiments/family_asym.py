"""Finite-k families in M of model components at several levels (probe of the analytic asymptotics' uniformity).
For a model component (level r, σ in the q = 2 strip coordinate, model constant C) the M family is c_k = -3/4 + iπ/(2k
+ σ_M(k)), σ_M → -2 conj(σ + 7/2 + 3πi/8) (mod 2), period r (2k + 1) + const.  Members are found top-down in k: at
large k from the model prediction (period window, the const fixed by the first hit), then by Neville extrapolation of
σ_M(k) in 1/k.  Prints jobs for bulb_batch ("name 0 c_re c_im 0 P").

  python3 family_asym.py components.txt > jobs   (lines: name r σ_re σ_im C)"""
import cmath, math, sys
import numpy as np
KS = [2048, 1448, 1024, 724, 512, 362, 256, 181, 128, 108, 90, 76, 64, 54, 45, 38, 32, 27, 22, 19, 16, 13, 11, 9, 8, 7, 6, 5, 4]
OFF = -3.5 - 3j * math.pi / 8
def newton(c, P, iters=60):   # vectorized over candidate periods P (array) at seeds c (array)
    c = np.array(c, dtype=complex); P = np.array(P); ok = np.zeros(len(c), bool)
    for _ in range(iters):
        z = np.zeros_like(c); dz = np.zeros_like(c); zP = np.zeros_like(c); dP = np.zeros_like(c)
        with np.errstate(all='ignore'):
            for i in range(1, P.max() + 1):
                dz = 2 * z * dz + 1; z = z * z + c
                m = P == i
                if m.any(): zP[m] = z[m]; dP[m] = dz[m]
            st = zP / dP
        st[~np.isfinite(st)] = 0; c = c - st
        ok = np.abs(st) < 1e-15 * np.abs(c)
        if ok.all(): break
    return c, ok
def exact_period(c, P):
    z = 0j
    for k in range(1, P):
        z = z * z + c
        if abs(z) < 1e-9 and P % k == 0: return False
    return True
def sM_of(c, k): return 1j * math.pi / (c + 0.75) - 2 * k
def neville(xs, ys, x):
    P = list(ys)
    for m in range(1, len(xs)):
        for i in range(len(xs) - m): P[i] = ((x - xs[i + m]) * P[i] + (xs[i] - x) * P[i + 1]) / (xs[i] - xs[i + m])
    return P[0]
for l in open(sys.argv[1]):
    nm, r, sr, si, C = l.split(); r = int(r); C = float(C); s = complex(float(sr), float(si))
    s0 = -2 * (s - OFF).conjugate(); shift = math.floor(s0.real / 2) * 2; s0 -= shift
    rad = 0.72 * math.sqrt(C)          # the component's radius in σ_M units
    found = []; const = None
    for k in KS:
        if len(found) >= 3:
            pred = neville([1 / kk for kk, _, _ in found[-6:]], [sm for _, sm, _ in found[-6:]], 1 / k)
        else:
            pred = s0 + 2 * 0.35 / k   # e(k) ≈ -0.35/k in σ_gen
        c0 = -0.75 + 1j * math.pi / (2 * k + pred)
        Ps = [r * (2 * k + 1) + d for d in (range(-12 * r - 20, 12 * r + 21) if const is None else [const])]
        cs, ok = newton([c0] * len(Ps), Ps)
        best = None
        for c, o, P in zip(cs, ok, Ps):
            if not o: continue
            d = abs(sM_of(c, k) - pred)
            if d < 3 * rad * (1 + (300 / k if len(found) < 3 else 0)) and exact_period(c, P) and (best is None or d < best[0]): best = (d, c, P)
        if best is None:
            print('# %s k %d: no member (pred σ_M %.6f%+.6fi)' % (nm, k, pred.real, pred.imag), file=sys.stderr)
            if found: break
            continue
        d, c, P = best; const = P - r * (2 * k + 1)
        found.append((k, sM_of(c, k), P))
        print('%s_%d 0 %.17g %.17g 0 %d' % (nm, k, c.real, c.imag, P), flush=True)
        print('# %s k %d P %d miss %.2e radii' % (nm, k, P, d / rad), file=sys.stderr)
