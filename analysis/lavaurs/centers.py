"""Single-transit centers of the Lavaurs model: σ with F^n(g_σ(v)) = 1/2 (the critical point), v = -1/4 the
critical value, by vectorized Newton from a grid.  (Prototype, double precision.)"""
import numpy as np
from lavaurs import A_COEF, C_LOG, phi_a_entry

def series_np(w, upper):
    L = np.where(upper, np.log(np.where(w.imag > 0, w, -w)), 0)
    s = 0.25 / (w * w) + 0.25 / w + C_LOG * L
    d = -0.5 / w**3 - 0.25 / (w * w) + C_LOG / w
    p = np.ones_like(w)
    for j, a in enumerate(A_COEF, 1):
        d = d + j * a * p; p = p * w; s = s + a * p
    return s, d

def psi_plus_np(zeta, R=100.0):
    m = np.maximum(0, np.ceil(zeta.real + np.maximum(R, 2 * np.abs(zeta.imag)))).astype(int)
    zl = zeta - m
    w = 1j / (2 * np.sqrt(-zl + 0j))
    for _ in range(40):
        s, d = series_np(w, True)
        w = w - (s - zl) / d
    s, d = series_np(w, True)
    dw = 1 / d
    for i in range(2 * int(m.max())):
        act = i < 2 * m
        dw = np.where(act, dw * (-1 + 2 * w), dw); w = np.where(act, -w + w * w, w)
    return w, dw

def psi_np(zeta, petal):
    if petal > 0: return psi_plus_np(zeta)
    w, dw = psi_plus_np(zeta - 0.5)
    return -w + w * w, (-1 + 2 * w) * dw

def endpoint(sigma, n, petal, zeta0):
    w, dw = psi_np(zeta0 + sigma, petal)
    for _ in range(n): dw = (-1 + 2 * w) * dw; w = -w + w * w
    return w, dw

def find(n, petal, zeta0, grid, iters=40):
    s = grid.copy()
    for _ in range(iters):
        with np.errstate(all='ignore'):
            w, dw = endpoint(s, n, petal, zeta0)
            step = (w - 0.5) / dw
        s = np.where(np.isfinite(step) & (np.abs(step) < 1), s - step, np.nan)
    with np.errstate(all='ignore'):
        w, dw = endpoint(s, n, petal, zeta0)
    ok = np.isfinite(w) & (np.abs(w - 0.5) < 1e-11)
    return s[ok], dw[ok]

if __name__ == '__main__':
    zeta0, dz0, petal_v = phi_a_entry(-0.25 + 0j)
    print('Φ_a(v) = %s, entering petal %d' % (zeta0, petal_v))
    re, im = np.meshgrid(np.linspace(-1, 1, 41), np.linspace(-2.5, 2.5, 101))
    grid = (re + 1j * im).ravel()
    for cross in (False, True):
        petal = -petal_v if cross else petal_v
        for n in range(0, 4):
            s, d = find(n, petal, zeta0, grid)
            u = np.unique(np.round(s, 9))
            sizes = sorted(((1 / abs(dd), ss) for ss, dd in zip(s, d)), reverse=True)
            seen = []; best = []
            for sz, ss in sizes:
                if all(abs(ss - t) > 1e-6 for t in seen): seen.append(ss); best.append((sz, ss))
            print('cross=%s n=%d: %d centers; largest 1/|∂endpoint/∂σ|: %s' % (cross, n, len(u),
                  '  '.join('%.3g@%.5f%+.5fi' % (sz, ss.real, ss.imag) for sz, ss in best[:5])))
