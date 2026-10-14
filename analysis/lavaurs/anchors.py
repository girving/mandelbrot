"""Tunings of a satellite source (the limb's bulb, a Farey bulb) by M's primitive components, located through anchors.

A satellite's copy of M is strongly distorted in σ, so guesses from its multiplier map fail.  But its satellite tree
can be placed exactly: the p/q satellite of a component is rooted at its multiplier map's point e^{2πi p/q}
(lavaurs_area --satellites), and on the M side the same recursion from the cardioid gives the matching centers c.
These (c, σ) pairs, plus known tunings, anchor the copy map χ: c ↦ σ; for a primitive X a local affine fit
χ(c) ≈ a + b c + d c̄ over the nearest anchors predicts U*X, and Newton (r = p r_U, n = p n_U + p - 1) finds it.

  python3 anchors.py --source "bulb 1 1 -1.0074583370365449 0.16135210336429348" --prims prims.txt [--qprod 12]
prints "U*X r n re im" for each found tuning (verify with lavaurs_area: cusp ≈ 0)."""
import argparse, cmath, math, os, subprocess, sys
import numpy as np

BIN = os.path.join(os.path.dirname(__file__), '../../build/release/lavaurs_area')

# ---- the M side: satellites of components of M
def fpz(c, p, z=0j):
    dz = 1
    for _ in range(p): dz = 2 * z * dz; z = z * z + c
    return z, dz

def mult_point(c0, P, mu, K=40):
    """The parameter on the component (center c0, period P) whose attracting cycle has multiplier mu"""
    c, z = c0, 0j
    for k in range(1, K + 1):
        m = mu * k / K
        for _ in range(60):
            w, dwz, dwc, d2 = z, 1, 0, 0
            for _ in range(P):
                d2 = 2 * (dwz * dwz + w * d2); dwc = 2 * w * dwc + 1; dwz = 2 * w * dwz; w = w * w + c
            h = 1e-8
            w2, dz2 = z, 1
            for _ in range(P): dz2 = 2 * w2 * dz2; w2 = w2 * w2 + (c + h)
            dmu = (dz2 - dwz) / h
            J = np.array([[dwz - 1, dwc], [d2, dmu]], dtype=complex)
            dz_, dc_ = np.linalg.solve(J, np.array([w - z, dwz - m], dtype=complex))
            z -= dz_; c -= dc_
            if abs(dz_) + abs(dc_) < 1e-15: break
    return c

def center(c, p):
    for _ in range(100):
        z, dz = 0j, 0j
        for _ in range(p): dz = 2 * z * dz + 1; z = z * z + c
        step = z / dz; c -= step
        if abs(step) < 1e-15: return c
    return None

def m_satellite(c0, P, p, q):
    th = 2 * math.pi * p / q
    root = mult_point(c0, P, cmath.exp(1j * th)); inner = mult_point(c0, P, 0.95 * cmath.exp(1j * th))
    for f in (1, 2, 4, 8, 16):
        c = center(root + f * (root - inner), P * q)
        if c is not None and abs(c - c0) > 1e-9:
            # exact period P q
            z = 0j
            per = None
            for i in range(1, P * q + 1):
                z = z * z + c
                if abs(z) < 1e-9: per = i; break
            if per == P * q: return c
    return None

# ---- the σ side
def sigma_satellites(name, r, n, s, qmax):
    out = subprocess.run([BIN, '--satellites'], input='%s %d %d %.17g %.17g %d\n' % (name, r, n, s.real, s.imag, qmax),
                         capture_output=True, text=True, timeout=3600).stdout
    res = {}
    for l in out.splitlines():
        f = l.split()
        if 'failed' in l: continue
        pq = f[0].split('/')[-1]
        res[pq] = (int(f[1]), int(f[2]), complex(float(f[3]), float(f[4])))
    return res

def build(src, qprod, qmax):
    """Anchor pairs (c, σ, path) for the satellite tree of src = (name, r, n, σ) down to denominator product qprod"""
    name, r, n, s = src
    anchors = [(0j, s, '')]
    frontier = [('', 0j, 1, r, n, s, 1)]   # (path, c, period in M, r, n, σ, product of q)
    while frontier:
        path, c, P, rr, nn, ss, qp = frontier.pop()
        qm = min(qmax, qprod // qp)
        if qm < 2: continue
        sats = sigma_satellites(name + path, rr, nn, ss, qm)
        for pq, (r2, n2, s2) in sats.items():
            p_, q_ = map(int, pq.split(':'))
            c2 = m_satellite(c, P, p_, q_)
            if c2 is None: continue
            anchors.append((c2, s2, path + '/' + pq))
            frontier.append((path + '/' + pq, c2, P * q_, r2, n2, s2, qp * q_))
    return anchors

def predict(anchors, cX, k=8):
    d = sorted(anchors, key=lambda a: abs(a[0] - cX))[:k]
    A = np.array([[1, a[0], a[0].conjugate()] for a in d], dtype=complex)
    w = np.array([1 / (abs(a[0] - cX) + 1e-3) for a in d])
    coef, *_ = np.linalg.lstsq(A * w[:, None], np.array([a[1] for a in d]) * w, rcond=None)
    return coef[0] + coef[1] * cX + coef[2] * cX.conjugate()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--source', required=True)
    ap.add_argument('--prims', required=True)
    ap.add_argument('--known', default='')        # extra anchors "c_re c_im s_re s_im"
    ap.add_argument('--qprod', type=int, default=12)
    ap.add_argument('--qmax', type=int, default=6)
    ap.add_argument('--anchors-out', default='')
    a = ap.parse_args()
    f = a.source.split(); src = (f[0], int(f[1]), int(f[2]), complex(float(f[3]), float(f[4])))
    anchors = build(src, a.qprod, a.qmax)
    if a.known:
        for l in open(a.known):
            g = l.split(); anchors.append((complex(float(g[0]), float(g[1])), complex(float(g[2]), float(g[3])), 'known'))
    print('# %d anchors' % len(anchors), file=sys.stderr)
    if a.anchors_out:
        with open(a.anchors_out, 'w') as fo:
            for c, s, path in anchors: fo.write('%.15g %.15g %.17g %.17g %s\n' % (c.real, c.imag, s.real, s.imag, path or '.'))
    name, r, n, s = src
    for l in open(a.prims):
        g = l.split(); p = int(g[0]); cX = complex(float(g[1]), float(g[2]))
        guess = predict(anchors, cX)
        print('%s*P%d(%.4f%+.4fi) %d %d %.17g %.17g' % (name, p, cX.real, cX.imag, r * p, p * n + p - 1, guess.real, guess.imag))

if __name__ == '__main__':
    main()
