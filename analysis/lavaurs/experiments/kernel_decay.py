"""Depth profile of the universal decoration kernel of the quadratic family (layer 3).  python3 kernel_decay.py"""
# m_j(W) = Σ_{P_j(c')=W} |(g^j)'(c')|^-2 |P_j'(c')|^-2 (quadratic family): roots by parameter-ray tracing, then Newton
import numpy as np, math, cmath
def ev(c, j):
    z = c.copy(); dz = np.ones_like(c); pr = np.ones_like(c)
    for _ in range(j): pr *= 2 * z; dz = 2 * z * dz + 1; z = z * z + c
    return z, dz, pr   # z = P_j(c) (P_0 = c), dz = P_j', pr = Π_{i<j} 2 P_i
def roots(W, j):
    N = 2 ** j; k = np.arange(N); a0 = cmath.phase(W) / (2 * math.pi)
    c = None; t = math.log(100.0); tend = math.log(abs(W)) / N * 1.0
    while True:
        n = max(0, math.ceil(math.log2(12 / t))); n = min(n, j)
        ph = ((a0 + k) * 2.0 ** n / N) % 1.0     # fine for j ≤ 30
        tgt = np.exp(2.0 ** n * t + 2j * math.pi * ph)
        if c is None: c = np.exp(t + 2j * math.pi * (a0 + k) / N)
        for _ in range(6):
            z, dz, _ = ev(c, n); c = c - (z - tgt) / dz
        if t <= tend: break
        t = max(tend, t * 2 ** -0.125)
    for _ in range(30):
        z, dz, _ = ev(c, j); c = c - (z - W) / dz
    return c
for W in [3.0, 5.0, 10.0, 2 + 3j]:
    out = []
    for j in range(0, 17):
        c = roots(W, j) if j else np.array([complex(W)])
        z, dz, pr = ev(c, j)
        nd = len(np.unique(np.round(c, 8))); res = np.max(np.abs(z - W))
        if nd != 2 ** j or res > 1e-6: print('  W', W, 'j', j, 'distinct', nd, 'resid %.1e' % res)
        out.append(np.sum(1 / (np.abs(pr) ** 2 * np.abs(dz) ** 2)))
    print('W %-6s' % W, ' '.join('%.2e' % m for m in out)); print('   ratio', ' '.join('%.3f' % (out[i + 1] / out[i]) for i in range(len(out) - 1)), flush=True)
