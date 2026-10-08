"""Seahorse-valley families of primitive copies at the 1/2 root c = −3/4, in the cardioid limbs k/(2k+1).

  python3 seahorse.py kmax > jobs     (bulb_batch P = 0 jobs, keys S1_k, S2_k, S3_k)

In the limb k/(2k+1) (q = 2k+1, wake words (01)^{k-1}001 / (01)^{k-1}010) the period ≤ 16 catalog shows three
primitive families continuing in k: S1 (period 2k+3, angles (01)^{k-1}00101 / (01)^{k-1}00110, the limb's largest
primitive), S2 (2k+2, (01)^{k-1}0011 / (01)^{k-1}0100) and S3 (2k+4, (01)^k 0001 / (01)^k 0010).  They approach −3/4
through the valley with δ = c + 3/4 ≈ iπ/(2k + σ).  Centers: complex Newton on f_c^p(0) = 0 from a quadratic
extrapolation of 1/δ in k, accepted if the step from the prediction is small against the spacing (continuity)."""
import sys

def center(c, p):
    for _ in range(100):
        z, dz = 0j, 0j
        for i in range(p): dz = 2 * z * dz + 1; z = z * z + c
        step = z / dz; c -= step
        if abs(step) < 1e-16 * abs(c + 0.75): break
    return c

SEEDS = {  # catalog centers (6 digits) for k = 2..4, from the period ≤ 16 roots
    'S1': (3, [-0.530828+0.668289j, -0.651014+0.478030j, -0.695131+0.368354j]),
    'S2': (2, [-0.596892+0.662981j, -0.690943+0.465350j, -0.721002+0.356015j]),
    'S3': (4, [-0.592466+0.621349j, -0.688501+0.447764j, -0.719768+0.346602j]),
}

def family(name, kmax):
    extra, seeds = SEEDS[name]
    out = {}
    for k, c in zip((2, 3, 4), seeds): out[k] = center(c, 2 * k + extra)
    for k in range(5, kmax + 1):
        u = [1 / (out[k - j] + 0.75) for j in (1, 2, 3)]
        pred = 1 / (3 * u[0] - 3 * u[1] + u[2]) - 0.75
        c = center(pred, 2 * k + extra)
        jump = abs(c - pred) / abs(out[k - 1] - out[k - 2])
        if jump > 0.05: sys.exit('%s k=%d: Newton moved %.3f spacings from the prediction' % (name, k, jump))
        out[k] = c
    return out

kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 120
for name in SEEDS:
    for k, c in sorted(family(name, kmax).items()):
        print('%s_%d 0 %.17g %.17g 0 %d' % (name, k, c.real, c.imag, 2 * k + SEEDS[name][0]))
