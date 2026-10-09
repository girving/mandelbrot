"""The kneading partition of f(z) = z² - 3/4 (the parabolic map at the 1/2 root), for limb classification.

ν(z) = 1 iff z lies on the α side of R_{1/6} ∪ [-1/2, 1/2] ∪ R_{2/3} (angles in (1/6, 2/3), the upper-left region):
the partition of every parameter in the limbs near 1/2 converges to it.  A gate passage near -1/2 is a run of 1s, and
the limb denominator is q = 1 + the first run of 1s in the critical orbit's kneading.

  side(z): 1 or 0 (point-in-polygon against the traced rays)."""
import cmath, math

C = -0.75
ER = 1e4

def dyn_ray(theta, depth=40, S=16):
    """Points of the dynamic ray of angle θ for z² + C, outward to inward (|z| ~ 100 down to near the landing point)"""
    pts = []
    z = cmath.rect(math.sqrt(ER), 2 * math.pi * theta)  # level 1: f(z) ≈ z² has radius ER
    for n in range(1, depth + 1):
        ang = 2 * math.pi * ((2**n * theta) % 1)
        for j in range(1, S + 1):
            target = cmath.rect(ER ** (2 ** (-j / S)) if n > 1 or j > 0 else ER, ang)
            # f^n(z) = ER^(2^(-j/S)) e^{i ang}: potential decreasing toward level n+1
            for _ in range(60):
                w, dw = z, 1.0
                for _ in range(n): dw = 2 * w * dw; w = w * w + C
                step = (w - target) / dw
                z -= step
                if abs(step) < 1e-15 * (1 + abs(z)): break
            pts.append(z)
        # next level: f^{n+1}(z) = ER e^{i 2^{n+1} θ 2π} at the same point (continuity)
    return pts

def _polygon():
    r16 = dyn_ray(1 / 6)
    r23 = dyn_ray(2 / 3)
    far = 200.0
    poly = [cmath.rect(far, 2 * math.pi / 6)] + r16 + [0.5, -0.5] + list(reversed(r23))[::-1][::-1]
    # r23 runs outward→inward; we need inward→outward after -1/2
    poly = [cmath.rect(far, 2 * math.pi / 6)] + r16 + [0.5 + 0j, -0.5 + 0j] + list(reversed(r23)) + [cmath.rect(far, 2 * math.pi * 2 / 3)]
    # close through the far arc from 240° back to 60° counterclockwise through 180° (decreasing in angle? the α side
    # holds angles in (60°, 240°), so the arc runs 240° → 180° → 60°)
    for d in range(239, 60, -1): poly.append(cmath.rect(far, math.radians(d)))
    return poly

_POLY = None

def side(z):
    global _POLY
    if _POLY is None: _POLY = _polygon()
    # ray casting in +x direction
    inside = False
    p = _POLY
    for i in range(len(p)):
        a, b = p[i], p[(i + 1) % len(p)]
        if (a.imag > z.imag) != (b.imag > z.imag):
            x = a.real + (z.imag - a.imag) * (b.real - a.real) / (b.imag - a.imag)
            if x > z.real: inside = not inside
    return 1 if inside else 0

if __name__ == '__main__':
    r = dyn_ray(2 / 3)
    print('R_2/3 from', r[0], 'to', r[-1], len(r))
    r = dyn_ray(1 / 6)
    print('R_1/6 from', r[0], 'to', r[-1])
    for z in (-1.0, -0.866, 0.866, 0.34j, -0.34j, -0.6 - 0.1j, -0.4 - 0.05j, 0.2 + 0.3j, 1.5j, -1.5j):
        print(z, side(complex(z)))
