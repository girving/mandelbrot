"""Which limb a general-root Lavaurs component lies in (limb_class.py for any p/q root).

At the root c0 = λ/2 - λ²/4 of the p/q bulb (λ = e^{2πip/q}), a side of the root is a parameter angle θ of the root
(the lower or upper of the bulb's two root angles); its limbs [CF(p/q), k] → p/q from that side.  The kneading
partition is R_{θ/2} ∪ A ∪ R_{(θ+1)/2}: of the two dynamic rays one lands at α (it is in θ's cycle), the other at -α;
A is the lift under Φ_a of a vertical half-line {Φ_a(crit) ± it} through the critical point's basin component, whose
two branches end at α and -α (f is even).  ν(z) = 1 on the side holding the angles (θ/2, (θ+1)/2) (the critical
value's), and points within R_NEAR of α (gate passages) are 1.  The critical orbit's explicit points come from
`glavaurs orbit`; the first 0 fixes the limb: (m, offset) with m the transit it follows and the offset from x_m, as
in limb_class.py; no 0 before the critical point: the component is its limb's bulb (m = r).

  python3 gl_limb.py p q gate theta arcsign [--test]"""
import argparse, cmath, math, os, subprocess, sys

BIN = os.path.join(os.path.dirname(__file__), '../../build/release/glavaurs')
ER = 1e4
R_NEAR = 0.15

class Partition:
    def __init__(self, p, q, gate, theta, arcsign):
        self.p, self.q, self.gate, self.theta = p, q, gate, theta
        self.lam = cmath.exp(2j * math.pi * p / q)
        self.c0 = self.lam / 2 - self.lam ** 2 / 4
        self.alpha = self.lam / 2
        self.env = dict(os.environ, GL_SIDE=str(gate))
        out = subprocess.run([BIN, str(p), str(q), 'arc', str(arcsign)], capture_output=True, text=True, env=self.env).stdout
        br = {1: [], -1: []}
        for l in out.splitlines():
            f = l.split(); br[int(f[1])].append(complex(float(f[2]), float(f[3])) + self.alpha)   # z coordinates
        # which branch ends at α (z = α) and which at -α
        ends = {b: pts[-1] for b, pts in br.items()}
        to_alpha = min(br, key=lambda b: abs(ends[b] - self.alpha))
        self.arc_alpha, self.arc_malpha = br[to_alpha], br[-to_alpha]
        assert abs(ends[to_alpha] - self.alpha) < 0.05 and abs(ends[-to_alpha] + self.alpha) < 0.05, ends
        a1, a2 = theta / 2, (theta + 1) / 2
        r1, r2 = self.dyn_ray(a1), self.dyn_ray(a2)
        # the ray landing at α joins the α branch, the other the -α branch
        if abs(r1[-1] - self.alpha) < abs(r2[-1] - self.alpha): arc1, arc2 = self.arc_alpha, self.arc_malpha
        else: arc1, arc2 = self.arc_malpha, self.arc_alpha
        far = 200.0
        poly = [cmath.rect(far, 2 * math.pi * a1)] + r1 + list(reversed(arc1)) + [0j] + arc2 + list(reversed(r2)) + \
               [cmath.rect(far, 2 * math.pi * a2)]
        # the far arc from a2 back to a1 through the angles in (a1, a2)
        nd = 720
        for d in range(1, nd):
            t = a2 - (a2 - a1) * d / nd
            poly.append(cmath.rect(far, 2 * math.pi * t))
        self.poly = poly

    def dyn_ray(self, theta, depth=40, S=16):
        """Points of the dynamic ray of angle θ for z² + c0, outward to inward"""
        C = self.c0
        pts = []
        z = cmath.rect(math.sqrt(ER), 2 * math.pi * theta)
        for n in range(1, depth + 1):
            ang = 2 * math.pi * ((2 ** n * theta) % 1)
            for j in range(1, S + 1):
                target = cmath.rect(ER ** (2 ** (-j / S)), ang)
                for _ in range(60):
                    w, dw = z, 1.0
                    for _ in range(n): dw = 2 * w * dw; w = w * w + C
                    step = (w - target) / dw
                    z -= step
                    if abs(step) < 1e-15 * (1 + abs(z)): break
                pts.append(z)
        return pts

    def inside(self, z):
        inside = False
        p = self.poly
        for i in range(len(p)):
            a, b = p[i], p[(i + 1) % len(p)]
            if (a.imag > z.imag) != (b.imag > z.imag):
                x = a.real + (z.imag - a.imag) * (b.real - a.real) / (b.imag - a.imag)
                if x > z.real: inside = not inside
        return inside

    def ksym(self, w):
        """ν of the point w (glavaurs coordinates w = z - α)"""
        if abs(w) < R_NEAR: return 1
        return 1 if self.inside(w + self.alpha) else 0

    def classify(self, comps, J=40):
        """comps: [(name, r, n, σ)] -> {name: (m, offset) | ('bulb', r) | None}"""
        text = ''.join('%s %d %d %.17g %.17g\n' % (nm, r, n, s.real, s.imag) for nm, r, n, s in comps)
        out = subprocess.run([BIN, str(self.p), str(self.q), 'orbit', str(J)], input=text, capture_output=True,
                             text=True, env=self.env).stdout
        pts = {}
        for l in out.splitlines():
            f = l.split()
            if f[1] == 'failed': pts[f[0]] = None; continue
            pts.setdefault(f[0], []).append((f[1], int(f[2]), int(f[3]), complex(float(f[4]), float(f[5]))))
        res = {}
        for nm, r, n, s in comps:
            seq = pts.get(nm)
            if seq is None: res[nm] = None; continue
            res[nm] = ('bulb', r)
            for kind, T, off, w in seq:
                if self.ksym(w) == 0:
                    res[nm] = (T, off) if kind != 'e' or T > 0 else ('pre', off)
                    break
        return res

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('p', type=int); ap.add_argument('q', type=int); ap.add_argument('gate', type=int)
    ap.add_argument('theta', type=float); ap.add_argument('arcsign', type=float)
    ap.add_argument('--comps', default='')   # file "name r n re im"
    a = ap.parse_args()
    P = Partition(a.p, a.q, a.gate, a.theta, a.arcsign)
    comps = []
    for l in open(a.comps):
        f = l.split(); comps.append((f[0], int(f[1]), int(f[2]), complex(float(f[3]), float(f[4]))))
    res = P.classify(comps)
    for nm, r, n, s in comps: print('%-28s r=%d n=%d %s' % (nm, r, n, res[nm]))

if __name__ == '__main__':
    main()
