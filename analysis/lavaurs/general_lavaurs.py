"""Lavaurs model at the root of the p/q bulb (prototype, double precision).

w = z - α_0, f(w) = λ w + w² (λ = e^{2πi p/q}); critical point w = -λ/2, critical value v = -λ²/4.  f^q(w) = w + A w^{q+1}
+ … with A = -1/(q a_{-q}) from the Fatou series Φ (general_fatou.py), Φ(f(w)) = Φ(w) + 1/q.  Attracting directions:
A w^q < 0 (q of them, permuted by f), repelling: A w^q > 0.  On each petal L(w) = log(w^q)/q takes a branch continuous
on that petal; the branches are fixed by following f: Φ(f(w)) = Φ(w) + 1/q.

  Phi_a(w): iterate f into an attracting petal (|w| < R0, arg(-A w^q) small), evaluate the series, subtract n/q.
  Psi(ζ, j): the repelling parametrization of petal j: solve Φ(u) = ζ - m on petal j (Re ζ - m ≪ 0) and apply f^{q m}
           (Ψ_j(ζ + 1) = f^q(Ψ_j(ζ))); Ψ_{j+1}(ζ + 1/q) = f(Ψ_j(ζ)) up to the petal labelling."""
import cmath, math
from general_fatou import solve

class Model:
    def __init__(self, p, q, N=10, R0=0.02):
        self.p, self.q = p, q
        self.a, self.beta, self.lam, res = solve(p, q, N)
        self.A = -1 / (q * self.a[-q])
        self.R0 = R0
        lam = self.lam
        self.v = -lam * lam / 4

    def series(self, w, branch):
        """Φ and Φ' by the series, L(w) = (log(w^q) + 2πi branch)/q"""
        q = self.q
        s = self.beta * (cmath.log(w ** q) + 2j * math.pi * branch) / q
        d = self.beta / w
        for j, c in self.a.items():
            s += c * w ** j; d += j * c * w ** (j - 1)
        return s, d

    def petal(self, w, kind):
        """Index of the attracting (kind = -1) or repelling (+1) petal containing direction w, or None if not in a
        sector of half-width 0.6 × (π/q) around it"""
        q = self.q
        t = cmath.phase(self.A * w ** q * kind)   # 0 on the petal's axis
        if abs(t) > 0.6 * math.pi: return None
        # petal index: which of the q axes; axes of A w^q kind > 0: arg w = (2πk - arg(A kind))/q
        base = -cmath.phase(self.A * kind) / q
        k = round((cmath.phase(w) - base) / (2 * math.pi / q)) % q
        return k

    def branch(self, w, kind, k):
        """The branch number so that L is continuous on petal k: fix L's value on the petal axis to log|w| + i arg_axis"""
        q = self.q
        base = -cmath.phase(self.A * kind) / q + 2 * math.pi * k / q   # axis argument of petal k
        # log(w^q) principal has imag in (-π, π]; we want imag(q * (arg w)) with arg w near base
        target = q * (base + ((cmath.phase(w) - base + math.pi) % (2 * math.pi) - math.pi))
        return round((target - cmath.log(w ** q).imag) / (2 * math.pi))

    def phi_a(self, w, max_steps=200000):
        """Attracting coordinate with derivative, and the entering petal index"""
        dw = 1.0 + 0j
        for n in range(max_steps):
            if abs(w) < self.R0:
                k = self.petal(w, -1)
                if k is not None:
                    s, d = self.series(w, self.branch(w, -1, k))
                    return s - n / self.q, d * dw, k
            dw *= self.lam + 2 * w; w = self.lam * w + w * w
            if abs(w) > 10: return None
        return None

    def psi(self, zeta, k):
        """Repelling parametrization of petal k with derivative"""
        q = self.q
        if not (abs(zeta) < 1e6): return None
        m = max(0, math.ceil(zeta.real + max(400.0, 2 * abs(zeta.imag))))
        zl = zeta - m
        # leading term: a_{-q} w^{-q} ≈ zl → w^q = a_{-q}/zl; pick the root in repelling petal k
        base = -cmath.phase(self.A) / q + 2 * math.pi * k / q
        r = (self.a[-q] / zl) ** (1 / q)
        roots = [r * cmath.exp(2j * math.pi * t / q) for t in range(q)]
        w = min(roots, key=lambda u: abs(((cmath.phase(u) - base + math.pi) % (2 * math.pi)) - math.pi))
        br = self.branch(w, 1, k)
        for _ in range(40):
            s, d = self.series(w, br)
            step = (s - zl) / d; w -= step
            if abs(step) < 1e-16 * abs(w): break
        s, d = self.series(w, br)
        dw = 1 / d
        for _ in range(q * m):
            dw *= self.lam + 2 * w; w = self.lam * w + w * w
            if abs(w) > 10: return None
        return w, dw

if __name__ == '__main__':
    import sys
    p, q = (int(x) for x in sys.argv[1:3]) if len(sys.argv) > 2 else (1, 3)
    M = Model(p, q)
    print('p/q = %d/%d: A = %s, v = %s' % (p, q, M.A, M.v))
    # Φ_a(f(w)) - Φ_a(w) = 1/q on orbit points far from 0
    w = M.v
    a0 = M.phi_a(w); a1 = M.phi_a(M.lam * w + w * w)
    print('Φ_a(v) = %s (petal %d);  Φ_a(f(v)) - Φ_a(v) = %s (expect %.6f)' % (a0[0], a0[2], a1[0] - a0[0], 1 / q))
    # Ψ consistency: f^q(Ψ_k(ζ)) = Ψ_k(ζ + 1), and Φ(Ψ_k(ζ)) = ζ (near the petal)
    for k in range(q):
        z = complex(-3, 0.5)
        u0, _ = M.psi(z, k); u1, _ = M.psi(z + 1, k)
        w = u0
        for _ in range(q): w = M.lam * w + w * w
        print('petal %d: |f^q(Ψ(ζ)) - Ψ(ζ+1)| = %.2e, Ψ(ζ) = %s' % (k, abs(w - u1), u0))

def single_centers(M, n, box=((-0.5, 0.5), (-1.0, 2.5)), grid=24, exit_shift=0):
    """Single-transit centers with excursion n: σ with f^n(Ψ_s(ζ0 + σ)) = -λ/2 (s the entering petal), by Newton from a
    grid; returns [(σ, size estimate C = (π²/... ) A_card/|A D|² in σ units)]"""
    e = M.phi_a(M.v)
    z0, dz0, s = e
    s = (s + exit_shift) % M.q
    wc = -M.lam / 2
    out = []
    def resid(sig):
        r = M.psi(z0 + sig, s)
        if r is None: return None
        x, d = r
        for _ in range(n):
            d *= M.lam + 2 * x; x = M.lam * x + x * x
            if abs(x) > 10: return None
        return x - wc, d
    for a in range(grid):
        for b in range(grid):
            sig = complex(box[0][0] + (box[0][1] - box[0][0]) * (a + 0.5) / grid, box[1][0] + (box[1][1] - box[1][0]) * (b + 0.5) / grid)
            ok = False
            for _ in range(40):
                rd = resid(sig)
                if rd is None or rd[1] == 0 or not (abs(rd[1]) < 1e300): break
                r, d = rd
                step = r / d; sig -= step
                if abs(step) < 1e-10: ok = True; break
                if abs(sig) > 10: break
            if not ok: continue
            if any(abs(sig - o[0]) < 1e-8 for o in out): continue
            # normal form: R(w) = f^n(Ψ(Φ_a(f(w)) + σ)) ≈ crit + D δσ + A (w - crit)²; D = ∂/∂σ = (f^n Ψ)', A = D Φ_a'(v)
            rd = resid(sig)
            if rd is None: continue
            _, D = rd
            Aq = D * dz0
            area = (3 * math.pi / 8) / abs(Aq * D) ** 2
            out.append((sig, area))
    return out
