"""Which limb a Lavaurs-model component lies in: the Farey fraction t = b/m of its rotation number in the strip.

At limb index k the critical orbit of an r-transit component (phase σ, final excursion n) is, in the representative
convention, z_1 = v, x_i = z_{1 + (2k+1) i} with x_i = Ψ(p_i) (p_1 = ζ0 + σ, p_{i+1} = Φ_a(x_i) + σ), and
z_P = 0 at P = r(2k+1) + n + 1.  Gate steps are kneading 1s (α side), so the first 0 of the kneading fixes the
first internal-address step q = (position of the first 0), the limb's denominator.  The limbs in the strip have
q = m(2k+1) + 2b for t = b/m (the bulb's limb: m = 1, b = 0), so the first 0 sits at offset 2b - 1 from x_m
(offset -1: the step just before x_1).  No 0 before the critical point: the component is the bulb of its limb.

Kneading sides come from kneading_side.py (R_{1/6} ∪ [-1/2, 1/2] ∪ R_{2/3}); points within R_NEAR of the parabolic
point count as gate steps (α side), the convention the known kneadings confirm.

  classify(sigma, n, r) -> (m, offset of the first 0 from x_m) or ('bulb', None); b = (offset + 1) / 2"""
import cmath
from lavaurs import F, psi, phi_a, phi_a_entry, R_LOC
from kneading_side import side

R_NEAR = 0.15   # |w| below this: deep in the gate (α side)
V = -0.25 + 0j  # critical value in w = z + 1/2

def ksym(w):
    """Kneading symbol of the point w (in w = z + 1/2 coordinates): points near the parabolic point (|w| < R_NEAR) are
    gate steps, α side (the repelling-petal lanes sit below R_{2/3}'s height Im ζ = 2.16 for every component)"""
    if abs(w) < R_NEAR: return 1
    return side(w - 0.5)

def classify(sigma, n, r, J=40):
    """(m, offset) of the first kneading 0 (offset relative to x_m: -j for the exit steps, l for the excursion), or
    ('bulb', None) if there is none before the critical point, or None on failure"""
    e = phi_a_entry(V)
    if e is None: return None
    p, _, petal = e
    p = p + sigma
    for i in range(1, r + 1):
        for j in range(J, 0, -1):  # exit steps before x_i: Ψ_{petal (-1)^j}(p - j/2)
            w, _ = psi(p - j / 2, petal * (-1) ** j)
            if ksym(w) == 0: return (i, -j)
        x, _ = psi(p, petal)
        if ksym(x) == 0: return (i, 0)
        w = x
        for l in range(1, (n if i == r else 10**6) + 1):
            if i < r and abs(w) < R_LOC and abs(w.imag) < 0.7 * abs(w.real): break
            w = F(w)
            if i == r and l == n: break  # the critical point
            if ksym(w) == 0: return (i, l)
            if abs(w) > 10: return None
        if i == r: return ('bulb', None)
        e = phi_a_entry(x)
        if e is None: return None
        q, _, petal = e
        p = q + sigma
    return ('bulb', None)

if __name__ == '__main__':
    tests = [
        ('S1 (bulb limb)', complex(0.17957995250940523, 0.63460217927426621), 1, 1),
        ('B_1 bulb', complex(-1.0074583370365449, 0.16135210336429348), 1, 1),
        ('D2 (bulb limb label)', complex(0.060750994692719239, 0.60275462968096571), 1, 2),
        ('copy B j2_0 (bulb limb)', complex(0.27484965437325415, 0.30741571267975182), 1, 2),
        ('H232 (label)', complex(-1.1995036363, 0.4784921871), 1, 2),
        ('B_2 mediant bulb', complex(-0.5009213918880588, 0.041786807060279318), 1, 2),
        ('copy A j2_0 (mediant)', complex(-1.4412721738, 0.2073657335), 7, 2),
        ('B2~bulb R-8 (mediant)', complex(-1.4451, 0.1581), 9, 3),
        ('B2~j2_0 L-2 (m=3, 25/53)', complex(-0.3016, 0.1020), 3, 3),
        ('B2~j2_0 L-4 (m=3, 26/55)', complex(-0.6449, 0.1013), 5, 3),
        ('B_3 bulb a', complex(-0.33337459396213631, 0.018701548138881138), 1, 3),
        ('B_3 bulb b', complex(-0.66712135122894356, 0.018696181092554307), 3, 3),
        ('B3~bulb L-4 (m=3 heavy)', complex(-0.6400, 0.0709), 5, 4),
        ('B_4 bulb (t=1/4)', complex(-0.24988, 0.01051), 1, 4),
        ('H9386 (bulb limb label)', complex(-0.23531991264558841, 0.35738769346047089), 1, 2),
        ('9/20 top NRP P=29 (m=2)', complex(-0.56394, 0.13584), 1, 3),
        ('? 7.1e-6 (B3 source)', complex(-0.36220, 0.06510), 1, 4),
        ('? 1.40e-5', complex(-1.21357, 0.24598), 7, 4),
    ]
    for name, s, n, r in tests:
        print('%-26s %s' % (name, classify(s, n, r)))

def kneading(sigma, n, r, J=40):
    """The component's kneading as blocks: [(transit i, explicit symbols before the gate (excursion after x_{i-1}),
    explicit symbols after it (exit steps, then x_i))], i = 1..r, plus the final excursion; gate interiors are 1s.
    Returns (blocks, final) with symbols as strings, or None"""
    e = phi_a_entry(V)
    if e is None: return None
    p, _, petal = e
    p = p + sigma
    blocks, pre = [], ''
    for i in range(1, r + 1):
        ex = ''.join(str(ksym(psi(p - j / 2, petal * (-1) ** j)[0])) for j in range(J, 0, -1))
        x, _ = psi(p, petal)
        ex += str(ksym(x))
        blocks.append((pre, ex))
        w, pre = x, ''
        for l in range(1, (n if i == r else 10**6) + 1):
            if i < r and abs(w) < R_LOC and abs(w.imag) < 0.7 * abs(w.real): break
            w = F(w)
            if i == r and l == n: break
            pre += str(ksym(w))
            if abs(w) > 10: return None
        if i == r: return blocks, pre
        e = phi_a_entry(x)
        if e is None: return None
        q, _, petal = e
        p = q + sigma
    return None

def tail(sigma, n, r, J=40):
    """(m, b, tokens): the limb t = b/m and the kneading after the first 0 as k-free tokens: '0' and runs of ones
    '1^(2k g + c)' written 'g:c' (g gate passages inside the run), normalized to the component's own limb index (for
    m > 1 the offset's b is reduced into [1, m), which shifts the index; for m = 1 to b = 0).  Positions: x_j sits at
    (2k+1) j from x_0 = v; block i's excursion symbols are x_{i-1} + s, its exit symbols x_i - J..x_i."""
    kn = kneading(sigma, n, r, J)
    if kn is None: return None
    blocks, final = kn
    seq = []   # (T, s, symbol): position (2k+1) T + s
    for i, (pre, ex) in enumerate(blocks, 1):
        seq += [(i - 1, s, c) for s, c in enumerate(pre, 1)]
        seq += [(i, t - J, c) for t, c in enumerate(ex)]
    seq += [(len(blocks), s, c) for s, c in enumerate(final, 1)]
    first = next((idx for idx, v in enumerate(seq) if v[2] == '0'), None)
    if first is None: return None
    T0, s0, _ = seq[first]
    m, off = (T0, s0) if s0 > 0 else (T0, s0)
    if s0 <= 0: m, off = T0, s0           # exit steps: offset relative to x_{T0}
    else: m, off = T0, s0                 # excursion after x_{T0}
    if m < 1 or (off + 1) % 2: return None
    braw = (off + 1) // 2
    b = braw % m if m > 1 else 0
    if m > 1 and b == 0: return None
    shift = (braw - b) // m if m > 1 else braw
    toks, run = [], None   # run = (start position) of the current 1-run
    def close(run_start, end):
        g = end[0] - run_start[0]
        c = g + end[1] - run_start[1] + 1      # length (2k+1) g + ds + 1 = 2k g + c
        toks.append('%d:%d' % (g, c - 2 * shift * g))
    prev = (T0, s0)
    for T, s_, c in seq[first + 1:]:
        if c == '1':
            if run is None: run = (T, s_)
            last = (T, s_)
        else:
            if run is not None:
                close(run, last); run = None
            # gap between explicit symbols is ones; if the previous explicit symbol was a 0 and this is a 0 in a
            # different block, the gate interior in between is a run
            if prev is not None and (T, s_) != prev and toks and toks[-1] == '0' and T > prev[0]:
                pass
            toks.append('0')
        prev = (T, s_)
    if run is not None: close(run, last)
    return m, b, shift, toks
