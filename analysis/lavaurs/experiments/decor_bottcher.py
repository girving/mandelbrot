"""Test of the renormalization picture of a satellite's decorations: in the island's renormalized coordinate c (U*M = M),
a decoration is a c outside M whose renormalized critical orbit leaves the quadratic-like domain after j returns and
then reaches a target by the outer dynamics, so its exit point is w = Φ_M(c)^(2^j): decorations of the same outer atom
at depths j and j + 1 satisfy Φ_M(c_{j+1})² = Φ_M(c_j).  Reads census_islands' "decorations of satellite" blocks, prints
ζ = Φ_M(c) (by continuation along the ray from |c| = 100, picking the 2^n-th root continuously) and for each decoration
the nearest match of ζ² among the decorations one return (r_U levels) shallower.

  python3 decor_bottcher.py census_islands.log"""
import cmath, math, re, sys
def phi_M(c, prev):
    """Φ_M(c) = φ_c(c) = lim z_n^(1/2^n), the root continuous with prev"""
    z = c
    for n in range(1, 2000):
        z = z * z + c
        if abs(z) > 1e12: break
    else: return None
    L = cmath.log(z); N = 2 ** n
    k = round((N * cmath.phase(prev) - L.imag) / (2 * math.pi)) if prev is not None else 0
    return cmath.exp((L + 2j * math.pi * k) / N)
def phi_path(c, steps=400):
    far = c * (100 / abs(c)); prev = far
    for t in range(1, steps + 1):
        cc = far + (c - far) * (t / steps)
        # geometric spacing toward c
        cc = c + (far - c) * (1 - t / steps) ** 3
        p = phi_M(cc, prev)
        if p is None: return None
        prev = p
    return prev
blocks, cur = [], None
for l in open(sys.argv[1]):
    m = re.match(r'decorations of satellite (\S+) \(L(\d+), C ([0-9.e+-]+)', l)
    if m: cur = (m.group(1), int(m.group(2)), []); blocks.append(cur); continue
    m = re.match(r'\s+c\s+([+-][0-9.]+)\s*([+-][0-9.]+)i\s+area ([0-9.e+-]+)\s+depth (\d+)', l)
    if m and cur: cur[2].append((complex(float(m.group(1)), float(m.group(2))), float(m.group(3)), int(m.group(4))))
for name, r, dec in blocks:
    if not dec: continue
    print('satellite %s (L%d): %d decorations' % (name, r, len(dec)))
    zs = [(c, a, d, phi_path(c)) for c, a, d in dec]
    for c, a, d, z in zs:
        if z is None: print('  c %+.4f%+.4fi depth %d: inside M?' % (c.real, c.imag, d)); continue
        line = '  c %+8.4f%+8.4fi area %.3e depth %d  |ζ| %.4f arg/2π %.4f' % (c.real, c.imag, a, d, abs(z), (cmath.phase(z) / (2 * math.pi)) % 1)
        shallow = [(c2, a2, d2, z2) for c2, a2, d2, z2 in zs if z2 is not None and d2 == d - r]
        if shallow:
            best = min(shallow, key=lambda t: abs(z * z - t[3]))
            line += '   ζ² vs depth-%d ζ: |Δ| %.4f (|ζ²| %.3f, ζ_j %.3f) area ratio %.3f' % (d - r, abs(z * z - best[3]), abs(z * z), abs(best[3]), a / best[1])
        print(line)
