"""Shared combinatorics for the general-root pipeline (gl_census.py, gl_limbsum.py).

root_angles(p, q): the two parameter angles of the root of the p/q bulb, as fractions: in the period-q cycle of
doubling with rotation number p/q, the unique adjacent pair 1/(2^q - 1) apart (1/3, 2/3; 1/7, 2/7; 1/15, 2/15;
9/31, 10/31).  gate_theta(p, q, gate): the kneading partition's angle for a gate: the lower root angle for gate +1,
the upper for gate -1 (checked at q = 2, 3, 4).
own_half(...): the half of a gate's σ-plane that is its own: the one whose bulb the gate's kneading classifier labels a
bulb (combinatorial; at q = 2, 3 this is also the bulb with excursion q - 1, but not at q = 4)."""
import os, subprocess
from fractions import Fraction

BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../build/release/glavaurs')

def root_angles(p, q):
    N = 2 ** q - 1
    for a in range(1, N):
        orbit, x = [], a
        for _ in range(q):
            orbit.append(x); x = (2 * x) % N
        if x != a or len(set(orbit)) != q: continue
        srt = sorted(orbit); pos = {v: i for i, v in enumerate(srt)}
        if all(pos[(2 * v) % N] == (pos[v] + p) % q for v in srt):
            for i in range(q):
                lo, hi = srt[i], srt[(i + 1) % q]
                if (hi - lo) % N == 1: return Fraction(lo, N), Fraction(hi, N)
    raise ValueError('no rotation %d/%d cycle' % (p, q))

def gate_theta(p, q, gate):
    lo, hi = root_angles(p, q)
    return lo if gate > 0 else hi

def own_half(p, q, gate, bulbs, threads=2):
    """bulbs: [(name, n, σ)] for the two bulbs; returns (cut Im σ, upper?) for the gate's own half"""
    theta = gate_theta(p, q, gate)
    text = ''.join('%s 1 %d %.17g %.17g\n' % (nm, n, s.real, s.imag) for nm, n, s in bulbs)
    out = subprocess.run([BIN, str(p), str(q), 'classify', '%.17g' % float(theta), str(threads)], input=text,
                         capture_output=True, text=True, env=dict(os.environ, GL_SIDE=str(gate))).stdout.split('\n')
    own = [b for b, l in zip(bulbs, out) if l.split()[1:2] == ['bulb']]
    assert len(own) == 1, (bulbs, out)
    cut = (bulbs[0][2].imag + bulbs[1][2].imag) / 2
    return cut, own[0][2].imag > cut

if __name__ == '__main__':
    for p, q in ((1, 2), (1, 3), (2, 5), (1, 4), (3, 7)):
        print(p, q, root_angles(p, q))
