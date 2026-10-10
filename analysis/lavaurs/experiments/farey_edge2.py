"""The Farey-edge bulbs B_m (limb t = 1/m of the q = 2 strip, gate -1) by Newton from their predicted positions: in the
representative Re σ ≈ (1 - 1/m) - 5, n = 8m + 1 (checked against the recursive farey_edge.py for m = 2..6), Im σ ≈
-0.167/m² (refined by extrapolating the previous members).  Prints "m Bm m n re im C cusp m^4 C" (the edge.txt format).

  python3 farey_edge2.py [--mmax 40]"""
import argparse, os, subprocess
B = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../build/release/glavaurs')
ap = argparse.ArgumentParser(); ap.add_argument('--mmax', type=int, default=40); a = ap.parse_args()
env = dict(os.environ, GL_SIDE='-1')
dev = []   # (m, σ - (1 - 1/m - 5)) scaled by m²
for m in range(2, a.mmax + 1):
    base = complex(1 - 1 / m - 5, 0)
    if len(dev) >= 2:   # extrapolate m² (σ - base) linearly in 1/m
        (m1, d1), (m2, d2) = dev[-2], dev[-1]
        dm = d2 + (d2 - d1) * (1 / m - 1 / m2) / (1 / m2 - 1 / m1)
    else: dm = complex(0, -0.167)
    g = base + dm / m ** 2
    res = subprocess.run([B, '1', '2', 'area', '1'], input='B%d %d %d %.17g %.17g\n' % (m, m, 8 * m + 1, g.real, g.imag),
                         capture_output=True, text=True, env=env).stdout.split()
    if len(res) < 9 or res[3] == 'failed': print('m %d failed' % m, flush=True); break
    s = complex(float(res[3]), float(res[4])); C = float(res[6]); cusp = float(res[8])
    dev.append((m, (s - base) * m ** 2))
    print('%d B%d %d %d %.17g %.17g %.6e %.3f %.4f' % (m, m, m, 8 * m + 1, s.real, s.imag, C, cusp, m ** 4 * C), flush=True)
