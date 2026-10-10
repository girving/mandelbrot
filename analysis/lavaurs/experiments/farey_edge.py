"""The Farey edge of the q = 2 strip: the bulbs B_m of the limbs t = 1/m (CF [0; 2, k, m], k → ∞), found recursively
as the satellite child of B_{m-1} (an r = m - 1 source) over the big bulb B (the target), at C ≈ 0.24 m^-4 (gate -1).
Prints "m name r n re im C cusp m^4 C" per bulb.

  python3 farey_edge.py single.txt [--mmax 20]"""
import argparse, math, os, subprocess, sys
B = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../build/release/glavaurs')
K = 2.4674011
ap = argparse.ArgumentParser(); ap.add_argument('single'); ap.add_argument('--mmax', type=int, default=20)
ap.add_argument('--out', default=''); a = ap.parse_args()
env = dict(os.environ, GL_SIDE='-1', CHILDREN_STARTS='26')
bulb = None
for l in open(a.single):
    f = l.split()
    if f[0].startswith('sat'):
        s = complex(float(f[1]), float(f[2]))
        if s.imag > -2.16: bulb = (int(f[0][3:]), s, float(f[3]))
nB, sB, CB = bulb
rad = lambda C: math.sqrt(C / (K * math.pi))
cur = (1, nB, sB, CB)
rows = [(1, 'B', 1, nB, sB, CB, 0.97)]
out = open(a.out, 'w') if a.out else None
for m in range(2, a.mmax + 1):
    r0, n0, s0, C0 = cur
    jc = -(n0 + 1)
    line = 'P %d %d %.17g %.17g %.6g %d %.17g %.17g %d %d\n' % (m, n0, s0.real, s0.imag, rad(C0), nB, sB.real, sB.imag, jc - 16, jc + 4)
    kids = []
    for l in subprocess.run([B, '1', '2', 'children', '2'], input=line, capture_output=True, text=True, env=env).stdout.splitlines():
        f = l.split(); j = int(f[0].rsplit('|', 1)[1]); w = float(f[3])
        if w >= 1: continue
        kids.append((CB * 1.5 * w, complex(float(f[1]), float(f[2])), nB - j))   # bulb target: C_nf = 1.5 C
    target = 0.24 * m ** -4
    kids = sorted(kids, key=lambda k: abs(math.log(k[0] / target)))[:12]
    res = subprocess.run([B, '1', '2', 'area', '2'], input=''.join('H%d %d %d %.17g %.17g\n' % (i, m, k[2], k[1].real, k[1].imag) for i, k in enumerate(kids)),
                         capture_output=True, text=True, env=env).stdout
    best = None
    for l in res.splitlines():
        f = l.split()
        if f[3] == 'failed': continue
        C, cusp = float(f[6]), float(f[8]); k = kids[int(f[0][1:])]
        if cusp > 0.5 and (best is None or abs(math.log(C / target)) < abs(math.log(best[0] / target))):
            best = (C, complex(float(f[3]), float(f[4])), k[2], cusp)
    if best is None:
        print('m %d: no satellite child found' % m); break
    C, s, n, cusp = best
    cur = (m, n, s, C)
    rows.append((m, 'B%d' % m, m, n, s, C, cusp))
    msg = '%d B%d %d %d %.17g %.17g %.6e %.3f %.4f' % (m, m, m, n, s.real, s.imag, C, cusp, m ** 4 * C)
    print(msg, flush=True)
    if out: out.write(msg + '\n'); out.flush()
