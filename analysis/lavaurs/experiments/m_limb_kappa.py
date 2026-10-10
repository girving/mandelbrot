"""κ of M's limbs [0; 2, k, m] (the q = 2 strip's Farey limbs t = 1/m at finite k) by exact enumeration: every NRP of the
wake with extra period ≤ J (limb_families wake: Lavaurs' algorithm on angle words, satellites and tunings excluded),
centres and size estimates (limb_families size), exact areas (bulb_batch) for the --exact heaviest, the rest from the
size estimate calibrated on those; over the limb's bulb (the cardioid's p/q bulb).  Prints per m: p/q, bulb area, Σ NRP
area / bulb (exact part, estimated part), heaviest ratio and its period.

  python3 m_limb_kappa.py --k 3 --mmax 10 --J 12 [--exact 200]"""
import argparse, math, os, subprocess, sys
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../build/release')
ap = argparse.ArgumentParser(); ap.add_argument('--k', type=int, default=3); ap.add_argument('--mmin', type=int, default=2)
ap.add_argument('--mmax', type=int, default=10); ap.add_argument('--J', type=int, default=12); ap.add_argument('--exact', type=int, default=200)
ap.add_argument('--threads', type=int, default=2); a = ap.parse_args()
def run(prog, args, text=''):
    return subprocess.run([os.path.join(D, prog)] + args, input=text, capture_output=True, text=True,
                          env=dict(os.environ, MANDELBROT_THREADS=str(a.threads))).stdout
k = a.k; strip = '01' * (k - 1)
print('m   p/q        bulb area    κ (exact + est)          heaviest/bulb (period)  NRPs', flush=True)
for m in range(a.mmin, a.mmax + 1):
    p, q = m * k + 1, 2 * m * k + m + 2
    g = math.gcd(p, q); p //= g; q //= g
    wake = [l.split() for l in run('limb_families', ['wake', str(a.J), str(p), str(q)]).splitlines() if l.strip()]
    lines = ''.join('W%d %s %s %d\n' % (i, w[1][len(strip):], w[2][len(strip):], k) for i, w in enumerate(wake))
    size = {}
    for l in run('limb_families', ['size', str(a.threads)], lines).splitlines():
        f = l.split(); name = f[0].rsplit('_', 1)[0]
        size[name] = (complex(float(f[1]) - 0.75, float(f[3])), int(f[5]), float(f[6]))
    bulb = run('bulb_batch', ['--N', '64'], 'B 1 0 0 %d %d\n' % (p, q)).split()
    bulb_area = float(bulb[5])
    top = sorted(size.items(), key=lambda kv: -kv[1][2])
    ex = top[:a.exact]
    res = run('bulb_batch', ['--N', '64'], ''.join('%s 0 %.17g %.17g 0 %d\n' % (nm, v[0].real, v[0].imag, v[1]) for nm, v in ex))
    area = {}
    for l in res.splitlines():
        f = l.split()
        if f[3] != 'failed': area[f[0]] = float(f[5])
    # calibrate area / |s|^2 on the exact ones (median)
    ratios = sorted(area[nm] / v[2] for nm, v in ex if nm in area)
    cal = ratios[len(ratios) // 2] if ratios else 0
    exact_sum = sum(area.values())
    est_sum = sum(cal * v[2] for nm, v in top[a.exact:]) + sum(cal * v[2] for nm, v in ex if nm not in area)
    hv = max(area.items(), key=lambda kv: kv[1])
    print('%2d %4d/%-5d %.6e  %.6e + %.2e   %.4e (%d)   %d' % (m, p, q, bulb_area, exact_sum / bulb_area, est_sum / bulb_area,
          hv[1] / bulb_area, size[hv[0]][1], len(size)), flush=True)
