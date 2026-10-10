"""Single-transit layer at the p/q root from the general model: the backward tree (glavaurs tree; each component once,
as (n mod q, σ mod 1)), exact areas and cusp (glavaurs area) for components with C_est > --exact, satellites (cusp ≥
0.01: double precision; e.g. the limbs' bulbs) dropped.  One gate's σ-plane mod 1 holds the families of both sides of
the root (limbs [CF(p/q), k] and the other side; at q = 2 mirror pairs), so the total is S_0 summed over both sides.

  python3 gl_single.py p q side [--cmin 1e-13] [--exact 1e-10] [--depth 40]
prints the two-sided S_0, the bulbs' constants, and the heaviest components."""
import argparse, os, subprocess, sys

BIN = os.path.join(os.path.dirname(__file__), '../../build/release/glavaurs')

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('p', type=int); ap.add_argument('q', type=int); ap.add_argument('side', type=int)
    ap.add_argument('--cmin', type=float, default=1e-13)
    ap.add_argument('--exact', type=float, default=1e-10)
    ap.add_argument('--depth', type=int, default=40)
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--out', default='')
    a = ap.parse_args()
    env = dict(os.environ, GL_SIDE=str(a.side))
    out = subprocess.run([BIN, str(a.p), str(a.q), 'tree', str(a.cmin), str(a.depth)], capture_output=True, text=True, env=env).stdout
    seen = {}
    for l in out.splitlines():
        f = l.split(); n = int(f[2]); s = complex(float(f[3]), float(f[4])); C = float(f[5])
        k = (n, round(s.real % 1, 6) % 1, round(s.imag, 6))
        if k not in seen or n < seen[k][1]: seen[k] = (C, n, s)
    comps = sorted(seen.values(), key=lambda x: -x[0])
    heavy = [c for c in comps if c[0] > a.exact]
    text = ''.join('H%d 1 %d %.17g %.17g\n' % (i, n, s.real, s.imag) for i, (C, n, s) in enumerate(heavy))
    res = subprocess.run([BIN, str(a.p), str(a.q), 'area', str(a.threads)], input=text, capture_output=True, text=True, env=env).stdout
    exact = {}
    for l in res.splitlines():
        f = l.split()
        if f[3] == 'failed': continue
        exact[int(f[0][1:])] = (float(f[6]), float(f[8]), complex(float(f[3]), float(f[4])), float(f[7]))
    S0 = 0.0; sats = []; failed = 0; est_part = 0.0; conv_err = 0.0
    rows = []
    for i, (C, n, s) in enumerate(comps):
        if i < len(heavy):
            if i not in exact: failed += 1; S0 += C; continue
            Cx, cusp, c, conv = exact[i]
            if cusp >= 0.01: sats.append((Cx, cusp, c)); continue   # double-precision cusp noise ~1e-9..1e-6
            S0 += Cx; conv_err += Cx * abs(conv); rows.append((Cx, n, c, 'exact'))
        else:
            S0 += C; est_part += C; rows.append((C, n, s, 'est'))
    print('%d/%d side %+d: %d distinct single-transit components (C_est > %g), %d exact (%d failed, counted by estimate)' % (
        a.p, a.q, a.side, len(comps), a.cmin, len(heavy), failed))
    print('  satellites: ' + ', '.join('C %.6e (cusp %.2f) @%.4f%+.4fi' % (C, cu, c.real, c.imag) for C, cu, c in sats))
    print('  S_0 (both sides) = %.10e  (estimated light part %.3e, Σ C |conv| %.1e)' % (S0, est_part, conv_err))
    for C, n, c, how in rows[:8]: print('    %.6e n=%d σ %.5f%+.5fi %s' % (C, n, c.real, c.imag, how))
    if a.out:
        with open(a.out, 'w') as fo:
            for C, n, c, how in rows: fo.write('%d %.17g %.17g %.10e %s\n' % (n, c.real, c.imag, C, how))
            for C, cu, c in sats: fo.write('sat %.17g %.17g %.10e %.3f\n' % (c.real, c.imag, C, cu))

if __name__ == '__main__':
    main()
