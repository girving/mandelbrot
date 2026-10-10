"""Per-limb masses of a gl_census level (glavaurs classify): limb t = b/m from (m, offset), b = (offset + 1)/q mod m;
satellites listed with their limb (bulbs: b = (n + 1)/q mod m).

  python3 gl_limbsum.py p q gate theta level_file [--threads N]"""
import argparse, os, subprocess, sys
from collections import defaultdict
from fractions import Fraction

BIN = os.path.join(os.path.dirname(__file__), '../../build/release/glavaurs')

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('p', type=int); ap.add_argument('q', type=int); ap.add_argument('gate', type=int)
    ap.add_argument('theta', type=float); ap.add_argument('level')
    ap.add_argument('--threads', type=int, default=4)
    a = ap.parse_args()
    q = a.q
    rows = [l.split() for l in open(a.level)]
    text = ''.join('%s %s %s %s %s\n' % (f[0], f[1], f[2], f[3], f[4]) for f in rows)
    out = subprocess.run([BIN, str(a.p), str(q), 'classify', str(a.theta), str(a.threads)], input=text,
                         capture_output=True, text=True, env=dict(os.environ, GL_SIDE=str(a.gate))).stdout.splitlines()
    mass, cnt, sats, bad = defaultdict(float), defaultdict(int), [], defaultdict(float)
    for f, l in zip(rows, out):
        g = l.split(); C = float(f[5]); sat = f[6] == '1'; n = int(f[2])
        if g[1] == 'bulb':
            m = int(g[2]); b = ((n + 1) // q) % m if (n + 1) % q == 0 else None
            if sat: sats.append((C, 'bulb of %s' % (Fraction(b, m) if b is not None else '?/%d' % m)))
            else: bad['primitive without a 0'] += C
            continue
        if g[1] in ('failed', 'pre'): bad[g[1]] += C; continue
        m, off = int(g[1]), int(g[2])
        if (off + 1) % q: bad['offset %% q'] += C; continue
        b = ((off + 1) // q) % m if m > 1 else 0
        t = Fraction(b, m) if m > 1 else Fraction(0)
        if sat: sats.append((C, 'satellite in %s' % t)); continue
        mass[t] += C; cnt[t] += 1
    for t in sorted(mass, key=lambda t: (t.denominator, t)): print('  limb %-6s %5d components  %.6e' % (t, cnt[t], mass[t]))
    print('  total %.6e' % sum(mass.values()))
    for k, v in bad.items(): print('  unclassified (%s): %.3e' % (k, v))
    for C, s in sorted(sats, reverse=True)[:12]: print('  sat %.4e %s' % (C, s))

if __name__ == '__main__':
    main()
