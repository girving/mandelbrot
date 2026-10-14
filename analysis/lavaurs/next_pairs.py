"""Island pairs for the next transit level: sources are the heaviest r-transit components (cusp-classified
lavaurs_area output; primitive, not tuned; satellites optional), targets the limb's bulb and the heaviest single-transit
families.  Output lines "T<σ>~target n_u u_re u_im n_c c_re c_im" for lavaurs_area --island with LAVAURS_R = r + 1.

  python3 next_pairs.py r classified[,…] census tune nsrc ntgt [satellites] > pairs"""
import sys, gzip

def opened(p): return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)

def key(c): return (round(c.real % 1.0, 8) % 1.0, round(c.imag, 8))

if __name__ == '__main__':
    r, nsrc, ntgt = int(sys.argv[1]), int(sys.argv[5]), int(sys.argv[6])
    sats = len(sys.argv) > 7 and sys.argv[7] == 'satellites'
    tuned = set()
    for l in opened(sys.argv[4]):
        f = l.split()
        if '*' in f[0] and f[3] != 'failed': tuned.add(key(complex(float(f[3]), float(f[4]))))
    src = {}
    for p in sys.argv[2].split(','):
        for l in opened(p):
            f = l.split()
            if f[3] == 'failed' or int(f[1]) != r: continue
            c = complex(float(f[3]), float(f[4])); k = key(c)
            if k in tuned or (float(f[9]) >= 1e-8 and not sats): continue
            src.setdefault(k, (float(f[7]), int(f[2]), c))
    tgt = [('bulb', 1, complex(-1.0074583370365449, 0.16135210336429348))]
    fams = []
    for l in opened(sys.argv[3]):
        f = l.split()
        if f[3] != 'failed': fams.append((float(f[7]), f[0], int(f[2]), complex(float(f[3]), float(f[4]))))
    tgt += [(nm, n, c) for C, nm, n, c in sorted(fams, reverse=True)[:ntgt]]
    for C, n, c in sorted(src.values(), key=lambda v: -v[0])[:nsrc]:
        for nm, nc, t in tgt:
            print('T%.7f%+.7f~%s %d %.17g %.17g %d %.17g %.17g' % (c.real, c.imag, nm, n, c.real, c.imag, nc, t.real, t.imag))
