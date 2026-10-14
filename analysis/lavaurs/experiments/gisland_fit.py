"""Representability of an island's state over its whole Voronoi cell (glavaurs gisland output): least-squares fits
of G by holomorphic polynomials of degree N in s (y = κ_c + s^2) and in σ - σ_c, on half the samples, with the
max relative error on the other half.  Geometric decay in N = the state is O(n)-dimensional.
  python3 gisland_fit.py gisland.txt"""
import sys, collections
import numpy as np
pts = collections.defaultdict(list)
for l in open(sys.argv[1]):
    f = l.split()
    if len(f) < 7: continue
    pts[f[0]].append([complex(float(f[1]), float(f[2])), complex(float(f[3]), float(f[4])), complex(float(f[5]), float(f[6]))])
for nm, v in pts.items():
    v = np.array(v); s, sg, G = v[:, 0], v[:, 1], v[:, 2]
    if '--check' in sys.argv: G = np.exp(2 * s) / (3 - s)   # sanity: analytic, pole at s = 3
    sc = sg[np.argmin(abs(s))]   # σ near the cell centre
    tr = np.random.default_rng(1).random(len(s)) < 0.5; te = ~tr   # random split (samples come ordered by sheet)
    out = []
    for lab, z in (('s', s), ('σ', sg - sg.mean())):
        zs = z / abs(z).max()
        errs = []
        for N in (4, 8, 16, 24, 32, 48):
            A = np.vander(zs, N + 1, increasing=True)
            c, *_ = np.linalg.lstsq(A[tr], G[tr], rcond=None)
            errs.append(abs(A[te] @ c - G[te]).max() / abs(G).max())
        out.append('%s: %s' % (lab, ' '.join('%.1e' % e for e in errs)))
    print('%-5s  max|G|/min|G| %.1e   N = 4 8 16 24 32 48 ->  %s' % (nm, abs(G).max() / abs(G).min(), '   '.join(out)))
