"""Match the Lavaurs model's single-transit centers with the M data's family phases σ_M (δ = iπ/(2k + σ_M)), up to
σ_model = a σ_M (or a conj σ_M) + b and the model's period 1/2 in Re σ (one more excursion step)."""
import sys
import numpy as np
from lavaurs import phi_a_entry
from centers import find

zeta0, _, petal_v = phi_a_entry(-0.25 + 0j)
fams = sorted(((float(r[2]), complex(float(r[4]), float(r[5]))) for r in (l.split() for l in open(sys.argv[1]))), reverse=True)[:200]
re, im = np.meshgrid(np.linspace(-0.5, 0.5, 41), np.linspace(-6, 6, 241))
grid = (re + 1j * im).ravel()
best = None
for cross in (False, True):
    petal = -petal_v if cross else petal_v
    cs = []
    for n in range(0, 9):
        s, d = find(n, petal, zeta0, grid)
        cs.extend(s)
    cs = np.array(cs)
    # reduce mod 1/2 in Re into [0, 1/2)
    red = lambda z: (z.real % 0.5) + 1j * z.imag
    cr = red(cs)
    top = sorted(set(np.round(cr, 6)), key=lambda z: z.imag)
    for a, conj in ((0.5, False), (-0.5, False), (0.5, True), (-0.5, True)):
        S1 = fams[0][1]
        for anchor in cr[:4000:7]:
            b = anchor - a * (np.conj(S1) if conj else S1)
            pred = red(a * (np.conj([f[1] for f in fams]) if conj else np.array([f[1] for f in fams])) + b)
            dist = np.array([np.min(np.abs(np.minimum(np.abs(cr - p), np.abs(cr - p + 0.5)))) for p in pred[:40]])
            score = np.sum(dist < 1e-4)
            if best is None or score > best[0]: best = (score, cross, a, conj, b, dist)
print('best: %d of 40 families matched (cross=%s, a=%s, conj=%s, b=%s); median distance %.1e' %
      (best[0], best[1], best[2], best[3], best[4], np.median(best[5])))
