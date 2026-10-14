# The final result (notes/hybrid-area.md, "Direct tail check"): T(2^32) measured directly by run tail-d13 to 2^36,
# combined with the production extrapolation corrected by its backtest, and subtracted from production's A(2^32)
#   python analysis/tail_direct.py
import math

import tree_results
from final_tail import fit_tail, octaves

prod, tail = tree_results.prod(), tree_results.tail()
i32 = prod.K - 1
assert prod.ks[i32] == 2**32 - 8 and tail.ks[18] == 2**32
A32, sA = prod.A(i32), prod.sA(i32)

# Independent samples agree on A(2^32) and on k·D(k) where the runs overlap
tA, tsA = tail.A(18), tail.sA(18)
print(f'A(2^32): production {A32:.13f} ± {sA:.1e}, tail-d13 {tA:.13f} ± {tsA:.1e} '
      f'({(tA - A32) / math.hypot(sA, tsA):+.2f}σ)')
for i in range(12, 22):
    k = tail.ks[i]
    p = f'production {k * prod.D(i):.4f}({k * prod.sD(i):.4f})' if i < i32 else ''
    print(f'  k·D(k) at 2^{round(math.log2(k))}: tail-d13 {k * tail.D(i):.4f}({k * tail.sD(i):.4f})  {p}')

# Direct: octaves 2^32..2^36 measured, and the rest past 2^36 with k·D(k) ≈ f̄ = 1.22..1.27 (the fitted values)
direct = sum(tail.D(i) for i in range(18, 22))
sdirect = math.sqrt(sum(tail.varD(i) for i in range(18, 22)))
beyond = 1.245 * 2.0**-35  # Σ_{j≥36} f̄ / 2^j
Td, sTd = direct + beyond, math.hypot(sdirect, 0.025 * 2.0**-35)
print(f'octaves 2^32..2^36: {direct:.4e} ± {sdirect:.1e}; beyond 2^36 {beyond:.1e}; T(2^32) direct {Td:.3e} ± {sTd:.1e}')

# Extrapolated: the n = 6 fits at 2^32, corrected by the backtest biases (a + b/j underpredicts by ~0.75%, the
# quadratic overpredicts by ~0.2%), which bracket T to ~±5e-12
js, f, sf = octaves(prod)
T1 = fit_tail(js[-6:], f[-6:], sf[-6:], 'a+b/j', 32)[0] * 1.0075
T2 = fit_tail(js[-6:], f[-6:], sf[-6:], 'quadratic', 32)[0] * 0.998
Te, sTe = (T1 + T2) / 2, 5e-12
print(f'T(2^32) extrapolated: corrected fits {T1:.3e} {T2:.3e} → {Te:.3e} ± {sTe:.0e} '
      f'(direct {(Td - Te) / math.hypot(sTd, sTe):+.1f}σ)')

# Combined by inverse variance
w1, w2 = 1 / sTd**2, 1 / sTe**2
T, sT = (w1 * Td + w2 * Te) / (w1 + w2), 1 / math.sqrt(w1 + w2)
mu, err = A32 - T, math.hypot(sA, sT)
print(f'T(2^32) = {T:.4e} ± {sT:.1e}')
print(f'mu = {mu:.13f} ± {err:.2e} (1σ; statistics {sA:.1e}, tail {sT:.1e}), ± {1.96 * err:.1e} at 95%')
