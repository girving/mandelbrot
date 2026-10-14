# Backtest the tail extrapolation on the production run: cut at 2^c, fit k·D(k) over the n octaves before the
# cut, predict A(2^c) - A(2^32) (the octaves up to 2^32), and compare with the measured value
#   python analysis/tail_backtest.py
import math

import tree_results
from final_tail import MODELS, fit_tail, octaves

r = tree_results.prod()
js, f, sf = octaves(r)
print('cut  n   model      predicted        measured         error   (error / measured)   meas err')
for c in (26, 27, 28, 29, 30):
    i = c - 14
    meas = sum(r.D(t) for t in range(i, r.K - 1))  # A(2^c) - A(2^32)
    # Error of the measured difference: octaves are nearly independent (a sample escapes in one octave)
    sm = math.sqrt(sum(r.varD(t) for t in range(i, r.K - 1)))
    for n in (6, 8, 10):
        for model in MODELS:
            T, _ = fit_tail(js[i - n:i], f[i - n:i], sf[i - n:i], model, c)
            pred = T - fit_tail(js[i - n:i], f[i - n:i], sf[i - n:i], model, 32)[0]  # Octaves c..31 only
            print(f'2^{c} {n:2d} {model:9s} {pred:.6e}  {meas:.6e}  {pred - meas:+.2e}  ({(pred - meas) / meas:+.2%})  {sm:.1e}')
