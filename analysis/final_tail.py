# Tail extrapolation for the production run (notes/hybrid-area.md, "Result"): fit k·D(k) over the last octaves
# below 2^32 and sum the fit past 2^32 to estimate T(2^32), the area of exterior points escaping after 2^32 steps
#   python analysis/final_tail.py
import math

import numpy as np

import tree_results

MODELS = {'a+b/j': lambda j: [1, 1 / j], 'quadratic': lambda j: [1, j, j * j]}


def fit_tail(js, f, sf, model, c):
    """Weighted fit of k·D(k) = f over octaves js, summed from octave c on: T(2^c) and its fit error"""
    cols = MODELS[model]
    W = 1 / sf
    X = np.array([cols(j) for j in js])
    co = np.linalg.lstsq(X * W[:, None], f * W, rcond=None)[0]
    cov = np.linalg.inv((X * W[:, None]).T @ (X * W[:, None]))
    J = range(c, 1000)
    T = sum(max(0.0, float(np.dot(co, cols(j)))) / 2.0**j for j in J)
    g = np.array([sum(cols(j)[i] / 2.0**j for j in J) for i in range(len(co))])
    return T, math.sqrt(g @ cov @ g)


def octaves(r):
    """Octaves j, k·D(k) and its errors for k = ks[i] = 2^j below the last k"""
    n = r.K - 1
    js = np.array([round(math.log2(r.ks[i])) for i in range(n)])
    f = np.array([r.ks[i] * r.D(i) for i in range(n)])
    sf = np.array([r.ks[i] * r.sD(i) for i in range(n)])
    return js, f, sf


if __name__ == '__main__':
    r = tree_results.prod()
    A32, sA = r.A(r.K - 1), r.sA(r.K - 1)
    print(f'A(2^32) = {A32:.13f} ± {sA:.2e}')
    js, f, sf = octaves(r)
    print('k·D(k) j=22..31:', ' '.join(f'{x:.4f}' for x in f[-10:]))
    print('       err     :', ' '.join(f'{x:.4f}' for x in sf[-10:]))
    res = []
    for n in (6, 8, 10):
        for model in MODELS:
            T, sT = fit_tail(js[-n:], f[-n:], sf[-n:], model, 32)
            res.append((model, n, T, sT))
            print(f'  {model:9s} last {n:2d} octaves: T = {T:.4e} ± {sT:.1e}')
    Ts = [x[2] for x in res]
    Tm = float(np.median(Ts))
    spread = (max(Ts) - min(Ts)) / 2
    sfit = float(np.median([x[3] for x in res]))
    mu = A32 - Tm
    err = math.sqrt(sA**2 + spread**2 + sfit**2)
    print(f'T(2^32) = {Tm:.4e} (model spread ±{spread:.1e}, fit ±{sfit:.1e})')
    print(f'mu = {mu:.13f} ± {err:.1e}  (stat {sA:.1e}, tail model {spread:.1e}, tail fit {sfit:.1e})')
    print(f'published 1.5065918849: difference {mu - 1.5065918849:+.2e}')
