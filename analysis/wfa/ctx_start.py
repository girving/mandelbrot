"""Does context enter only through the start vector?  F(r1; r2) = α(r1)ᵀ N(r2 digits) β(last), with the cardioid
level's N and β; fit α(r1) per parent on half of its children, test on the rest"""
import os, sys, math
import numpy as np
from fractions import Fraction
exec(open('wfa_hp_sum.py').read().split("\nif __name__ == '__main__'")[0])
T = os.environ['BULB_DATA'] + '/tree/'
K = int(sys.argv[1]) if len(sys.argv) > 1 else 40
al, Nb, Bt = build(K)
par = [tuple(map(int, l.split()[:2])) for l in open(T + 'parents16.txt')]
rng = np.random.default_rng(3)
res_fit, res_test, res_card = [], [], []
for (p1, q1) in par:
    X, y = [], []
    for line in open(T + 'd2/%d-%d.out' % (p1, q1)):
        t = line.split()
        if t[2] == 'failed': continue
        w = cf(Fraction(int(t[0]), int(t[1])))
        v = np.eye(K)
        for a in w[:-1]: v = v @ Nb(a)
        X.append(v @ Bt(w[-1])); y.append(float(t[5]) + float(t[10]))
    X, y = np.array(X), np.array(y)
    tr = rng.random(len(y)) < 0.5
    a1 = np.linalg.lstsq(X[tr], y[tr], rcond=None)[0]
    res_fit.append(np.abs(X[tr] @ a1 - y[tr]).max()); res_test.append(np.abs(X[~tr] @ a1 - y[~tr]).max())
    res_card.append(np.abs(X @ al - y).max())
print('K=%d, %d parents: max |F - α(r1)ᵀNβ| fit median %.1e, held-out median %.1e max %.1e;  with the cardioid α: median %.1e'
      % (K, len(par), np.median(res_fit), np.median(res_test), max(res_test), np.median(res_card)))
