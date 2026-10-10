import os, sys, math
import numpy as np
exec(open('wfa_hp_sum.py').read().split("\nif __name__ == '__main__'")[0])
for n, Bmax in ((40, 2000), (56, 2000), (40, 4000)):
    S_model, part = model_sum(40, n=n, Bmax=Bmax)
    print('K=40 n=%d Bmax=%d: S_model %.16f  estimate(Q=1000) %.16f' % (n, Bmax, S_model, S_exact(1000) + S_model - math.fsum(part)))
