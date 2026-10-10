"""Extrapolate C = lim k^4 area for a family from M-side areas (bulb_batch / bulb_batch3 output, keys NAME_k): least
squares in Decimal (100 digits) over C + Σ_{j ≤ m} b_j k^-j (optionally + (log k)/k^j, j = 2..6) for k ≥ K0, against
the model constant Cm.  python3 s1_family_fit.py results/s1_family_m_e3.txt"""
import sys
from decimal import Decimal as D, getcontext
getcontext().prec = 100
rows=[l.split() for l in open(sys.argv[1])]
data=[(int(r[0].split('_')[1]), D(r[5])+D(r[10])+(D(r[12]) if len(r)>12 else 0)) for r in rows]
G=[(k, A*D(k)**4) for k,A in data]
Cm=D('1.72497406349896482647335138261525654565985967371417839746019e-3')
def lsq(M,y):
    n=len(M[0])
    N=[[sum(M[r][i]*M[r][j] for r in range(len(M))) for j in range(n)]+[sum(M[r][i]*y[r] for r in range(len(M)))] for i in range(n)]
    for i in range(n):
        p=max(range(i,n),key=lambda r:abs(N[r][i])); N[i],N[p]=N[p],N[i]
        for r in range(n):
            if r!=i and N[r][i]:
                f=N[r][i]/N[i][i]; N[r]=[a-f*b for a,b in zip(N[r],N[i])]
    return [N[i][n]/N[i][i] for i in range(n)]
for logs in (False, True):
  for K0 in (32, 64, 128):
    for m in (10, 14, 18, 22):
        pts=[(k,g) for k,g in G if k>=K0]
        def row(k):
            kk=D(k); r=[D(1)]+[kk**(-j) for j in range(1,m+1)]
            if logs: r += [kk**(-j)*kk.ln() for j in range(2,7)]
            return r
        M=[row(k) for k,_ in pts]
        if len(pts) < len(M[0]) + 4: continue
        c=lsq(M,[g for _,g in pts])
        res=max(abs(sum(a*b for a,b in zip(r,c))-g)/g for r,(k,g) in zip(M,pts))
        print('logs %d K0 %4d m %2d  C - C_model %+.3e (rel %+.1e)  max rel resid %.1e' % (logs, K0, m, c[0]-Cm, (c[0]-Cm)/Cm, res))
