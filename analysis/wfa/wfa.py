import os
import numpy as np, pickle
from fractions import Fraction
D=os.environ['BULB_DATA'] + '/'
def load(*names):
    F={}
    for n in names:
        for line in open(D+n):
            t=line.split()
            if t[2]=='failed': continue
            p,q=int(t[0]),int(t[1]); F[Fraction(p,q)]=F[Fraction(q-p,q)]=float(t[5])
    return F
def val(w):
    x=Fraction(0)
    for a in reversed(w): x=1/(a+x)
    return x
def cf(x):
    a=[]; p,q=x.numerator,x.denominator
    while q: a.append(p//q); p,q=q,p%q
    return a[1:]
F=load('words.out','hankel.out')
U,S=pickle.load(open(os.environ['BULB_DATA'] + '/hankel_sets.pkl','rb'))
H=np.array([[F[val(list(u)+list(s))] for s in S] for u in U])
sv=np.linalg.svd(H,compute_uv=False)
print("Hankel %dx%d singular values:"%H.shape," ".join("%.1e"%x for x in sv[:16]))
Uidx={u:i for i,u in enumerate(U)}; Sidx={s:i for i,s in enumerate(S)}
Fv=load('fade.out','modes.out','large.out','card.out')
def model(K):
    Uu,Ss,Vt=np.linalg.svd(H,full_matrices=False)
    P=Uu[:,:K]*Ss[:K]; Q=Vt[:K]
    N={}
    for b in range(1,15):
        rows=[(Uidx[u],Uidx[u+(b,)]) for u in U if u+(b,) in Uidx]
        if len(rows)<K: continue
        i0=[r[0] for r in rows]; i1=[r[1] for r in rows]
        N[b]=np.linalg.lstsq(P[i0],P[i1],rcond=None)[0]
    alpha=P[Uidx[()]]
    base={s[0]:Q[:,Sidx[s]] for s in S if len(s)==1}
    def pred(w):
        if w[-1] not in base or any(a not in N for a in w[:-1]): return None
        v=alpha.copy()
        for a in w[:-1]: v=v@N[a]
        return v@base[w[-1]]
    return pred
for K in (4,6,8,10,12):
    pred=model(K)
    errs=[]
    for x,f in Fv.items():
        w=cf(x)
        if len(w)<2: continue
        p=pred(w)
        if p is None: continue
        if (tuple(w[1:]) in Sidx) and (w[0],) in Uidx: continue   # skip words in the Hankel block itself
        errs.append((x.denominator,abs(p-f)))
    e=np.array(errs)
    out=[]
    for lo,hi in [(0,200),(200,1000),(1000,5000),(5000,50000)]:
        m=(e[:,0]>lo)&(e[:,0]<=hi)
        if m.sum(): out.append("q∈(%d,%d] n=%d med %.1e max %.1e"%(lo,hi,m.sum(),np.median(e[m,1]),e[m,1].max()))
    print("K=%2d validation |F_pred - F|: %s"%(K," | ".join(out)))
