"""Spectral learning with N(b) from suffix closure: H[u, b·s] = P N(b) Q[:, s] for s with b·s in the suffix set"""
import os
import numpy as np, pickle, sys
from fractions import Fraction
exec(open('wfa.py').read().split("F=load")[0])
exec(open('hankel_words.py').read().split("def words")[0])
F=load('words.out','hankel.out','hankel20.out','large2.out','large2_retry.out')
U,S=pickle.load(open(os.environ['BULB_DATA'] + '/hankel20_sets.pkl','rb'))
H=np.array([[F[val(list(u)+list(s))] for s in S] for u in U])
Uidx={u:i for i,u in enumerate(U)}; Sidx={s:i for i,s in enumerate(S)}
Fv=load('fade.out','modes.out')
val_words=[(cf(x),f) for x,f in Fv.items() if len(cf(x))>=4 and all(a<=4 for a in cf(x)[:-1]) and cf(x)[-1]<=60]
Uu,Ss,Vt=np.linalg.svd(H,full_matrices=False)
def learn(K, bmax=4, how='suffix'):
    P=Uu[:,:K]*Ss[:K]; Q=Vt[:K]; N={}
    Pp=np.linalg.pinv(P)
    for b in range(1,bmax+1):
        if how=='suffix':
            cols=[(Sidx[s],Sidx[(b,)+s]) for s in S if (b,)+s in Sidx]
            Qs=Q[:,[c[0] for c in cols]]
            N[b]=Pp@H[:,[c[1] for c in cols]]@np.linalg.pinv(Qs)
        else:
            rows=[(Uidx[u],Uidx[u+(b,)]) for u in U if u+(b,) in Uidx]
            if len(rows)<K: return None
            N[b]=np.linalg.lstsq(P[[r[0] for r in rows]],P[[r[1] for r in rows]],rcond=None)[0]
    al=P[Uidx[()]]; be={s[0]:Q[:,Sidx[s]] for s in S if len(s)==1}
    return al,N,be
def err(model):
    al,N,be=model; e=[]
    for w,f in val_words:
        v=al.copy()
        for a in w[:-1]: v=v@N[a]
        e.append(abs(v@be[w[-1]]-f))
    return np.median(e), np.max(e)
print("%d held-out long words; Hankel %dx%d"%(len(val_words),*H.shape))
for K in range(6,41,2):
    m=learn(K,how='suffix'); r=learn(K,how='rows')
    print("K=%2d  suffix closure: median %.1e max %.1e   %s"%(K,*err(m),("row closure: median %.1e max %.1e"%err(r)) if r else ""))
