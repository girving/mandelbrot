import os
import numpy as np, pickle
from fractions import Fraction
exec(open('wfa.py').read().split("F=load")[0])
exec(open('hankel_words.py').read().split("def words")[0])
F=load('words.out','hankel.out','hankel20.out','large2.out','large2_retry.out')
U,S=pickle.load(open(os.environ['BULB_DATA'] + '/hankel20_sets.pkl','rb'))
H=np.array([[F[val(list(u)+list(s))] for s in S] for u in U])
Uidx={u:i for i,u in enumerate(U)}; Sidx={s:i for i,s in enumerate(S)}
print("Hankel %dx%d singular values:"%H.shape," ".join("%.1e"%x for x in np.linalg.svd(H,compute_uv=False)[:28]))
Fv=load('fade.out','modes.out')
val_words=[(cf(x),f) for x,f in Fv.items() if len(cf(x))>=4 and all(a<=4 for a in cf(x)[:-1]) and cf(x)[-1]<=60]
print(len(val_words),"held-out long words (interior digits ≤ 4, q up to 40000)")
for qmax in (10,14,20):
    Us=[u for u in U if cont(list(u))<=qmax]; ui={u:i for i,u in enumerate(Us)}
    Uu,Ss,Vt=np.linalg.svd(H[[Uidx[u] for u in Us]],full_matrices=False)
    res=[]
    for K in range(4,31,2):
        P=Uu[:,:K]*Ss[:K]; Q=Vt[:K]; N={}
        for b in range(1,5):
            rows=[(ui[u],ui[u+(b,)]) for u in Us if u+(b,) in ui]
            if len(rows)<K: N=None; break
            N[b]=np.linalg.lstsq(P[[r[0] for r in rows]],P[[r[1] for r in rows]],rcond=None)[0]
        if N is None: break
        al=P[ui[()]]; be={s[0]:Q[:,Sidx[s]] for s in S if len(s)==1}
        e=[]
        for w,f in val_words:
            v=al.copy()
            for a in w[:-1]: v=v@N[a]
            e.append(abs(v@be[w[-1]]-f))
        res.append("K%d:%.1e"%(K,np.median(e)))
    print("prefixes q_u ≤ %2d (%3d rows): median |error| %s"%(qmax,len(Us)," ".join(res)))
