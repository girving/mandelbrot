import os
import numpy as np, pickle
from fractions import Fraction
exec(open('wfa.py').read().split("Fv=load")[0])
exec(open('hankel_words.py').read().split("def words")[0])
K=8
Uu,Ss,Vt=np.linalg.svd(H,full_matrices=False); Q=Vt[:K]
# Columns: greedy pivoted selection of 24 suffixes with q_s ≤ 12, maximizing spanned volume of Q's columns
cand=[j for j,s in enumerate(S) if cont(list(s))<=12]
sel=[]; R=Q[:,cand].copy()
for _ in range(24):
    j=int(np.argmax(np.linalg.norm(R,axis=0))); sel.append(cand[j])
    v=R[:,j]/np.linalg.norm(R[:,j]); R=R-np.outer(v,v@R)
S0=[S[j] for j in sel]
U0=sorted([u for u in U],key=lambda u:(cont(list(u)) if u else 0))[:16]
U1=sorted([u for u in U],key=lambda u:(cont(list(u)) if u else 0))[:24]
bvals=list(range(15,41))+list(range(45,101,5))+[120,140,170,200,240,280,330,400]
cvals=list(range(61,121))+list(range(125,201,5))+list(range(220,401,20))
need=set()
for u in U0:
    for b in bvals:
        for s in S0: need.add(val(list(u)+[b]+list(s)))
for u in U1:
    for c in cvals: need.add(val(list(u)+[c]))
pickle.dump((S0,U0,U1,bvals,cvals),open(os.environ['BULB_DATA'] + '/large_sets.pkl','wb'))
have=set(F)
seen=set()
with open(D+'large2.txt','w') as f:
    for x in need:
        if x in have or x in seen or 1-x in seen: continue
        seen.add(x); print(x.numerator,x.denominator,file=f)
print(len(seen),"bulbs, Σq",sum(x.denominator for x in seen),"max q",max(x.denominator for x in seen))
print("S0:",S0); print("U0:",U0)
