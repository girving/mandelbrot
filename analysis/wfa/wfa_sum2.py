import os
import numpy as np, pickle, sys
from fractions import Fraction
from math import gcd
from numpy.polynomial import chebyshev as C
exec(open('wfa.py').read().split("Fv=load")[0])
K=int(sys.argv[1]) if len(sys.argv)>1 else 8
Uu,Ss,Vt=np.linalg.svd(H,full_matrices=False); P=Uu[:,:K]*Ss[:K]; Q=Vt[:K]
N={}
for b in range(1,15):
    rows=[(Uidx[u],Uidx[u+(b,)]) for u in U if u+(b,) in Uidx]
    N[b]=np.linalg.lstsq(P[[r[0] for r in rows]],P[[r[1] for r in rows]],rcond=None)[0]
alpha=P[Uidx[()]]
beta={s[0]:Q[:,Sidx[s]] for s in S if len(s)==1}
# Large digits from targeted bulbs: N(b) for b ≤ 400 from prefixes U0·b on columns S0; β(c) for c ≤ 400 from rows U1
F.update(load('large2.out','large2_retry.out'))
S0,U0,U1,bvals,cvals=pickle.load(open(os.environ['BULB_DATA'] + '/large_sets.pkl','rb'))
S0=list(dict.fromkeys(S0)); Qs=Q[:,[Sidx[s] for s in S0]]; Qp=np.linalg.pinv(Qs)
P0=P[[Uidx[u] for u in U0]]
for b in bvals:
    Pb=np.array([[F[val(list(u)+[b]+list(s))] for s in S0] for u in U0])@Qp
    N[b]=np.linalg.lstsq(P0,Pb,rcond=None)[0]
P1p=np.linalg.pinv(P[[Uidx[u] for u in U1]])
for c in cvals: beta[c]=P1p@np.array([F[val(list(u)+[c])] for u in U1])
def interp(table, k, tail):
    # Interpolate in 1/k between measured digits; beyond the largest, extrapolate with a fit in 1/k
    ks=np.array(sorted(table)); 
    if k in table: return table[k]
    if k<ks[-1]:
        i=np.searchsorted(ks,k); k0,k1=ks[i-1],ks[i]; t=(1/k-1/k0)/(1/k1-1/k0)
        return (1-t)*table[k0]+t*table[k1]
    return tail(k)
bs=np.array([b for b in bvals if b>=100]); X=np.vstack([(1/bs)**j for j in range(3)]).T
Nc=np.linalg.lstsq(X,np.array([N[b].ravel() for b in bs]),rcond=None)[0]
def Nb(b): return interp(N,b,lambda b:(np.array([(1/b)**j for j in range(3)])@Nc).reshape(K,K))
cs=np.array([c for c in cvals if c>=100]); Xc=np.vstack([(1/cs)**j for j in range(4)]).T
Bc=np.linalg.lstsq(Xc,np.array([beta[c] for c in cs]),rcond=None)[0]
def Bt(c): return interp(beta,c,lambda c:np.array([(1/c)**j for j in range(4)])@Bc)
# Model F for any word, and exact sums
def Fmodel(w):
    v=alpha.copy()
    for a in w[:-1]: v=v@Nb(a)
    return v@Bt(w[-1])
def g(x): return np.pi*np.sin(np.pi*x)**2
# Resolvent: V(x) = g(x) α + Σ_b (b+x)^-4 N(b)^T V(1/(b+x)), Chebyshev collocation, K-vector valued
n=40; xs=(1-np.cos(np.pi*(np.arange(n)+0.5)/n))/2; Bmax=20000
Bmat=C.chebvander(2*xs-1,n-1)
L=np.zeros((n*K,n*K))
for b in range(1,Bmax+1):
    w=(b+xs)**-4.0; Vb=C.chebvander(2/(b+xs)-1,n-1); Nt=Nb(b).T
    L+=np.kron(Nt, w[:,None]*Vb)   # block (i,j) of K×K: N^T[i,j] * (w V)
A=np.kron(np.eye(K),Bmat)-L
rhs=np.concatenate([alpha[i]*g(xs) for i in range(K)])
coef=np.linalg.solve(A,rhs).reshape(K,n)
def V(x): return np.array([C.chebval(2*x-1,coef[i]) for i in range(K)])
Smodel=sum(c**-4.0*Bt(c)@V(1/c) for c in range(2,200001))
print("K=%d: S_model (resolvent, b ≤ %d) = %.13f"%(K,Bmax,Smodel))
# Partial sums over q ≤ Qc: exact data and model
Fx=load('words.out','hankel.out')
rows=[]
for q in range(2,201):
    for p in range(1,q):
        if gcd(p,q)!=1: continue
        x=Fraction(p,q); w=cf(x); W=np.pi*np.sin(np.pi*p/q)**2/q**4
        rows.append((q,W*Fx[x],W*Fmodel(w)))
r=np.array(rows)
for Qc in (32,64,100,128,160,200):
    m=r[:,0]<=Qc; Se,Sm=r[m,1].sum(),r[m,2].sum()
    print("  Q=%3d: S_exact(Q) %.13f  S_model(Q) %.13f  head error %.1e   estimate S_exact + S_model - S_model(Q) = %.13f"%(Qc,Se,Sm,Sm-Se,Se+Smodel-Sm))
