import os
from fractions import Fraction
import itertools, pickle
def val(w):
    x=Fraction(0)
    for a in reversed(w): x=1/(a+x)
    return x
def cont(w):  # continuant: denominator scale of a digit block
    q0,q1=0,1
    for a in w: q0,q1=q1,a*q1+q0
    return q1
def words(maxq, canonical):
    out=[]
    def rec(w):
        if w and (not canonical or w[-1]>=2): out.append(tuple(w))
        for a in itertools.count(1):
            if cont(w+[a])>maxq: break
            rec(w+[a])
    rec([])
    return out
U=[()]+[u for u in words(14,False)]
S=[s for s in words(60,True)]
pairs=set()
for u in U:
    for s in S: pairs.add(val(list(u)+list(s)))
print(len(U),"prefixes",len(S),"suffixes",len(pairs),"bulbs, max q",max(x.denominator for x in pairs))
pickle.dump((U,S),open(os.environ['BULB_DATA'] + '/hankel_sets.pkl','wb'))
D=os.environ['BULB_DATA'] + '/'
have=set()
for line in open(D+'words.out'):
    t=line.split()
    if t[2]!='failed': have.add(Fraction(int(t[0]),int(t[1]))); have.add(1-Fraction(int(t[0]),int(t[1])))
need=[x for x in pairs if x not in have]
seen=set()
with open(D+'hankel.txt','w') as f:
    for x in need:
        if x in seen or 1-x in seen: continue
        seen.add(x); print(x.numerator,x.denominator,file=f)
print(len(seen),"new bulbs; Σq",sum(x.denominator for x in seen))
