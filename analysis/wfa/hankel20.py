import os
from fractions import Fraction
import pickle
exec(open('hankel_words.py').read().split("U=[()]")[0])
U=[()]+words(20,False); S=words(60,True)
D=os.environ['BULB_DATA'] + '/'
have=set()
for n in ('words.out','hankel.out','large2.out','large2_retry.out'):
    for line in open(D+n):
        t=line.split()
        if t[2]!='failed': x=Fraction(int(t[0]),int(t[1])); have.add(x); have.add(1-x)
need=set()
for u in U:
    for s in S:
        x=val(list(u)+list(s))
        if x not in have and 1-x not in need: need.add(x)
pickle.dump((U,S),open(os.environ['BULB_DATA'] + '/hankel20_sets.pkl','wb'))
with open(D+'hankel20.txt','w') as f:
    for x in need: print(x.numerator,x.denominator,file=f)
print(len(U),"prefixes",len(S),"suffixes;",len(need),"new bulbs, Σq",sum(x.denominator for x in need),"max q",max(x.denominator for x in need))
