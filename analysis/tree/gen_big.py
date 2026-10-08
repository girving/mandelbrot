"""Context-rich block for the joint automaton: rows r1#u (all r1 with q1 ≤ 20; u with q_u ≤ 4, u may be empty),
columns digit suffixes s (q_s ≤ 40) and separator futures #s' (q_s' ≤ 20)."""
import os, sys
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from store import Store, val
from gen import S40, rationals, words
st = Store()
U4 = [()] + [u for u in words(4, False) if val(list(u)) != 1]
ctx = rationals(20)
req = []
for r1 in ctx:
    for u in U4:
        for s in S40: req.append((r1, val(list(u) + list(s))))
        if u:
            for s2 in rationals(20): req.append((r1, val(list(u)), s2))
st.request(req)
n = sum(len(v) for v in st.wanted.values()); print(len(st.wanted), 'parents,', n, 'bulbs'); sys.stdout.flush()
st.compute()
print('done')
