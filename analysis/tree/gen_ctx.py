"""Context rows for separator learning: F(r1; s) for every r1 with q1 ≤ 20 and F(r1, r2; s) for a set of depth-2
contexts, over all digit suffixes s (q_s ≤ 40)."""
import os, sys
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from store import Store, val
from gen import S40, rationals
st = Store()
kids = [val(list(s)) for s in S40]
ctx1 = [(r1,) for r1 in rationals(20)]
ctx2 = [(r1, r2) for r1 in (Fraction(1, 2), Fraction(1, 3), Fraction(2, 5), Fraction(1, 4), Fraction(2, 7))
        for r2 in rationals(6)]
st.request([c + (k,) for c in ctx1 + ctx2 for k in kids])
n = sum(len(v) for v in st.wanted.values()); print(len(st.wanted), 'parents,', n, 'bulbs'); sys.stdout.flush()
st.compute()
print('done')
