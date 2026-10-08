"""High-precision bulb total: exact double-double head (q ≤ Q) + automaton control variate for the tail.

  BULB_DATA=dir python3 wfa_hp_sum.py K [K ...]

Automaton F(w) = αᵀ N(a1)⋯N(a_{k-1}) β(a_k), learned from the Hankel block H[u, s] = F(u·s) (256 prefixes q_u ≤ 20,
1101 suffixes q_s ≤ 60):
  N(1..3): suffix closure, H[u, b·s] = P N(b) Q[:, s];
  N(4..64): shifted blocks F(u·b·s) on a 56 × 56 basis of short words (stage2);
  N(80..256): the same shifted blocks at b = 80, 96, 128, 160, 192, 256 (stage3); cubic interpolation in 1/b between
    measured digits, and a quadratic fit in 1/b over b ≥ 96 beyond 256 (those words weigh ~1e-8 in all);
  β(c): Hankel columns for c ≤ 60, a fit in 1/c beyond (large last digits are analytic in 1/c).
Sum over all words: V(x) = π sin²(πx) α + Σ_b (b+x)^-4 N(b)ᵀ V(1/(b+x)) by Chebyshev collocation, then
S_model = Σ_c c^-4 β(c)ᵀ V(1/c), and S ≈ S_exact(Q) + S_model - S_model(Q).
"""
import os, sys, pickle, math
import numpy as np
from fractions import Fraction
from math import gcd
from numpy.polynomial import chebyshev as C
exec(open('wfa_hp.py').read().split("if __name__")[0])
D = os.environ['BULB_DATA'] + '/'

F = load('all_e.out', 'stage2.out', 'stage3.out')
area = {}
for n in ('all_e.out', 'stage2.out', 'stage3.out'):
    for line in open(D + n):
        t = line.split()
        if t[2] != 'failed':
            p, q = int(t[0]), int(t[1]); area[Fraction(p, q)] = area[Fraction(q - p, q)] = (float(t[4]), float(t[9]))
U, S = pickle.load(open(D + 'hankel20_sets.pkl', 'rb'))
Ub, Sb = pickle.load(open(D + 'basis_44_14.pkl', 'rb'))
bvals = pickle.load(open(D + 'stage2_b.pkl', 'rb')) + pickle.load(open(D + 'stage3_b.pkl', 'rb'))
S0, U0, U1, bvals2, cvals2 = pickle.load(open(D + 'large_sets.pkl', 'rb'))
S0 = list(dict.fromkeys(S0))
H = np.array([[F[val(list(u) + list(s))] for s in S] for u in U])
Uidx = {u: i for i, u in enumerate(U)}; Sidx = {s: i for i, s in enumerate(S)}
Uu, Ss, Vt = np.linalg.svd(H, full_matrices=False)

# Exact head by q, both half planes: Σ_{p<q, gcd=1} area(p/q), accumulated with fsum
QMAX = 1000
head = np.zeros(QMAX + 1)
for q in range(2, QMAX + 1):
    terms = []
    for p in range(1, q):
        if gcd(p, q) != 1: continue
        x = Fraction(min(p, q - p), q)
        terms += list(area[x])
    head[q] = math.fsum(terms)
def S_exact(Q): return math.fsum(head[:Q + 1])

def build(K):
    P = Uu[:, :K] * Ss[:K]; Q = Vt[:K]; Pp = np.linalg.pinv(P)
    N = {}
    for b in (1, 2, 3):
        cols = [(Sidx[s], Sidx[(b,) + s]) for s in S if (b,) + s in Sidx]
        N[b] = Pp @ H[:, [c[1] for c in cols]] @ np.linalg.pinv(Q[:, [c[0] for c in cols]])
    Pb = np.linalg.pinv(P[[Uidx[u] for u in Ub]]); Qb = np.linalg.pinv(Q[:, [Sidx[s] for s in Sb]])
    for b in bvals:
        Hb = np.array([[F[val(list(u) + [b] + list(s))] for s in Sb] for u in Ub])
        N[b] = Pb @ Hb @ Qb
    ks = np.array(sorted(N))
    big = ks[ks >= 96]
    X = np.vstack([(1 / big) ** j for j in range(3)]).T
    Nc = np.linalg.lstsq(X, np.array([N[b].ravel() for b in big]), rcond=None)[0]
    def Nb(b):
        if b in N: return N[b]
        if b < ks[-1]:
            # Cubic Lagrange interpolation in t = 1/b through the 4 nearest measured digits
            i = int(np.searchsorted(ks, b)); lo = max(0, min(i - 2, len(ks) - 4)); pts = ks[lo:lo + 4]
            t = 1 / b; ts = 1 / pts.astype(float); out = 0
            for j in range(4):
                lj = np.prod([(t - ts[m]) / (ts[j] - ts[m]) for m in range(4) if m != j])
                out = out + lj * N[pts[j]]
            return out
        return (np.array([(1 / b) ** j for j in range(3)]) @ Nc).reshape(K, K)
    al = P[Uidx[()]]
    beta = {s[0]: Q[:, Sidx[s]] for s in S if len(s) == 1}
    cs = np.arange(20, 61)
    Xc = np.vstack([(1 / cs) ** j for j in range(9)]).T
    Bc = np.linalg.lstsq(Xc, np.array([beta[c] for c in cs]), rcond=None)[0]
    def Bt(c): return beta[c] if c in beta else np.array([(1 / c) ** j for j in range(9)]) @ Bc
    return al, Nb, Bt

def model_sum(K, n=40, Bmax=2000):
    al, Nb, Bt = build(K)
    xs = (1 - np.cos(np.pi * (np.arange(n) + 0.5) / n)) / 2
    Bmat = C.chebvander(2 * xs - 1, n - 1)
    bs = np.arange(1, Bmax + 1)
    # L = Σ_b N(b)ᵀ ⊗ W_b, W_b[i, j] = (b + x_i)^-4 T_j(2/(b + x_i) - 1), as one matrix product over b
    Wst = np.array([((b + xs) ** -4.0)[:, None] * C.chebvander(2 / (b + xs) - 1, n - 1) for b in bs])
    Nst = np.array([Nb(int(b)).T for b in bs])
    L = (Nst.reshape(len(bs), K * K).T @ Wst.reshape(len(bs), n * n)).reshape(K, K, n, n)
    L = L.transpose(0, 2, 1, 3).reshape(K * n, K * n)
    # Tail b > Bmax: N(b) ≈ N(Bmax), V(1/(b+x)) ≈ V(0): Σ_{b>Bmax} (b+x)^-4 ≈ (Bmax + 1/2 + x)^-3 / 3
    T0 = C.chebvander(np.array([-1.0]), n - 1)[0]
    tail = ((Bmax + 0.5 + xs) ** -3 / 3)[:, None] * T0[None, :]
    L += np.kron(Nb(Bmax).T, tail)
    A = np.kron(np.eye(K), Bmat) - L
    rhs = np.concatenate([al[i] * np.pi * np.sin(np.pi * xs) ** 2 for i in range(K)])
    coef = np.linalg.solve(A, rhs).reshape(K, n)
    cmax = 200000
    c = np.arange(2, cmax + 1)
    Vc = np.array([C.chebval(2 / c - 1, coef[i]) for i in range(K)])
    Bc = np.array([Bt(int(k)) for k in c]).T
    terms = c ** -4.0 * np.einsum('ic,ic->c', Bc, Vc)
    S_model = math.fsum(terms) + (cmax + 0.5) ** -3 / 3 * (Bt(cmax) @ np.array([C.chebval(-1.0, coef[i]) for i in range(K)]))
    # Model partial sums by q
    part = np.zeros(QMAX + 1)
    for q in range(2, QMAX + 1):
        tq = []
        for p in range(1, q):
            if gcd(p, q) != 1: continue
            w = cf(Fraction(p, q)); v = al.copy()
            for a in w[:-1]: v = v @ Nb(a)
            tq.append(np.pi * math.sin(math.pi * p / q) ** 2 / q ** 4 * (v @ Bt(w[-1])))
        part[q] = math.fsum(tq)
    return S_model, part

if __name__ == '__main__':
    for K in map(int, sys.argv[1:]):
        S_model, part = model_sum(K)
        print('K=%d: S_model = %.16f' % (K, S_model))
        for Q in (200, 300, 400, 500, 600, 700, 800, 900, 1000):
            Se, Sm = S_exact(Q), math.fsum(part[:Q + 1])
            print('  Q=%4d: S_exact %.16f  head error %+.1e  estimate %.16f' % (Q, Se, Sm - Se, Se + S_model - Sm))
        sys.stdout.flush()
