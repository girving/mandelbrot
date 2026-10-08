"""Pick a well-conditioned basis of short prefixes U' and suffixes S' for the shifted blocks H_b[u, s] = F(u·b·s)"""
import os, pickle
import numpy as np
from fractions import Fraction
exec(open('wfa_hp.py').read().split("if __name__")[0])
exec(open('hankel_words.py').read().split("def words")[0])
F = load('all_e.out')
U, S = pickle.load(open(D + 'hankel20_sets.pkl', 'rb'))
H = np.array([[F[val(list(u) + list(s))] for s in S] for u in U])
Uu, Ss, Vt = np.linalg.svd(H, full_matrices=False)
def pivot(M, cand, m):
    """Greedy pivoted QR: m columns of M among cand"""
    R = M[:, cand].copy(); sel = []
    for _ in range(m):
        j = int(np.argmax(np.linalg.norm(R, axis=0))); sel.append(cand[j])
        v = R[:, j] / np.linalg.norm(R[:, j]); R = R - np.outer(v, v @ R); R[:, j] = 0
    return sel
for K, qmax, m in [(44, 12, 50), (44, 14, 56), (44, 16, 60)]:
    P = Uu[:, :K] * Ss[:K]; Q = Vt[:K]
    ucand = [i for i, u in enumerate(U) if cont(list(u)) <= qmax]
    scand = [j for j, s in enumerate(S) if cont(list(s)) <= qmax]
    # Scale columns of P^T / Q by singular values so pivoting favors the small directions equally
    us = pivot((Uu[:, :K]).T, ucand, m); ss = pivot(Vt[:K], scand, m)
    cu = np.linalg.svd(Uu[us, :K], compute_uv=False); cs = np.linalg.svd(Vt[:K][:, ss], compute_uv=False)
    print('K=%d, q ≤ %d candidates (%d prefixes, %d suffixes), %d each: smallest singular value of basis rows %.1e, cols %.1e'
          % (K, qmax, len(ucand), len(scand), m, cu[-1], cs[-1]))
    pickle.dump(([U[i] for i in us], [S[j] for j in ss]), open(D + 'basis_%d_%d.pkl' % (K, qmax), 'wb'))
