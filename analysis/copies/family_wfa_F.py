"""WFA for the normalized family areas F = area / (A_card |s|²) at a fixed k, over kneading tails.

  python3 family_size.py census.txt.gz > F.txt;  limb_families keys J > keys  (K0 = census limb)
  python3 family_wfa_F.py F.txt keys K0 k [a b]

h(t) = F_t(k) - 1 (O(1e-3), exactly computed: no extrapolation in k), spectral learning on the block |u| ≤ a, |s| ≤ b,
tested on all tails of length a + b + 2.  Unlike the raw constants (1e-3..1e-30, power laws in run lengths), F is
O(1) and smooth, and the size estimate s carries the orbit's derivatives."""
import sys
import numpy as np
from family_wfa import kneading, words

if __name__ == '__main__':
    Fpath, keys, K0, K = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    a, b = (int(sys.argv[5]), int(sys.argv[6])) if len(sys.argv) > 6 else (5, 6)
    F = {}
    for l in open(Fpath):
        n, f, _ = l.split(); j, i, k = n.split('_')
        if int(k) == K: F[(int(j[1:]), int(i))] = float(f)
    head, pre = '1' * (2 * K0) + '0', '01' * (K0 - 1)
    G = {}
    for l in open(keys):
        j, i, lo, hi = l.split()
        if (int(j), int(i)) in F: G[kneading(pre + lo)[len(head):]] = F[(int(j), int(i))] - 1
    U, V = words(a), words(b)
    missing = [u + s for u in U for s in V for x in ('', '0', '1') if u + x + s not in G]
    if missing: sys.exit('missing tails (failed jobs at this k?): %s' % missing[:3])
    H = np.array([[G[u + s] for s in V] for u in U])
    Hx = {x: np.array([[G[u + x + s] for s in V] for u in U]) for x in '01'}
    hU = np.array([G[u] for u in U])
    Uu, Sv, Vt = np.linalg.svd(H, full_matrices=False)
    print('h = F(k=%d) - 1: range %.2e..%.2e; Hankel %dx%d σ/σ0 = %s' % (K, min(G.values()), max(G.values()), *H.shape,
          ' '.join('%.0e' % x for x in Sv[:30:2] / Sv[0])))
    tests = [w for w in G if len(w) == a + b + 2]
    for r in (2, 4, 8, 12, 16, 20, 24, 32):
        P = Uu[:, :r] * Sv[:r]; Q = Vt[:r]; Pp, Qp = np.linalg.pinv(P), np.linalg.pinv(Q)
        A = {x: Pp @ Hx[x] @ Qp for x in '01'}; al = P[0]; be = Pp @ hU
        def pred(w):
            v = al
            for x in w: v = v @ A[x]
            return v @ be
        e = np.abs([pred(w) - G[w] for w in tests])
        print('  r=%2d: |error in F| on %d tails: median %.1e, 99%% %.1e, max %.1e' % (r, len(tests), np.median(e), np.quantile(e, .99), e.max()))
