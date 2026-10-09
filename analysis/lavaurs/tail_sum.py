"""The 1/2-root layer's tail Σ_{k≥K} (k^-4-scaled area per strip) from M data and Lavaurs-model constants.

Per strip k (all limbs with rotation numbers in [θ_k, θ_{k+1}), one mirror side), the NRP area is
G(k)/k^4 with G(k) → S_tot.  Pieces:
  - one transit, census words: G_0(k) from the M census (k = 16..64), fitted with the model constant fixed;
  - two-transit labels present at all k = 16..96: G_1(k) from M, same fit;
  - everything else (labels missing some k, island-only two-transit components, r ≥ 3): model constants with the
    labels' relative 1/k coefficient b1/D (an estimate: one transit has b1/S = 0.37, two-transit labels 1.57), with
    a ±100% error on the 1/k term;
  - one-transit words beyond the census (6.0666e-8): C k^-4 with a ±100% error (long digits are absent at small k);
  - unseen transits r ≥ 4: an estimate passed in, with its own error.
The layer's contribution to μ(M) is 2 R Σ (both mirror sides, R the copy factor); here R is left out.

  python3 tail_sum.py census.txt.gz diag_area.txt.gz model_diag_valid.txt rest_mass rest_unseen unseen_err"""
import sys, gzip, numpy as np
from collections import defaultdict

def hzeta(s, K, N=200):
    M = K + N
    t = sum(k**-s for k in range(K, M))
    return t + M**(1-s)/(s-1) + 0.5*M**-s + s*M**(-s-1)/12 - s*(s+1)*(s+2)*M**(-s-3)/720

def fit_tails(ks, G, S, Ks, orders=(6, 7, 8)):
    """Σ_{k≥K} G(k)/k^4 for G = S + Σ b_j k^-j fitted on (ks, G); returns (value, spread over orders, b1)"""
    out = []
    for m in orders:
        A = np.array([[k**-j for j in range(1, m+1)] for k in ks])
        b, *_ = np.linalg.lstsq(A, np.array(G) - S, rcond=None)
        out.append(([S*hzeta(4, K) + sum(b[j-1]*hzeta(4+j, K) for j in range(1, m+1)) for K in Ks], b[0]))
    vals = np.array([o[0] for o in out])
    return vals[-1], np.abs(vals - vals[-1]).max(axis=0), out[-1][1]

if __name__ == '__main__':
    census, diag, valid = sys.argv[1:4]
    rest, unseen, unseen_err = map(float, sys.argv[4:7])
    Ks = (16, 32, 64)
    S0c = 1.893391954412168e-3    # model constants of the census words
    S0_beyond = 6.0666442e-8
    G0 = defaultdict(float)
    for l in gzip.open(census, 'rt'):
        f = l.split()
        if len(f) > 5 and f[3] != 'failed': k = int(f[0].rsplit('_', 1)[1]); G0[k] += float(f[5]) * k**4
    ks0 = sorted(G0)
    t0, e0, b0 = fit_tails(ks0, [G0[k] for k in ks0], S0c, Ks)
    C = {}
    for l in open(valid):
        f = l.split(); C[f[0]] = float(f[7])
    data = defaultdict(dict)
    for l in gzip.open(diag, 'rt'):
        f = l.split()
        if len(f) < 6 or f[3] == 'failed': continue
        nm, k = f[0].rsplit('_', 1); data[nm][int(k)] = float(f[5])
    ks1 = [16, 20, 24, 32, 40, 48, 64, 80, 96]
    full = [nm for nm in data if nm in C and all(k in data[nm] for k in ks1)]
    D = sum(C[nm] for nm in full)
    t1, e1, b1 = fit_tails(ks1, [sum(data[nm][k] for nm in full) * k**4 for k in ks1], D, Ks, orders=(5, 6, 7))
    beta = b1 / D
    labels_partial = sum(C.values()) - D
    other = labels_partial + rest + unseen
    print('one transit: census words S = %.12e (b1/S = %.3f); two-transit labels at all k: D = %.10e (b1/D = %.3f)' % (S0c, b0 / S0c, D, beta))
    print('rest by model constants: partial labels %.4e + given %.4e + unseen %.4e (± %.1e)' % (labels_partial, rest, unseen, unseen_err))
    print('S_tot (this layer, one mirror side, k→∞ constant) = %.6e' % (S0c + S0_beyond + D + other))
    for i, K in enumerate(Ks):
        z4, z5 = hzeta(4, K), hzeta(5, K)
        t_rest = other * (z4 + beta * z5)
        e_rest = other * beta * z5 + unseen_err * z4
        e_beyond = S0_beyond * z4
        tot = t0[i] + t1[i] + S0_beyond * z4 + t_rest
        err = e0[i] + e1[i] + e_beyond + e_rest
        print('K = %2d: Σ_{k≥K} = %.10e ± %.1e   (census %.10e ± %.0e, labels %.10e ± %.0e, rest %.6e ± %.1e, beyond %.1e)'
              % (K, tot, err, t0[i], e0[i], t1[i], e1[i], t_rest, e_rest, e_beyond))
        print('        2 Σ = %.10e ± %.1e  (× R for μ(M); R ≈ μ(M)/A_card = 1.279 gives %.4e)' % (2 * tot, 2 * err, 2 * tot * 1.2788))
