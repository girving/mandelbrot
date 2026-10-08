#!/usr/bin/env python3
"""Summarize family_owner output: who owns the escape-time tail, and whether the family expansion predicts it.

    family_shares.py out.txt [J0 J1]

Coverage: per escape octave j, the tail area times 2^j (flat if the tail is ∝ 1/n) and the cumulative share owned
by gates of period ≤ p, by deeper gates (max_period < p ≤ Qscan), and by none.
Prediction: g_q = owned n·D / A for the cardioid's p/q roots (averaged over octaves J0..J1, with time shifted by
the parent period P), then Σ_r A_r P_r g_q over each period's satellite roots against what they own.  Off-axis
roots' rows count their conjugates too (family_owner doubles the upper half plane), so they are halved here.
"""
import collections
import sys


def parse(path):
    lines = open(path).read().split('\n')
    w = float(lines[1].split()[3].rstrip(';'))
    roots, rows = [], {}
    for L in lines:
        if L.startswith('root'):
            head, cnt = L.split(':', 1)
            t = head.split()
            p, P, q = int(t[2]), int(t[4]), int(t[6])
            c = complex(float(t[8]), float(t[9][:-1]))
            roots.append((p, P, q, c, float(t[11]), list(map(int, cnt.split()))))
        elif L.startswith('deep') or L.startswith('no gate'):
            rows[L.split()[0]] = list(map(int, L.split(':')[1].split()))
    return w, roots, rows['deep'], rows['no']


def main():
    path = sys.argv[1]
    J0, J1 = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (10, 13)
    w, roots, deep, none = parse(path)
    byp = collections.defaultdict(lambda: [0] * len(deep))
    for p, P, q, c, A, cnt in roots:
        for i, x in enumerate(cnt):
            byp[p][i] += x
    print('octave j  tail·2^j  cumulative share by gate period ≤ p (%)                deep   none')
    for i in range(len(deep)):
        total = sum(byp[p][i] for p in byp) + deep[i] + none[i]
        if total < 30:
            break
        cum, shares = 0, []
        for p in sorted(byp):
            cum += byp[p][i]
            shares.append('%d:%.0f' % (p, 100 * cum / total))
        print('  %2d     %.4f   %-52s %3.0f%%   %3.0f%%' % (i + 6, total * w * 2 ** (i + 6), ' '.join(shares),
                                                                                                                 100 * deep[i] / total, 100 * none[i] / total))
    # Owned n·D per root, averaged over octaves J0..J1
    own = []
    for p, P, q, c, A, cnt in roots:
        half = 0.5 if abs(c.imag) > 1e-9 else 1.0
        v = sum(cnt[j - 6] * w * half * 2 ** j for j in range(J0, J1 + 1)) / (J1 - J0 + 1)
        own.append((p, P, q, A, v, sum(cnt[J0 - 6:J1 - 5])))
    gs = collections.defaultdict(list)
    for p, P, q, A, v, k in own:
        if P == 1 and q > 1 and k > 20:
            gs[q].append(v / A)
    g = {q: sum(x) / len(x) for q, x in gs.items()}
    print('g_q from cardioid roots (octaves %d..%d): %s' % (J0, J1, ' '.join('%d:%.4f' % (q, g[q]) for q in sorted(g))))
    pred = collections.defaultdict(lambda: [0, 0])
    for p, P, q, A, v, k in own:
        if P < p and q in g:
            pred[p][0] += v
            pred[p][1] += g[q] * A * P
    print('period  owned n·D  predicted Σ A P g_q  ratio')
    for p in sorted(pred):
        m, pr = pred[p]
        print('  %2d     %.4f     %.4f             %.2f' % (p, m, pr, m / pr))


if __name__ == '__main__':
    main()
