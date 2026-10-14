"""Formal Fatou coordinate of F(w) = -w + w² (z² - 3/4 at the fixed point z = -1/2, multiplier -1): the series
Φ(w) = A/w² + B/w + C L(w) + Σ_{j=1}^N a_j w^j with Φ(F(w)) = Φ(w) + 1/2, where L is log(w) on one petal of each pair
and log(-w) on the other (so L(F(w)) - L(w) = log(1 - w)).  Exact rational coefficients: orders w^-1, w^0 fix A = B
= 1/4, and orders 1..N+1 are a square linear system in (C, a_1..a_N) (coupled in pairs: (C, a_1), (a_2, a_3), ...).

  python3 fatou_series.py N      # prints A, B, C and a_1..a_N"""
import sys
from fractions import Fraction as Q

def mul(a, b, M):
    out = [Q(0)] * (M + 1)
    for i, x in enumerate(a):
        if x == 0: continue
        for j, y in enumerate(b[:M + 1 - i]): out[i + j] += x * y
    return out

def solve(rows, rhs):
    """Exact Gaussian elimination"""
    n = len(rows); A = [r[:] + [b] for r, b in zip(rows, rhs)]
    for i in range(n):
        p = next(r for r in range(i, n) if A[r][i] != 0); A[i], A[p] = A[p], A[i]
        for r in range(n):
            if r != i and A[r][i] != 0:
                f = A[r][i] / A[i][i]; A[r] = [x - f * y for x, y in zip(A[r], A[i])]
    return [A[i][n] / A[i][i] for i in range(n)]

def series(N):
    M = N + 3
    Fw = [Q(0), Q(-1), Q(1)] + [Q(0)] * (M - 2)
    # 1/F^2 - 1/w^2 = (1/w^2)((1-w)^-2 - 1) = Σ_{n≥-1} (n+3) w^n ;  1/F - 1/w = -(1/w)((1-w)^-1 + 1) = -2/w - Σ_{n≥0} w^n
    A = B = Q(1, 4)
    assert A * 2 + B * (-2) == 0 and A * 3 + B * (-1) == Q(1, 2)   # orders w^-1 and w^0
    known = {n: A * (n + 3) + B * (-1) for n in range(1, M + 1)}       # order n ≥ 1
    log1mw = {n: Q(-1, n) for n in range(1, M + 1)}
    Fpow = [[Q(1)] + [Q(0)] * M]
    for j in range(1, N + 1): Fpow.append(mul(Fpow[-1], Fw, M))
    dpow = {j: {n: Fpow[j][n] - (1 if n == j else 0) for n in range(1, M + 1)} for j in range(1, N + 1)}
    rows = [[log1mw[n]] + [dpow[j][n] for j in range(1, N + 1)] for n in range(1, N + 2)]
    rhs = [-known[n] for n in range(1, N + 2)]
    x = solve(rows, rhs)
    return A, B, x[0], x[1:]

if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    A, B, C, a = series(N)
    # consistency at the next order: residual of order N+2 shows the truncation, not inconsistency
    print('A =', A, ' B =', B, ' C =', C)
    for j, x in enumerate(a, 1): print('a_%d = %s  (%.17g)' % (j, x, float(x)))
