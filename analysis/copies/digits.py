"""Run-length digits of seahorse family kneading tails, and the digit word ↔ angle key correspondence.

A family's kneading sequence is 1^{2k} 0 t.  XOR each bit of t with its predecessor (t_0 = 0, the head's last bit):
runs of 1s become runs of 0s and alternations (10)^n become runs of 1s.  The run lengths of the XORed word, with its
first symbol, are the family's digits: (s, r1, r2, ...)."""
from family_wfa import kneading

def digits(t):
    d, prev = [], '0'
    for x in t: d.append('1' if x != prev else '0'); prev = x
    runs = []
    for x in d:
        if runs and runs[-1][0] == x: runs[-1][1] += 1
        else: runs.append([x, 1])
    return (runs[0][0] if runs else '-',) + tuple(r for _, r in runs)

def census(keys_path, K0):
    head, pre = '1' * (2 * K0) + '0', '01' * (K0 - 1)
    out = {}
    for l in open(keys_path):
        j, i, lo, hi = l.split(); t = kneading(pre + lo)[len(head):]
        out[digits(t)] = (int(j), int(i), lo, hi, t)
    return out

def insertion(a, b):
    """Positions i with b = a[:i] + (2 bits) + a[i:]"""
    if len(b) != len(a) + 2: return []
    return [i for i in range(len(a) + 1) if b[:i] + b[i + 2:] == a]

def extend(cs, d, pos, n):
    """Angle keys (lo, hi) for the digit word d with d[pos] replaced by n, from census members: the member with
    that digit at 2 or 3 (same parity as n) and its +2 neighbor give a 2-bit insertion, repeated (n - base)/2 times.
    Valid for digits ≥ 2 (a digit +2 from 2 or more is always one insertion in both words)."""
    base = 2 if n % 2 == 0 else 3
    d0 = d[:pos] + (base,) + d[pos + 1:]; d1 = d[:pos] + (base + 2,) + d[pos + 1:]
    if n <= base + 2: return cs[d[:pos] + (n,) + d[pos + 1:]][2:4]
    out = []
    for w0, w1 in ((cs[d0][2], cs[d1][2]), (cs[d0][3], cs[d1][3])):
        i = insertion(w0, w1)[0]; ins = w1[i:i + 2]
        out.append(w0[:i] + ins * ((n - base) // 2) + w0[i:])
    return tuple(out)

def keys_for(cs, d):
    """Angle keys for any digit word d: digits ≥ 4 reduced to 2 or 3 (same parity) give a base word in the census;
    each reduced digit's insertion (position and 2-bit pattern, from the base word and the base word with that digit
    +2) is applied (n - base)/2 times, all at once in the base words."""
    big = [p for p in range(1, len(d)) if d[p] >= 4]
    db = tuple(x if (p == 0 or x < 4) else (2 if x % 2 == 0 else 3) for p, x in enumerate(d))
    out = []
    for w in (2, 3):
        base = cs[db][w]
        ins = []
        for p in big:
            w1 = cs[db[:p] + (db[p] + 2,) + db[p + 1:]][w]
            i = insertion(base, w1)[0]
            ins.append((i, w1[i:i + 2] * ((d[p] - db[p]) // 2)))
        word = base
        for i, s in sorted(ins, reverse=True): word = word[:i] + s + word[i:]
        out.append(word)
    return tuple(out)
