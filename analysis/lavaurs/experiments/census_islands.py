"""Island statistics on the full strip census (q = 2, gate -1; a pipeline output directory: single.txt, c_r2..c_rL,
c_tuned_r*): every component above --cmin gets its piece itinerary (glavaurs itinerary, ITIN_Z=1), islands are keyed by
(pieces 1..r-1, critical point z_r), and every descendant is assigned to the islands whose key prefixes its itinerary
(the exact partition).  Per island U (level r): ε_U, descendant mass / C_U at levels r+1..r+3, the share of the
atoms one and two transits past the island with ρ = |σ_W - σ_U| / D > 1/2 (D = the post-island chain's closest
approach to v_U = ζ0 + σ_U mod the deck lattice), B90/B99 over the sub-island pieces at transit r+1 (critical pieces
carry their sheet, so the two sibling islands on one 2:1 piece are distinct), and the share
through a non-critical piece.  Globally: per level, the mass whose parent island is found, orphaned (non-critical
piece) or unmatched.

  python3 census_islands.py dir [--levels 6] [--cmin 1e-12] [--islands 300]"""
import argparse, math, os, subprocess, sys, time
from collections import defaultdict
B = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../build/release/glavaurs')
env = dict(os.environ, GL_SIDE='-1', ITIN_Z='1')
TAU = 2j * math.pi * 1.375 / 2
Z0 = complex(1.27578623595, 4.31968989869)
def lat(z): return min(abs((lambda w: w - round(w.real * 2) / 2)(z - k * TAU)) for k in range(-2, 3))
def rkey(z): return (round(z.real % 1, 6) % 1, round(z.imag, 6))
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dir'); ap.add_argument('--levels', type=int, default=6); ap.add_argument('--cmin', type=float, default=1e-12)
    ap.add_argument('--islands', type=int, default=300); ap.add_argument('--threads', type=int, default=8)
    a = ap.parse_args()
    t0 = time.time()
    comps = {}   # name -> (r, n, σ, C, sat)
    for l in open(os.path.join(a.dir, 'single.txt')):
        f = l.split(); sat = f[0].startswith('sat')
        comps['L1_%d' % len(comps)] = (1, int(f[0][3:]) if sat else int(f[0]), complex(float(f[1]), float(f[2])), float(f[3]), int(sat))
    tuned = set()
    for r in range(2, a.levels + 1):
        tf = os.path.join(a.dir, 'c_tuned_r%d.txt' % r)
        if os.path.exists(tf): tuned |= {l.split()[0] for l in open(tf)}
        for l in open(os.path.join(a.dir, 'c_r%d.txt' % r)):
            f = l.split(); C = float(f[5])
            if C > a.cmin and f[0] not in tuned: comps[f[0]] = (r, int(f[2]), complex(float(f[3]), float(f[4])), C, int(f[6]))
    names = list(comps)
    print('%d components above %.0e (tuned removed), %.0f s' % (len(names), a.cmin, time.time() - t0), flush=True)
    text = ''.join('%s %d %.17g %.17g\n' % (nm, comps[nm][0], comps[nm][2].real, comps[nm][2].imag) for nm in names)
    out = subprocess.run([B, '1', '2', 'itinerary', str(a.threads)], input=text, capture_output=True, text=True, env=env).stdout
    it = {}
    for l in out.splitlines():
        f = l.split()
        if 'fail' in f: continue
        pieces, zs, k = [], [], 1
        while k < len(f) and f[k] != 'z':
            pieces.append((f[k + 1], rkey(complex(float(f[k + 2]), float(f[k + 3])))))
            zs.append(complex(float(f[k + 4]), float(f[k + 5]))); k += 6
        zs.append(complex(float(f[k + 1]), float(f[k + 2])))
        it[f[0]] = (pieces, zs)
    print('itineraries: %d of %d (failed mass %.3e of %.3e), %.0f s' % (len(it), len(names), sum(comps[n][3] for n in names if n not in it and not comps[n][4]),
          sum(comps[n][3] for n in names if not comps[n][4]), time.time() - t0), flush=True)
    # island keys
    island = {}
    for nm, (pieces, zs) in it.items():
        r = comps[nm][0]
        island[tuple(pieces) + (('C', rkey(zs[-1])),)] = nm
    # parent assignment per level
    for R in range(2, a.levels + 1):
        found = orph = unm = 0.0
        for nm, (pieces, zs) in it.items():
            r, n, s, C, sat = comps[nm]
            if r != R or sat: continue
            if pieces[-1][0] == 'N': orph += C
            elif tuple(pieces) in island: found += C
            else: unm += C
        tot = found + orph + unm
        print('level %d: mass %.6e: parent island found %.4f%%, orphan (non-critical piece) %.4f%%, parent not in census %.4f%%' % (
            R, tot, 100 * found / tot, 100 * orph / tot, 100 * unm / tot), flush=True)
    # descendants of the heaviest islands
    desc = defaultdict(list)   # island name -> [(descendant name)]
    keys = {}
    for k, nm in island.items(): keys[nm] = k
    for nm, (pieces, zs) in it.items():
        r = comps[nm][0]
        for rr in range(1, r):
            k = tuple(pieces[:rr - 1]) + (('C', pieces[rr - 1][1]),) if pieces[rr - 1][0].startswith('C') else None
            if k is not None and k in island: desc[island[k]].append(nm)
    cand = sorted([nm for nm in island.values() if comps[nm][0] <= a.levels - 2], key=lambda nm: -comps[nm][3])[:a.islands]
    loc = subprocess.run([B, '1', '2', 'locate'], input=''.join('%s %d %.17g %.17g\n' % (nm, comps[nm][0], comps[nm][2].real, comps[nm][2].imag) for nm in cand),
                         capture_output=True, text=True, env=env).stdout
    eps = {}
    for l in loc.splitlines():
        f = l.split()
        if len(f) > 5: eps[f[0]] = 1 / abs(complex(float(f[4]), float(f[5])))
    rows = []
    for nm in cand:
        r, n, sU, CU, satU = comps[nm]
        vU = Z0 + sU
        m_by = defaultdict(float); tot = {1: 0.0, 2: 0.0}; bad = {1: 0.0, 2: 0.0}; badN = 0.0; badby = defaultdict(float)
        for w in desc[nm]:
            R, _, sW, CW, satW = comps[w]
            if satW: continue
            m_by[R - r] += CW
            mm = R - r - 1
            if mm not in (1, 2): continue
            pieces, zs = it[w]
            D = min(lat(z - vU) for z in zs[r + 1:])
            dd = sW - sU; dd = complex((dd.real + 0.5) % 1 - 0.5, dd.imag)
            rho = abs(dd) / D if D > 0 else float('inf')
            tot[mm] += CW
            if rho > 0.5:
                bad[mm] += CW; badby[pieces[r]] += CW
                if pieces[r][0] == 'N': badN += CW
        bm = sorted(badby.values(), reverse=True); acc = 0; b90 = b99 = 0; tb = sum(bm)
        for k, m in enumerate(bm, 1):
            acc += m
            if not b90 and acc >= 0.9 * tb: b90 = k
            if not b99 and acc >= 0.99 * tb: b99 = k
        e = eps.get(nm, float('nan'))
        rows.append((r, e, bad[1] / tot[1] if tot[1] else float('nan'), bad[2] / tot[2] if tot[2] else float('nan'), b90, b99, tot[1], tot[2]))
        print('L%d %-12s C %.3e sat %d eps %.4f  desc/C %s  bad m=1 %.2e m=2 %.2e (via N %.2e)  B90 %d B99 %d of %d' % (
            r, nm, CU, satU, e, ' '.join('%.3g' % (m_by[k] / CU) for k in (1, 2, 3)), rows[-1][2], rows[-1][3],
            badN / (tot[1] + tot[2]) if tot[1] + tot[2] else 0, b90, b99, len(bm)), flush=True)
    print('by ε (islands with atoms one transit past the island): mass-weighted bad share m=1, m=2; mean B90, B99')
    for lo, hi in ((0.3, 2), (0.1, 0.3), (0.03, 0.1), (0.01, 0.03), (0, 0.01)):
        sel = [x for x in rows if lo <= x[1] < hi and x[6] > 0]
        if not sel: continue
        t1 = sum(x[6] for x in sel); b1 = sum(x[2] * x[6] for x in sel)
        s2 = [x for x in sel if x[7] > 0]; t2 = sum(x[7] for x in s2); b2 = sum(x[3] * x[7] for x in s2)
        print('  ε %.2f-%.2f: %3d islands  bad m=1 %.2e  m=2 %s  B90 %.1f B99 %.1f' % (lo, hi, len(sel), b1 / t1,
              '%.2e' % (b2 / t2) if t2 else '-', sum(x[4] for x in sel) / len(sel), sum(x[5] for x in sel) / len(sel)))
    print('done, %.0f s' % (time.time() - t0))
if __name__ == '__main__':
    main()
