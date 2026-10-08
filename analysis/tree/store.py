"""A store of satellite-tree bulbs by address, and a scheduler that computes whatever addresses are requested.

An address is a tuple of rotation numbers (Fractions in (0,1)): (r1,) is the cardioid's r1 bulb, (r1, r2) its r2
bulb, and so on.  F(address) is bulb_areas' normalized area of the last bulb (relative to its parent's multiplier
map at the root).  Children of the cardioid come from the one-level data (all_e.out, stage2.out, stage3.out, cusp
coordinates); children of a deeper parent from tree/nodes/<parent>.out, one bulb_areas run per parent.

  store = Store(); store.request(addresses); store.compute(); store.F(address)
Words over the extended alphabet (digits and '#') map to addresses by `parse`; a word with an empty level is
invalid (F = 0)."""
import os, subprocess
from fractions import Fraction

ROOT = os.environ['BULB_DATA'] + '/'
NODES = ROOT + 'tree/nodes/'
BULB = os.path.dirname(os.path.abspath(__file__)) + '/../../build/release/bulb_areas'

def val(w):
    x = Fraction(0)
    for a in reversed(w): x = 1 / (a + x)
    return x

def parse(word):
    """Extended word (tuple of ints and '#') → address, or None if a level is empty or a word evaluates to 1"""
    levels, cur = [], []
    for s in word:
        if s == '#':
            levels.append(tuple(cur)); cur = []
        else:
            cur.append(s)
    levels.append(tuple(cur))
    if any(len(l) == 0 for l in levels): return None
    addr = tuple(val(list(l)) for l in levels)
    if any(x >= 1 for x in addr): return None
    return addr

def mirror(addr):
    return tuple(1 - x for x in addr)

def canonical(addr):
    """Representative of {addr, mirror(addr)}: first rotation number ≤ 1/2"""
    return mirror(addr) if addr and addr[0] > Fraction(1, 2) else addr

def name(addr):
    return '_'.join('%d-%d' % (x.numerator, x.denominator) for x in addr) if addr else 'root'

def read(path):
    out = {}
    for line in open(path):
        t = line.split()
        if t[2] != 'failed':
            lo = len(t) == 11  # Double-double low parts
            f = float(t[5]) + (float(t[10]) if lo else 0.0)
            area = float(t[4]) + (float(t[9]) if lo else 0.0)
            out[Fraction(int(t[0]), int(t[1]))] = (f, t[2], t[3], area)
    return out

class Store:
    def __init__(self):
        os.makedirs(NODES, exist_ok=True)
        self.kids = {}   # parent address → {child: (F, center_re, center_im)}
        root = {}
        for p in [ROOT + n for n in ('all_e.out', 'stage2.out', 'stage3.out')] + [NODES + 'root.out']:
            if os.path.exists(p):
                for x, v in read(p).items():
                    root[x] = v
                    if 1 - x not in root: root[1 - x] = (v[0], v[1], repr(-float(v[2])), v[3])
        self.kids[()] = root
        self.wanted = {}

    def load(self, parent):
        if parent not in self.kids:
            p = NODES + name(parent) + '.out'
            self.kids[parent] = read(p) if os.path.exists(p) else {}
        return self.kids[parent]

    def F(self, addr):
        """F of an address, using conjugate symmetry (1 - r1, 1 - r2, ...) when only the mirror is stored"""
        if addr is None: return 0.0
        for a in (canonical(addr), mirror(canonical(addr))):
            v = self.load(a[:-1]).get(a[-1])
            if v: return v[0]
        return None

    def A(self, addr):
        """Area of the bulb at an address (conjugate symmetry as in F)"""
        for a in (canonical(addr), mirror(canonical(addr))):
            v = self.load(a[:-1]).get(a[-1])
            if v: return v[3]
        return None

    def request(self, addrs):
        for a in addrs:
            if a is None or self.F(a) is not None: continue
            a = canonical(a)
            while a:
                parent = a[:-1]
                if a[-1] in self.load(parent): break
                self.wanted.setdefault(parent, set()).add(a[-1])
                a = parent   # the parent itself must exist (its center)

    def compute(self, threads=10):
        """Run bulb_areas for every parent with wanted children, shallowest first"""
        for parent in sorted(self.wanted, key=len):
            kids = sorted(self.wanted[parent] - set(self.load(parent)))
            if not kids: continue
            c = ('0', '0', '0') if parent == () else self.load(parent[:-1]).get(parent[-1])
            if c is None:  # The parent's own bulb failed
                print('store: no center for parent %s; skipping %d children' % (name(parent), len(kids)))
                continue
            period = 1
            for x in parent: period *= x.denominator
            p = NODES + name(parent)
            with open(p + '.in', 'w') as f:
                for x in kids: print(x.numerator, x.denominator, file=f)
            env = dict(os.environ, BULB_EXP='1', BULB_TOL='1e-9', MANDELBROT_THREADS=str(threads))
            with open(p + '.in') as fin, open(p + '.new', 'w') as fout:
                subprocess.run([BULB, str(period), c[1], c[2], '64'], stdin=fin, stdout=fout, env=env, check=True)
            with open(p + '.out', 'a') as f:
                f.write(open(p + '.new').read())
            os.remove(p + '.new')
            if parent == ():
                for x, v in read(p + '.out').items():
                    self.kids[()][x] = v
            else:
                self.kids.pop(parent, None); self.load(parent)
        self.wanted = {}
