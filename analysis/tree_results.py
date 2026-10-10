# Area estimates from escape_tree result files (one per shard, as written by --output): A(k), the area of points
# not escaping by step k, and D(k), the area escaping between consecutive ks, with their variances
import glob
import math
import os
import re


def _vec(text, name):
    for line in text.split('\n'):
        t = line.split()
        if t and t[0] == name:
            return list(map(int, t[2:]))
    raise ValueError(f'no {name} line')


class Results:
    def __init__(self, dir):
        files = sorted(glob.glob(os.path.join(dir, 'shard-*.txt')))
        texts = [open(f).read() for f in files]
        p = texts[0].split('\n')[1]  # Parameter line
        field = lambda name: re.search(rf'\b{name} (\S+)', p).group(1)
        assert re.search(r'\bbox -2 0.5 0 1.2\b', p), 'cell() assumes the default box'
        self.ks = list(map(int, p[p.index(' ks ') + 4:].split()))
        self.base, self.depth, self.m = int(field('base')), int(field('depth')), int(field('m'))
        assert len(files) == int(texts[0].split('\n')[2].split()[2]), f'{len(files)} shard files in {dir}'
        self.K = len(self.ks)
        ss = int(field('strata')) ** 2
        self.G = self.m // ss  # Groups of ss jittered samples per leaf, for variance estimates
        self._vscale = 1 / (ss * ss * self.G * (self.G - 1))
        sums = {n: [sum(v) for v in zip(*(_vec(t, n) for t in texts))] for n in ('area', 'diff', 'certified')}
        self.area, self.diff, self.cert = sums['area'], sums['diff'], sums['certified']

    def cell(self, d):
        # The default box -2 0.5 0 1.2 (mirrored in the real axis by the factors of 2 below)
        return (2.5 / (self.base << d)) * (1.2 / (self.base << d))

    def A(self, i):
        """Area not escaping by step ks[i]"""
        a = self.cell(self.depth)
        return 2 * (a * self.area[3 * i] / self.m +
                    sum(self.cell(d) * self.cert[d * self.K + i] for d in range(self.depth + 1)))

    def D(self, i):
        """Area escaping in (ks[i], ks[i + 1]]: A(i) - A(i + 1), from per-sample differences"""
        a, K = self.cell(self.depth), self.K
        return 2 * (a * self.diff[3 * i] / self.m +
                    sum(self.cell(d) * (self.cert[d * K + i] - self.cert[d * K + i + 1]) for d in range(self.depth + 1)))

    def _var(self, v, i):
        # Per leaf with group sums c_g: var = a^2 (Σ c_g^2 - (Σ c_g)^2 / G) / (ss^2 G (G - 1)), doubled area
        a = self.cell(self.depth)
        _, q, p = v[3 * i:3 * i + 3]
        return 4 * a * a * (q - p / self.G) * self._vscale

    def varA(self, i):
        return self._var(self.area, i)

    def varD(self, i):
        return self._var(self.diff, i)

    def sA(self, i):
        return math.sqrt(self.varA(i))

    def sD(self, i):
        return math.sqrt(self.varD(i))


HERE = os.path.dirname(os.path.abspath(__file__))
prod = lambda: Results(os.path.join(HERE, 'results', 'prod-d15'))
tail = lambda: Results(os.path.join(HERE, 'results', 'tail-d13'))
