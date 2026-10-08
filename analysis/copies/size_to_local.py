"""Convert `limb_families size` output (name δre_hi δre_lo δim_hi δim_lo p ...) to `bulb_batch --local` input
(name p c_re_hi c_re_lo c_im_hi c_im_lo), c = −3/4 + δ by exact two-sums."""
import sys

def two_sum(a, b):
    s = a + b; bb = s - a
    return s, (a - (s - bb)) + (b - bb)

for l in sys.stdin:
    r = l.split()
    rh, rl = two_sum(-0.75, float(r[1])); rl += float(r[2]); rh, rl = two_sum(rh, rl)
    print('%s %s %.17g %.17g %.17g %.17g' % (r[0], r[5], rh, rl, float(r[3]), float(r[4])))
