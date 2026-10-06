# Splash image of the method: the certified quadtree (interior cells filled, exterior cells outlined) and the Monte
# Carlo samples in the boundary leaves, colored by escape octave, with an inset zoom on the seahorse valley.
# Opaque black background, in several color schemes.  The data comes from two small CPU runs:
#   ks=$(python3 -c 'print(*[2**j for j in range(4, 21)])')
#   ./build/release/escape_tree --box -3.0 1.5 0 4.5 --base 45 --depth 6 --max-iter 1048584 \
#     --dump DATA/wide/leaves.bin --dump-cells DATA/wide/cells.bin $ks
#   ./build/release/escape_tree --box -0.7636 -0.7236 0.1114 0.1514 --base 4 --depth 5 --max-iter 1048584 \
#     --dump DATA/inset/leaves.bin --dump-cells DATA/inset/cells.bin $ks
#   python analysis/splash.py OUT.png [SCHEME [DATA]]   (3840x2160; SCHEMES below, default emerald)
import io
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image, ImageDraw

M, STRATA, SEED, NK = 16, 2, 1, 17   # Samples per leaf, strata, seed, thresholds (2^4 ... 2^20)
SCHEMES = {
    # name: interior fill, interior edges, exterior edges, sample colormap (and its range), text, frame
    'emerald': ((0.10, 0.75, 0.55), (0.90, 1.00, 0.95), (0.50, 0.56, 0.70), ('plasma', 0.15, 1.0), '#4a86ff', (95, 125, 165)),
    'ice': ((0.10, 0.32, 0.80), (0.75, 0.88, 1.00), (0.38, 0.45, 0.58),
            (LinearSegmentedColormap.from_list('ice', ['#2a1b6e', '#2c6fbb', '#5fd3e8', '#ffffff']), 0.0, 1.0), '#7fd8ff',
            (90, 120, 160)),
    'ember': ((0.55, 0.08, 0.12), (1.00, 0.78, 0.62), (0.48, 0.42, 0.40), ('inferno', 0.2, 1.0), '#ffb347', (150, 110, 90)),
    'aurora': ((0.00, 0.52, 0.58), (0.70, 1.00, 0.95), (0.40, 0.46, 0.58),
               (LinearSegmentedColormap.from_list('aurora', ['#1b3a8a', '#7b2cbf', '#f72585', '#ffd6ff']), 0.0, 1.0),
               '#9ff3ff', (90, 130, 150)),
    'gold': ((0.80, 0.60, 0.15), (1.00, 0.93, 0.62), (0.45, 0.45, 0.45),
             (LinearSegmentedColormap.from_list('gold', ['#3a3a3a', '#8a8a8a', '#e8d9a8', '#ffffff']), 0.0, 1.0), '#f2c94c',
             (140, 125, 90)),
}
scheme = sys.argv[2] if len(sys.argv) > 2 else 'emerald'
DATA = sys.argv[3] if len(sys.argv) > 3 else 'scratch/splash'
_i, _ie, _ee, (_cm, CM0, CM1), TEXT, FRAME = SCHEMES[scheme]
INTERIOR = np.array(_i, np.float32)
INTERIOR_EDGE = np.array(_ie, np.float32)
EXTERIOR_EDGE = np.array(_ee, np.float32)
cmap = matplotlib.colormaps[_cm] if isinstance(_cm, str) else _cm


def octave_color(o):
    """Escape octave o (log2 of the step) → RGB: fast escapes violet, slow ones yellow"""
    t = np.clip((np.asarray(o, np.float64) - 3) / (14 - 3), 0, 1)
    return np.asarray(cmap(np.atleast_1d(CM0 + (CM1 - CM0) * t)))[..., :3].reshape(t.shape + (3,)).astype(np.float32)


def octave_alpha(o):
    """Fast escapes are a faint wash; slow ones (the hard part) opaque"""
    return float(np.interp(o, [3, 5, 7, 9, 12], [0.05, 0.12, 0.35, 0.7, 1.0]))


def mix64(z):
    z = z + np.uint64(0x9e3779b97f4a7c15)
    z = (z ^ (z >> np.uint64(30))) * np.uint64(0xbf58476d1ce4e5b9)
    z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94d049bb133111eb)
    return z ^ (z >> np.uint64(31))


def uniform(seed, key, j):
    return (mix64(np.uint64(seed) ^ mix64(np.uint64(2) * key + np.uint64(j))) >> np.uint64(11)).astype(np.float64) * 2.0**-53


def render(prefix, x0, y0, side0, depth, view, px, ss=2, mirror=False, ext_alpha=(0.3, 0.3, 0.28, 0.24, 0.18, 0.12, 0.08),
           dot=(2, 3), lw_fine=1):
    """The tree of a run with base cells of side side0 from corner (x0, y0), over view (vx0, vx1, vy0, vy1) at px
    pixels per unit, supersampled ss×: returns a premultiplied float RGBA image"""
    vx0, vx1, vy0, vy1 = view
    W, H = int(round((vx1 - vx0) * px)), int(round((vy1 - vy0) * px))
    S = px * ss
    img = np.zeros((H * ss, W * ss, 4), np.float32)
    to_px = lambda x, y: ((x - vx0) * S, (vy1 - y) * S)
    signs = (1, -1) if mirror else (1,)

    def blend(i0, i1, j0, j1, rgb, a):
        i0, i1, j0, j1 = max(0, i0), min(img.shape[0], i1), max(0, j0), min(img.shape[1], j1)
        if i0 >= i1 or j0 >= j1: return
        r = img[i0:i1, j0:j1]
        r *= (1 - a)
        r[..., :3] += a * rgb
        r[..., 3] += a

    def rect(x, y, s, fill, fa, edge, ea, lw):
        px0, py1 = to_px(x, y)
        px1, py0 = to_px(x + s, y + s)
        i0, i1, j0, j1 = int(round(py0)), int(round(py1)), int(round(px0)), int(round(px1))
        if fa > 0: blend(i0, i1, j0, j1, fill, fa)
        if ea > 0 and i1 - i0 > 2 * lw and j1 - j0 > 2 * lw:
            blend(i0, i0 + lw, j0, j1, edge, ea); blend(i1 - lw, i1, j0, j1, edge, ea)
            blend(i0 + lw, i1 - lw, j0, j0 + lw, edge, ea); blend(i0 + lw, i1 - lw, j1 - lw, j1, edge, ea)

    # Certified cells: interior (every threshold's mask bit set) filled, exterior outlined
    cells = np.fromfile(prefix + '/cells.bin', dtype=np.int32).reshape(-1, 4)
    full = (1 << NK) - 1
    for d in range(depth + 1):
        side = side0 / (1 << d)
        lw = ss if d <= 3 else lw_fine
        ea = ext_alpha[min(d, len(ext_alpha) - 1)]
        for _, ix, iy, mask in cells[cells[:, 0] == d]:
            x, y = x0 + ix * side, y0 + iy * side
            for sg in signs:
                yy = y if sg > 0 else -y - side
                if int(mask) & 0xffffffff == full: rect(x, yy, side, INTERIOR, 0.6, INTERIOR_EDGE, 0.45, lw)
                else: rect(x, yy, side, None, 0, EXTERIOR_EDGE, ea, lw)

    # Leaf samples, placed as SampleTask::start places them
    rec = np.fromfile(prefix + '/leaves.bin', dtype=np.uint32).reshape(-1, 2 + M)
    with np.errstate(over='ignore'):
        ix = rec[:, 0].astype(np.int64)[:, None]; iy = rec[:, 1].astype(np.int64)[:, None]
        s = np.arange(M, dtype=np.uint64)[None, :]
        key = mix64((ix.astype(np.uint64) & np.uint64(0xffffffff)) | (iy.astype(np.uint64) << np.uint64(32))) + s
        j = (np.arange(M) % (STRATA * STRATA))[None, :]
        leaf = side0 / (1 << depth)
        xs = x0 + (ix + (j % STRATA + uniform(SEED, key, 0)) / STRATA) * leaf
        ys = y0 + (iy + (j // STRATA + uniform(SEED, key, 1)) / STRATA) * leaf
    codes = rec[:, 2:]
    kind = (codes >> 30).ravel()
    octv = np.log2(np.maximum((codes & ((1 << 30) - 1)).astype(np.float64) * 4, 1)).ravel()
    xs, ys = xs.ravel(), ys.ravel()

    def dots(sel, color, alpha, r):
        for sg in signs:
            px_, py_ = to_px(xs[sel], sg * ys[sel])
            keep = (px_ >= r) & (px_ < img.shape[1] - r) & (py_ >= r) & (py_ < img.shape[0] - r)
            pxi, pyi = px_[keep].astype(int), py_[keep].astype(int)
            for dy in range(-r + 1, r):
                for dx in range(-r + 1, r):
                    if dx * dx + dy * dy > (r - 0.5) ** 2: continue
                    q = img[pyi + dy, pxi + dx]
                    q *= (1 - alpha)
                    q[:, :3] += alpha * color
                    q[:, 3] += alpha
                    img[pyi + dy, pxi + dx] = q

    dots(kind != 0, INTERIOR, 0.5, dot[0])
    oi = np.floor(octv).astype(int)
    for o in range(3, 21):
        sel = (kind == 0) & (oi == o)
        if sel.any(): dots(sel, octave_color(o + 0.5), octave_alpha(o + 0.5), dot[0] if o < 9 else dot[1])
    return img.reshape(H, ss, W, ss, 4).mean(axis=(1, 3))


def to_image(prem):
    a = prem[..., 3:4]
    rgb = np.where(a > 0, prem[..., :3] / np.maximum(a, 1e-6), 0)
    return Image.fromarray((np.concatenate([np.clip(rgb, 0, 1), np.clip(a, 0, 1)], -1) * 255).astype(np.uint8), 'RGBA')


def tex(s, size, color=TEXT):
    plt.rcParams['mathtext.fontset'] = 'cm'
    fig = plt.figure(figsize=(0.01, 0.01))
    fig.text(0, 0, s, fontsize=size, color=color)
    buf = io.BytesIO()
    fig.savefig(buf, dpi=200, transparent=True, bbox_inches='tight', pad_inches=0.02)
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).convert('RGBA')



out = sys.argv[1]
W, H = 3840, 2160
VIEW = (-2.92, 1.42, -1.220625, 1.220625)  # 16:9, centered on the set
PX = W / (VIEW[1] - VIEW[0])
main = to_image(render(DATA + '/wide', -3.0, 0.0, 4.5 / 45, 6, VIEW, PX, mirror=True,
                       ext_alpha=(0.3, 0.3, 0.26, 0.2, 0.13, 0.08, 0.05)))
assert main.size == (W, H), main.size
canvas = Image.new('RGBA', (W, H), (0, 0, 0, 255))
canvas.alpha_composite(main)
draw = ImageDraw.Draw(canvas)
frame = FRAME + (230,)

# Inset: the seahorse valley, deeper, in the empty top left (above the antenna, left of the period-2 bulb's filaments)
IB = (-0.7636, -0.7236, 0.1114, 0.1514)
ISIZE = 900
inset = to_image(render(DATA + '/inset', IB[0], IB[2], (IB[1] - IB[0]) / 4, 5, IB, ISIZE / (IB[1] - IB[0]),
                        ext_alpha=(0.35,), dot=(3, 4), lw_fine=2))
ix0, iy0 = 80, 80
sx0, sy0 = (IB[0] - VIEW[0]) * PX, (VIEW[3] - IB[3]) * PX
sx1, sy1 = (IB[1] - VIEW[0]) * PX, (VIEW[3] - IB[2]) * PX
draw.rectangle([sx0, sy0, sx1, sy1], outline=frame, width=4)
draw.line([(ix0 + ISIZE, iy0), (sx1, sy0)], fill=frame, width=3)
draw.line([(ix0 + ISIZE, iy0 + ISIZE), (sx1, sy1)], fill=frame, width=3)
panel = Image.new('RGBA', (ISIZE, ISIZE), (0, 0, 0, 255))
panel.alpha_composite(inset)
canvas.alpha_composite(panel, (ix0, iy0))
draw.rectangle([ix0, iy0, ix0 + ISIZE, iy0 + ISIZE], outline=frame, width=5)

# Legend, bottom right
lx, ly = 3130, 1600
sw, fs = 72, 21
box = lambda x, y, fill, edge: draw.rectangle([x, y, x + sw, y + sw], fill=fill, outline=edge, width=3)
box(lx, ly, tuple(int(255 * c) for c in INTERIOR) + (170,), tuple(int(255 * c) for c in INTERIOR_EDGE) + (140,))
box(lx, ly + 100, None, tuple(int(255 * c) for c in EXTERIOR_EDGE) + (220,))
for i, label in enumerate([r'$\mathrm{proven\ inside}$', r'$\mathrm{proven\ outside}$']):
    lt = tex(label, fs)
    canvas.alpha_composite(lt, (lx + sw + 26, ly + 100 * i + (sw - lt.height) // 2))
gy, gw, gh = ly + 330, 560, 42
for k in range(gw):
    o = 3 + 11 * k / (gw - 1)
    draw.line([(lx + k, gy), (lx + k, gy + gh)],
              fill=tuple(int(255 * v) for v in octave_color(o)) + (int(255 * max(0.35, octave_alpha(o))),))
lt = tex(r'$\mathrm{samples,\ by\ escape\ time}$', fs)
canvas.alpha_composite(lt, (lx, gy - lt.height - 12))
for o, xx in ((3, lx), (14, lx + gw)):
    lt = tex(r'$2^{%d}$' % o, fs)
    canvas.alpha_composite(lt, (xx - lt.width // 2, gy + gh + 10))

# The result, bottom left
t = tex(r'$\mu(M) = 1.506591883653 \pm 4.7 \times 10^{-11}$', 40)
canvas.alpha_composite(t, (100, H - t.height - 110))
canvas.convert('RGB').save(out)
print(out, canvas.size, scheme)
