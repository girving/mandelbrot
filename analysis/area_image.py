# The thread image: the repo logo (rendered at 6x) with μ(M) typeset in Computer Modern on one line, transparent,
# in a mid-luminance blue that reads on both white and black
#   python analysis/area_image.py OUT.png [fontsize [center]]
import io
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

out = sys.argv[1]
plt.rcParams['mathtext.fontset'] = 'cm'
color = '#4a86ff'  # Relative luminance ~0.25: contrast ~6 on black and ~3.5 on white

# The logo rendered at 6x by `./build/release/logo f-k27.npy 1536 64 scratch/logo-6x.png` (transparent background)
canvas = Image.open('scratch/logo-6x.png').convert('RGBA')
W, H = canvas.size


def tex(s, size):
    """Typeset s (mathtext) to a tight transparent image"""
    fig = plt.figure(figsize=(0.01, 0.01))
    fig.text(0, 0, s, fontsize=size, color=color)
    buf = io.BytesIO()
    fig.savefig(buf, dpi=200, transparent=True, bbox_inches='tight', pad_inches=0.02)
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).convert('RGBA')


# One line on the real axis, from inside the period-2 bulb through the pinch into the main cardioid
size = float(sys.argv[2]) if len(sys.argv) > 2 else 22
center = float(sys.argv[3]) if len(sys.argv) > 3 else 0.60  # Horizontal center, as a fraction of the width
t = tex(r'$\mu(M) = 1.506591883653 \pm 4.7 \times 10^{-11}$', size)
x = int(center * W - t.width / 2)
# Keep clear of the bulb's left edge and the cardioid's cusp on the real axis
assert x >= 0.275 * W and x + t.width <= 0.895 * W, (x / W, (x + t.width) / W)
canvas.alpha_composite(t, (x, (H - t.height) // 2))
canvas.save(out)
print(out, canvas.size, 'text', t.size, 'spans', round(x / W, 3), round((x + t.width) / W, 3))
