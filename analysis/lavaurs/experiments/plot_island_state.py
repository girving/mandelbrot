import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
S = sys.argv[1] if len(sys.argv) > 1 else './'   # dir with g.txt (gisland), kids.txt
z0 = complex(1.27578623595388763, 4.3196898986859651)
sU = complex(-0.28523807327430817, -0.95845771894112364); sc = complex(-0.28376405598, -0.958068553875)
rows = [l.split() for l in open(S + 'g.txt')]
s = np.array([complex(float(r[1]), float(r[2])) for r in rows]); sg = np.array([complex(float(r[3]), float(r[4])) for r in rows])
G = np.array([complex(float(r[5]), float(r[6])) for r in rows])
Sv = abs(G) ** 2 / abs((sg - sU) * (sg - sc)) ** 2
# consecutive lines are the two sheets of one y
y = s[0::2] ** 2; Ssum = Sv[0::2] + Sv[1::2]
nx, ny, hx, hy = 61, 121, 0.5, 2.1
X = np.round((y.real + hx) / (2 * hx) * (nx - 1)).astype(int); T = np.round((y.imag + hy) / (2 * hy) * (ny - 1)).astype(int)
img = np.full((ny, nx), np.nan); smooth = np.full((ny, nx), np.nan)
img[T, X] = np.log10(Ssum); smooth[T, X] = np.log10(Ssum * abs(y) ** 2)
# children: y = H(x) + σ - z0, relative to κ_c
k = [l.split() for l in open(S + 'kids.txt')]
kc = complex(float(k[-1][5]), float(k[-1][6])) + sc - z0
ky = np.array([complex(float(r[5]), float(r[6])) + complex(float(r[2]), float(r[3])) - z0 - kc for r in k[:-1]])
km = np.array([float(r[4]) for r in k[:-1]])
inc = (abs(ky.real) <= hx) & (abs(ky.imag) <= hy) & (km > 0)
fig, ax = plt.subplots(1, 3, figsize=(16, 9), gridspec_kw={'width_ratios': [1, 1, 1.25]})
ext = [-hx, hx, -hy, hy]
im = ax[0].imshow(img, origin='lower', extent=ext, cmap='viridis', aspect='auto')
ax[0].scatter(ky[inc].real, ky[inc].imag, s=4 + 3 * (np.log10(km[inc]) - np.log10(km[inc]).min()), c='red', edgecolors='white', linewidths=0.3, label='children (targets in the cell)')
ax[0].plot(0, 0, 'w+', ms=14, mew=2); ax[0].annotate('κ_c (fold: |y-κ_c|^-2)', (0.02, 0.05), color='white', fontsize=9)
ax[0].set_title('log10 S(y): the island state over its cell\n(both sheets summed; red: the children, C_W = C_t S(y_t))', fontsize=10)
ax[0].set_xlabel('Re(y - κ_c)'); ax[0].set_ylabel('Im(y - κ_c)'); ax[0].legend(loc='upper right', fontsize=8)
plt.colorbar(im, ax=ax[0], fraction=0.08)
im2 = ax[1].imshow(smooth, origin='lower', extent=ext, cmap='magma', aspect='auto')
ax[1].plot(0, 0, 'c+', ms=14, mew=2)
ax[1].set_title('log10 ( S(y) |y-κ_c|^2 ) = log10 (R0 + |y-κ_c| R1)\nthe smooth part once the fold singularity is divided out', fontsize=10)
ax[1].set_xlabel('Re(y - κ_c)')
plt.colorbar(im2, ax=ax[1], fraction=0.08)
sc3 = ax[2].scatter((sg - sc).real, (sg - sc).imag, c=np.log10(abs(G)), s=2, cmap='coolwarm')
ax[2].plot(0, 0, 'k+', ms=12, mew=2); ax[2].plot((sU - sc).real, (sU - sc).imag, 'kx', ms=9, mew=2)
ax[2].annotate('σ_c (Θ\' = 0)', (0, 0), xytext=(8, 8), textcoords='offset points', fontsize=9)
ax[2].annotate('σ_U (Λ = 0)', ((sU - sc).real, (sU - sc).imag), xytext=(8, -12), textcoords='offset points', fontsize=9)
ax[2].set_title('the island in the σ-plane, coloured by log10 |G|,\nG = (σ-σ_U)(σ-σ_c)/(Θ\'Λ): holomorphic, no poles left', fontsize=10)
ax[2].set_xlabel('Re(σ - σ_c)'); ax[2].set_ylabel('Im(σ - σ_c)'); ax[2].set_aspect('equal')
plt.colorbar(sc3, ax=ax[2], fraction=0.05)
fig.suptitle('Island of the level-1 component F7 (q = 2 strip): the function S_U and its singular/smooth split', fontsize=12)
plt.tight_layout(); plt.savefig(S + 'island_state.png', dpi=110)
print('children in cell', inc.sum(), 'of', len(ky), '; S range', Ssum.min(), Ssum.max())
