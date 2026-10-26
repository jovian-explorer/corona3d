"""Draws the corona3d logo (docs/logo.svg, docs/logo.png, docs/logo_wordmark.png): the Sun, a
plasma-level surface coloured by Alfven speed, closed coronal loops, and the shell of a coronal
mass ejection rising from the east limb."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch

NAVY = '#0B1B3A'


def draw(ax):
    ax.set_xlim(0, 512); ax.set_ylim(512, 0); ax.set_aspect('equal'); ax.axis('off')
    ax.add_patch(FancyBboxPatch((0, 0), 512, 512, boxstyle='round,pad=0,rounding_size=96', fc=NAVY, ec='none'))
    cx, cy, R = 272, 256, 112
    # plasma-level surface: a ring coloured by Alfven speed (Spectral_r, slow blue to fast red)
    th = np.linspace(0, 2 * np.pi, 361)
    va = 0.5 + 0.5 * np.sin(3 * th + 0.6) * np.cos(th)
    cmap = plt.get_cmap('Spectral_r')
    for a, b, v in zip(th[:-1], th[1:], va[:-1]):
        ra, rb = (R + 30 + 7 * np.sin(4 * t) for t in (a, b))
        ax.plot([cx + ra * np.cos(a), cx + rb * np.cos(b)], [cy + ra * np.sin(a), cy + rb * np.sin(b)],
                color=cmap(v), lw=14, solid_capstyle='round', zorder=2)
    ax.add_patch(Circle((cx, cy), R, fc='#9e9e9e', ec='none', zorder=3))
    # closed loops (arcades) over the west limb, both feet on the surface
    for h, w in ((52, 0.22), (84, 0.34), (116, 0.46)):
        t = np.linspace(-1, 1, 80)
        ang = w * t
        rr = R + h * np.sqrt(1 - t ** 2)
        ax.plot(cx + rr * np.cos(ang), cy + rr * np.sin(ang), color='#EAF0FA', lw=3.5, alpha=0.85, zorder=4)
    # CME shell on the east limb: a round front on two legs, with rings across its axis
    fx, fy, fr = cx - R - 72, cy, 64
    for sgn in (-1, 1):
        foot = (cx + R * np.cos(np.pi - sgn * 0.32), cy + R * np.sin(np.pi - sgn * 0.32))
        ax.plot([foot[0], fx + 6], [foot[1], fy + sgn * fr], color='#d65ae8', lw=6, zorder=5,
                solid_capstyle='round')
    t = np.linspace(np.pi / 2, 3 * np.pi / 2, 100)
    ax.plot(fx + 6 + (fr + 2) * np.cos(t) * 1.05, fy + fr * np.sin(t), color='#d65ae8', lw=6, zorder=5)
    for k in np.linspace(-0.75, 0.75, 6):
        h = fr * np.sqrt(1 - k * k)
        ax.add_patch(matplotlib.patches.Ellipse((fx + 6 + k * fr * 0.9, fy), 14, 2 * h, fc='none', ec='#d65ae8',
                                                lw=2, alpha=0.75, zorder=5))

fig = plt.figure(figsize=(5.12, 5.12), dpi=100)
ax = fig.add_axes([0, 0, 1, 1]); draw(ax)
fig.savefig('docs/logo.png', dpi=100, transparent=True); fig.savefig('docs/logo.svg', transparent=True); plt.close(fig)

fig = plt.figure(figsize=(12.8, 4.0), dpi=100)
ax = fig.add_axes([0.0, 0.0, 0.3125, 1.0]); draw(ax)
ax2 = fig.add_axes([0.33, 0.0, 0.67, 1.0]); ax2.axis('off')
ax2.text(0.0, 0.5, 'corona3d', fontsize=96, family='serif', color=NAVY, va='center', fontweight='bold')
fig.savefig('docs/logo_wordmark.png', dpi=100, transparent=True); plt.close(fig)
print('ok')
