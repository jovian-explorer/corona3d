"""LASCO C2 and STEREO-B COR1/COR2 base-difference images as planes in the three-dimensional reconstruction.

Each image is rotated to solar north up, base-differenced, binned, and placed on the plane
through Sun centre perpendicular to its observer's line of sight (x west, y north, R_sun), so
seen from that observer it overlays the model exactly.
"""
from datetime import datetime
from pathlib import Path

import astropy.units as u
import numpy as np
import pyvista as pv
import sunpy.map
from astropy.coordinates import SkyCoord
from astropy.io import fits

from .geometry import observer_basis

# inner and outer edge of the field of view [R_sun] (COR2 starts where COR1 ends, so the two
# coplanar STEREO-B images do not overlap); observer key in the configuration
INSTR = dict(c2=dict(r=(2.2, 6.0), obs="earth"),
             cor1=dict(r=(1.5, 4.0), obs="stereo_b"),
             cor2=dict(r=(4.0, 15.0), obs="stereo_b"))
NBIN = 256


def _time(f):
    h = fits.getheader(f)
    if "TIME-OBS" in h and h["TIME-OBS"]:
        return datetime.strptime(h["DATE-OBS"] + " " + h["TIME-OBS"][:8], "%Y/%m/%d %H:%M:%S")
    return datetime.fromisoformat(h["DATE-OBS"][:19])


def _frames(wl_dir, inst):
    """[(time, [files])] of total-brightness frames, earliest first."""
    d = Path(wl_dir)
    if inst == "c2":
        fs = [f for f in sorted(d.glob("2*.fts")) if f.stem.isdigit()]
        fs = [f for f in fs if fits.getval(f, "NAXIS1") == 1024
              and fits.getval(f, "POLAR") == "Clear" and fits.getval(f, "DETECTOR") == "C2"]
        return [(_time(f), [f]) for f in fs]
    if inst == "cor2":
        return [(_time(f), [f]) for f in sorted(d.glob("*_d4c2B.fts"))]
    groups = {}                                           # COR1: polarisation triplets
    for f in sorted(d.glob("*_s4c1B.fts")):
        groups.setdefault(f.stem[:13], []).append(f)
    return [(_time(g[0]), g) for g in groups.values() if len(g) == 3]


def _image(files):
    """Exposure-normalised, binned total-brightness map (COR1: tB = 2/3 sum of the three angles)."""
    maps = [sunpy.map.Map(f) for f in files]
    data = sum(m.data.astype(float) / m.exposure_time.to_value(u.s) for m in maps)
    if len(maps) == 3:
        data *= 2.0 / 3.0
    k = data.shape[0] // NBIN                             # bin to NBIN^2 (bases differ in size)
    return sunpy.map.Map(data, maps[0].meta).superpixel([k, k] * u.pix, func=np.mean)


def load_plane(inst, wl_dir, when, observer):
    frames = _frames(wl_dir, inst)
    if not frames:
        raise FileNotFoundError(f"no {inst} files in {wl_dir}")
    t = datetime.fromisoformat(when)
    base = _image(frames[0][1])
    img = _image(min(frames, key=lambda fr: abs(fr[0] - t))[1])
    diff = sunpy.map.Map(img.data - base.data, img.meta)
    diff = diff.rotate(order=1, missing=np.nan)
    c = diff.wcs.world_to_pixel(SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame=diff.coordinate_frame))
    s = diff.scale.axis1.to_value(u.arcsec / u.pix) / diff.rsun_obs.to_value(u.arcsec)
    ny, nx = diff.data.shape
    X, Y = np.meshgrid((np.arange(nx) - c[0]) * s, (np.arange(ny) - c[1]) * s)
    v = diff.data.copy()
    r = np.hypot(X, Y)
    lo, hi = INSTR[inst]["r"]
    v[(r < lo) | (r > hi)] = np.nan
    ex, ey, _ = observer_basis(*observer)
    P = X[..., None] * ex + Y[..., None] * ey
    grid = pv.StructuredGrid(P[..., 0, None], P[..., 1, None], P[..., 2, None])
    grid["wl"] = v[..., None].ravel(order="F")
    print(f"{inst}: {img.date.isot[:19]} - base {base.date.isot[:19]}")
    return grid
