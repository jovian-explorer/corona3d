"""LASCO C2 and STEREO (A or B) COR1/COR2 images as planes in the three-dimensional reconstruction.

Each frame is normalised by its exposure time, binned, compared with the first frame of the
download window (a ratio, or a difference), rotated to solar north up, and placed on the plane
through Sun centre perpendicular to its observer's line of sight (x west, y north, R_sun), so
that seen from that observer it overlays the model at its true scale and orientation.
"""
from datetime import datetime
from pathlib import Path

import astropy.units as u
import numpy as np
import sunpy.map
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time

from .geometry import observer_basis

# inner and outer edge of the field of view [R_sun] (COR2 starts where COR1 ends, so the two
# coplanar STEREO images do not overlap); observer key in the configuration
INSTR = dict(c2=dict(r=(2.2, 6.0), obs="earth"),
             cor1=dict(r=(1.5, 4.0), obs="stereo"),
             cor2=dict(r=(4.0, 15.0), obs="stereo"))
NBIN = 256


def utc(t):
    """Naive UTC datetime from an ISO time with or without a time zone."""
    return Time(t, scale="utc").to_datetime()


def _time(f):
    h = fits.getheader(f)
    if h.get("TIME-OBS"):
        return datetime.strptime(h["DATE-OBS"] + " " + h["TIME-OBS"][:8], "%Y/%m/%d %H:%M:%S")
    return datetime.fromisoformat(h["DATE-OBS"][:19])


def _frames(wl_dir, inst, sc="B"):
    """[(time, [files])] of total-brightness frames, earliest first."""
    d = Path(wl_dir)
    if inst == "c2":
        fs = []
        for f in sorted(d.glob("2*.fts")):
            h = fits.getheader(f) if f.stem.isdigit() else {}
            if h.get("DETECTOR") == "C2" and h.get("POLAR") == "Clear" and h.get("NAXIS1") == 1024:
                fs.append(f)
        return [(_time(f), [f]) for f in fs]
    if inst == "cor2":
        return [(_time(f), [f]) for f in sorted(d.glob(f"*_d4c2{sc}.fts"))]
    out, cur = [], []                                  # COR1: three exposures within 60 s
    for f in sorted(d.glob(f"*_[sn]4c1{sc}.fts"), key=_time):
        if cur and (_time(f) - _time(cur[0])).total_seconds() > 60:
            out.append(cur)
            cur = []
        cur.append(f)
    out.append(cur)
    return [(_time(g[0]), g) for g in out if len(g) == 3]


def _image(files, shape=None):
    """Exposure-normalised total-brightness map (COR1: tB = 2/3 sum of the three angles), block
    averaged to about NBIN pixels across, and resampled to `shape` if it differs from it."""
    maps = [sunpy.map.Map(f) for f in files]
    data = sum(m.data.astype(float) / m.exposure_time.to_value(u.s) for m in maps)
    if len(maps) == 3:
        data *= 2.0 / 3.0
    k = max(1, data.shape[0] // NBIN)
    m = sunpy.map.Map(data, maps[0].meta).superpixel([k, k] * u.pix, func=np.mean)
    if shape is not None and m.data.shape != shape:
        m = m.resample(u.Quantity(shape[::-1], u.pix))
    return m


WL_MAX_GAP_MIN = 30.0      # a frame farther than this from wl_time is not used


def nearest_times(frames, t):
    """Text naming the frame times of `frames` ([(time, files)]) just before and just after t."""
    before = [f[0] for f in frames if f[0] <= t]
    after = [f[0] for f in frames if f[0] > t]
    parts = ([f"{max(before):%Y-%m-%dT%H:%M:%S} (before)"] if before else []) +             ([f"{min(after):%Y-%m-%dT%H:%M:%S} (after)"] if after else [])
    return "nearest available frames: " + (" and ".join(parts) if parts else "none in the folder")


def nearest_frame(wl_dir, inst, sc, window, when):
    """(frame time, signed gap in minutes) of the frame of `inst` inside `window` that lies nearest to
    `when` and after the base frame (the first frame of the window), or None when the window has no
    frame after the base frame."""
    t0, t1 = utc(window[0]), utc(window[1])
    t = utc(when)
    frames = [f for f in _frames(wl_dir, inst, sc) if t0 <= f[0] <= t1][1:]
    if not frames:
        return None
    tf = min(frames, key=lambda f: abs(f[0] - t))[0]
    return tf, (tf - t).total_seconds() / 60.0


def load_plane(inst, wl_dir, when, observer, window, mode="ratio", sc="B", nearest=False):
    """Plane of one instrument at the frame nearest to `when`, compared with the first frame of
    `window` ([start, end]); returns (grid, frame time, base time). With `nearest` a frame farther
    than WL_MAX_GAP_MIN is used too, as long as it is the nearest one after the base frame (the
    caller labels it); without it such a frame is refused."""
    t0, t1 = utc(window[0]), utc(window[1])
    every = _frames(wl_dir, inst, sc)
    frames = [f for f in every if t0 <= f[0] <= t1]
    t = utc(when)
    if not frames:
        raise FileNotFoundError(f"no {inst} frames between {window[0]} and {window[1]} in {wl_dir} "
                                f"({nearest_times(every, t)}); run `corona3d fetch --wl`")
    if not t0 <= t <= t1:
        raise ValueError(f"wl_time {when} lies outside wl_window {window}")
    if nearest:
        if len(frames) < 2:
            raise ValueError(f"no {inst} frame after the base frame in {window[0]} to {window[1]} ({nearest_times(every, t)})")
        tf, files = min(frames[1:], key=lambda fr: abs(fr[0] - t))
    else:
        tf, files = min(frames, key=lambda fr: abs(fr[0] - t))
    if not nearest and abs((tf - t).total_seconds()) > 60 * WL_MAX_GAP_MIN:
        raise ValueError(f"no {inst} frame within {WL_MAX_GAP_MIN:g} min of wl_time {when} ({nearest_times(every, t)}); "
                         "set wl_time (or --t-ref) to one of them")
    if tf == frames[0][0]:
        raise ValueError(f"the {inst} frame nearest to {when} is the base frame itself ({tf}); "
                         "choose a later wl_time")
    base = _image(frames[0][1])
    img = _image(files, base.data.shape)
    if mode == "ratio":                                # relative change; cancels vignetting
        with np.errstate(divide="ignore", invalid="ignore"):
            d = np.where(base.data > 0, img.data / base.data - 1.0, np.nan)
    elif mode == "difference":
        d = img.data - base.data
    else:
        raise ValueError('wl_mode must be "ratio" or "difference"')
    diff = sunpy.map.Map(d, img.meta).rotate(order=1, missing=np.nan)
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
    from ._extras import need
    pv = need("render", "pyvista", "the coronagraph layer")
    grid = pv.StructuredGrid(P[..., 0, None], P[..., 1, None], P[..., 2, None])
    grid["wl"] = v[..., None].ravel(order="F")
    tb = frames[0][0]
    print(f"{inst}: {img.date.isot[:19]} vs base {base.date.isot[:19]} ({mode})")
    return grid, tf, tb
