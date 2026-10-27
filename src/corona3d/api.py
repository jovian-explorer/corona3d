"""Analysis without the figure: every layer of the reconstruction as arrays.

``build`` reads the working directory's configuration and returns a :class:`Reconstruction`,
which holds the density surface, the field lines, the CME shell, the radio sources, and the
observers as NumPy arrays. The layers are those of ``corona3d render`` (they are built by the
same code), and ``corona3d.exporting`` writes them to files.

Coordinate frame and units of every position: Heliographic Carrington Cartesian, centred on the
Sun, X towards Carrington longitude 0, Z towards solar north, in solar radii (R_sun). The frame
rotates with the Sun; the observers are given by their Carrington longitude and latitude at the
time of the reconstruction. Plane-of-sky positions (``x_rsun``, ``y_rsun``) are seen from Earth,
x towards the west and y towards solar north.
"""
from dataclasses import dataclass, field

import numpy as np

from .layers import LAYERS, ON
from .paths import config as load_config

__all__ = ["Reconstruction", "build", "load_config", "carrington_angles", "FRAME", "UNITS"]

FRAME = "Heliographic Carrington Cartesian, Sun-centred, X to longitude 0, Z to solar north"
UNITS = {"position": "R_sun", "angle": "deg", "va": "km/s", "ne": "cm^-3", "tb": "K"}
RADIO_COLUMNS = ("instrument", "freq_mhz", "time_utc", "x_rsun", "y_rsun", "x_arcsec", "y_arcsec", "tb_k",
                 "r_rsun", "height_rsun", "x", "y", "z", "lon_carrington_deg", "lat_deg",
                 "lon_stonyhurst_deg")


@dataclass
class Reconstruction:
    """The layers of one reconstruction (``None`` or empty when a layer is switched off).

    Attributes
    ----------
    config : dict
        The configuration used (defaults, config.json, and the parameters passed to ``build``).
    layers : dict
        Layer name -> bool, as used.
    surface : dict or None
        ``vertices`` (n, 3) [R_sun, Carrington Cartesian], ``faces`` (m, 3) indices of triangles,
        ``va_kms`` (n,) Alfven speed [km/s] and ``ne_cm3`` (n,) electron density [cm^-3] per
        vertex (the iso-surface of ``ne_level``, rotated by ``model_dlon``).
    lines : dict
        ``open``, ``closed``, ``streamer``: lists of (k, 3) arrays [R_sun, Carrington Cartesian],
        drawn within ``line_r_max``. The streamer lines are the closed lines through the CME shell.
    gcs : dict
        Name (``time``, and ``image`` when ``gcs_wl`` is configured) -> dict with ``time`` (ISO)
        and ``rings`` (n_axis, n_circle, 3) [R_sun, Carrington Cartesian].
    radio : dict
        Columns ``RADIO_COLUMNS`` as equal-length arrays, one row per source position. Rows of
        ``radio_sources`` of the configuration have ``tb_k`` = NaN and an empty time.
    observers : dict
        ``earth`` and ``stereo`` -> ``lon_deg``, ``lat_deg`` (Carrington) and ``direction``
        (unit vector from Sun centre towards the observer, Carrington Cartesian).
    whitelight : dict
        Instrument -> ``points`` (n, 3) [R_sun], ``wl`` (n,), ``dims`` (3,), ``time``, ``base_time``.
    times : dict
        Times used by the layers, as returned by ``reconstruct``.
    """
    config: dict
    layers: dict
    surface: dict = None
    lines: dict = field(default_factory=lambda: dict(open=[], closed=[], streamer=[]))
    gcs: dict = field(default_factory=dict)
    radio: dict = field(default_factory=dict)
    observers: dict = field(default_factory=dict)
    whitelight: dict = field(default_factory=dict)
    times: dict = field(default_factory=dict)


def carrington_angles(xyz):
    """Carrington longitude and latitude of points.

    Parameters
    ----------
    xyz : array_like, shape (n, 3)
        Positions [R_sun], Carrington Cartesian.

    Returns
    -------
    lon, lat, r : ndarray, shape (n,)
        Longitude in [0, 360) deg, latitude in deg, distance from Sun centre in R_sun.
    """
    xyz = np.atleast_2d(xyz)
    r = np.linalg.norm(xyz, axis=1)
    return (np.degrees(np.arctan2(xyz[:, 1], xyz[:, 0])) % 360.0,
            np.degrees(np.arcsin(np.clip(xyz[:, 2] / r, -1, 1))), r)


def _observers(p):
    from .geometry import unit_vec
    return {k: dict(lon_deg=float(p[k][0]), lat_deg=float(p[k][1]), direction=unit_vec(*p[k]))
            for k in ("earth", "stereo")}


def _radio_table(sc, p):
    """Radio positions with their sky-plane and Carrington coordinates (radio_sources first)."""
    from astropy import units as u
    from sunpy.coordinates import sun
    rows = []
    rs = sun.angular_radius(p["time"]).to_value(u.arcsec)
    for s, (f, xyz, _) in zip(p["radio_sources"], sc.get("radio", [])[:len(p["radio_sources"])]):
        k = 1.0 / rs if s.get("unit") == "arcsec" else 1.0
        for x, y, q in zip(np.atleast_1d(s["x"]) * k, np.atleast_1d(s["y"]) * k, xyz):
            rows.append(("user", f, "", x, y, np.nan, q))
    for b in sc.get("radio_blobs", []):
        for t, x, y, tb, q in zip(b["time"], b["x"], b["y"], b["tb"], b["xyz"]):
            rows.append((b["instrument"], b["freq"], t.isoformat(), x, y, tb, q))
    cols = {c: [] for c in RADIO_COLUMNS}
    if rows:
        xyz = np.array([r[6] for r in rows])
        lon, lat, r = carrington_angles(xyz)
        for (inst, f, t, x, y, tb, _), q, lo, la, rr in zip(rows, xyz, lon, lat, r):
            vals = (inst, f, t, x, y, x * rs, y * rs, tb, rr, rr - 1.0, *q, lo, la,
                    (lo - p["earth"][0] + 180.0) % 360.0 - 180.0)
            for c, v in zip(RADIO_COLUMNS, vals):
                cols[c].append(v)
    return {c: np.array(v, dtype=str if c in ("instrument", "time_utc") else float) for c, v in cols.items()}


def build(home=".", params=None, **layers):
    """Build the reconstruction of a working directory and return its layers as arrays.

    Parameters
    ----------
    home : str or Path
        Working directory (its ``config.json`` is read; the data are in ``data/`` below it).
    params : dict, optional
        Configuration keys that replace those of ``config.json`` (nested keys may be partial).
    **layers : bool
        Any of ``LAYERS`` as True or False; the defaults are those of ``corona3d render``
        (the coronagraph layers and ``radio`` are off). ``sphere`` has no data and is ignored.

    Returns
    -------
    Reconstruction
        Positions in R_sun (Carrington Cartesian), Alfven speed in km/s, density in cm^-3,
        angles in degrees; see :class:`Reconstruction` for the shapes.

    Raises
    ------
    ValueError
        For a layer name that is not in ``LAYERS``.
    FileNotFoundError
        When a layer needs data that are missing (``corona3d fetch``).
    """
    from . import reconstruction as rc
    from .paths import merge
    bad = set(layers) - set(LAYERS)
    if bad:
        raise ValueError(f"unknown layers {sorted(bad)}; choose from {', '.join(LAYERS)}")
    on = dict(ON, **layers)
    p = merge(load_config(home), params)
    sc = rc.reconstruct(p, dict(on, gcs=False), home)
    rec = Reconstruction(config=p, layers=on, times=sc["times"], observers=_observers(p))
    if "surface" in sc:
        s = sc["surface"]
        rec.surface = dict(vertices=np.asarray(s.points, float), faces=s.faces.reshape(-1, 4)[:, 1:].copy(),
                           va_kms=np.asarray(s["va"], float), ne_cm3=np.asarray(s["ne"], float))
    for k in rec.lines:
        rec.lines[k] = sc.get(k, [])
    if on["gcs"]:
        rec.gcs["time"] = dict(time=p["time"], rings=rc.shell_of(p["gcs"], p, p["time"]))
        if p.get("gcs_wl"):
            t = p["gcs_wl"]["time"]
            rec.gcs["image"] = dict(time=t, rings=rc.shell_of(p["gcs_wl"], p, t))
    if on["radio"]:
        rec.radio = _radio_table(sc, p)
    for inst in rc.WL:
        if on[inst]:
            g = sc[inst]
            rec.whitelight[inst] = dict(points=np.asarray(g.points, float), wl=np.asarray(g["wl"], float),
                                        dims=np.array(g.dimensions), time=sc["times"][inst][0].isoformat(),
                                        base_time=sc["times"][inst][1].isoformat())
    return rec
