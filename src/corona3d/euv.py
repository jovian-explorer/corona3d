"""Extreme-ultraviolet (EUV) images of the low corona and the apex height of the shell measured in them.

`corona3d fetch --euv` downloads SDO/AIA 193 and 211 A and STEREO/EUVI 195 A images (of the STEREO
spacecraft of the event) through the Virtual Solar Observatory, at one image every `euv_sample_min`
minutes within `euv_window`, into data/euv.

`corona3d gcs --euv` (and `corona3d event --euv`) then measure the apex height of the shell at the early
times of the burst in these images, where the coronagraphs do not yet see the CME. The shape of the shell
fitted in white light (aspect ratio, half-angle, direction, tilt) is kept, and at every EUV frame within
`EUV_SPAN_MIN` minutes of t_out the apex height is chosen on a grid of 0.01 R_sun from 1.05 to 3 R_sun
so that the mean intersection over union of the projected shell with the EUV envelopes of the views is
largest. The envelope of an EUV frame is formed as for the coronagraphs (base difference, robust
threshold, connected to the flare site, convex hull), in the annulus above the limb from 1.02 R_sun to
1.3 R_sun (AIA) or 1.7 R_sun (EUVI), so that the on-disc part of the EUV wave is left out. With two or
more frames, a straight line through the EUV heights gives the height at t_out with the standard error
of its intercept and the early speed of the apex; with one frame, the height of that frame is used, with
the half-width of the heights whose mean IoU lies within 0.02 of the best as its uncertainty.
"""
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

EUV_SPAN_MIN = 10.0                       # EUV frames within this many minutes of t_out are used
H_GRID = np.round(np.arange(1.05, 3.0001, 0.01), 3)
IOU_TOL = 0.02                            # heights within this IoU of the best give the uncertainty (one frame)
FOV = {"aia193": (1.02, 1.3, 1.1), "aia211": (1.02, 1.3, 1.1),
       "euvi195a": (1.02, 1.7, 1.15), "euvi195b": (1.02, 1.7, 1.15)}     # r_in, r_out, r_seed [R_sun]


DEFAULT_VIEWS = ("aia193", "euvi195a", "euvi195b")    # one EUV view per viewpoint; aia211 on request


def euv_dir(h):
    return Path(h) / "data" / "euv"


def fetch(h, cfg):
    """Download AIA 193 and 211 A and EUVI 195 A of the STEREO spacecraft of the event within euv_window
    (default: 20 min before to 15 min after the event time), one image every euv_sample_min minutes
    (default 5); returns the number of files that failed."""
    import astropy.units as u
    from sunpy.net import Fido, attrs as a
    from .whitelight import utc
    t = utc(cfg["time"])
    win = cfg.get("euv_window") or [(t - timedelta(minutes=20)).isoformat(), (t + timedelta(minutes=15)).isoformat()]
    d = euv_dir(h)
    d.mkdir(parents=True, exist_ok=True)
    tr, smp = a.Time(*win), a.Sample(float(cfg.get("euv_sample_min", 5.0)) * u.min)
    queries = [(a.Instrument("AIA"), a.Wavelength(193 * u.AA), smp), (a.Instrument("AIA"), a.Wavelength(211 * u.AA), smp),
               (a.Instrument("SECCHI"), a.Detector("EUVI"), a.Wavelength(195 * u.AA), a.Source(f"STEREO_{cfg['stereo_sc']}"), smp)]
    failed = 0
    for q in queries:
        try:
            r = Fido.search(tr, *q)
        except Exception as e:                       # an archive that does not answer: the others still run
            print(f"EUV search failed ({e}); run corona3d fetch --euv again")
            failed += 1
            continue
        print(r)
        res = Fido.fetch(r, path=str(d / "{file}"), max_conn=2, progress=False)
        for e in res.errors:
            print(f"download failed: {e}")
        failed += len(res.errors)
    return failed


def key_of(header):
    """aia193, aia211, euvi195a, or euvi195b for a FITS header; None for other images."""
    wl = int(round(float(header.get("WAVELNTH", 0) or 0)))
    tel, det = str(header.get("TELESCOP", "")).upper(), str(header.get("DETECTOR", "")).upper()
    if "AIA" in tel or "AIA" in str(header.get("INSTRUME", "")).upper():
        return f"aia{wl}" if wl in (193, 211) else None
    if det == "EUVI" and wl == 195:
        obs = str(header.get("OBSRVTRY", "")).upper()
        return "euvi195" + ("a" if obs.endswith("A") else "b")
    return None


def frames(folder):
    """{key: [(time, [file])]} of the EUV images in a folder, earliest first."""
    from astropy.io import fits
    out = {}
    for f in sorted(list(Path(folder).glob("*.fts")) + list(Path(folder).glob("*.fits"))):
        try:
            hd = fits.getheader(f, 1) if f.name.startswith("aia") else fits.getheader(f)
        except Exception:
            hd = fits.getheader(f)
        if "WAVELNTH" not in hd:
            hd = fits.getheader(f)
        k = key_of(hd)
        if k:
            out.setdefault(k, []).append((datetime.fromisoformat(str(hd["DATE-OBS"]).replace("/", "-")[:19]), [f]))
    return {k: sorted(v, key=lambda x: x[0]) for k, v in out.items()}


def best_height(score, grid=H_GRID, tol=IOU_TOL):
    """(height, mean IoU, (low, high)) of the grid height with the largest score(h), and the span of the
    heights whose score lies within tol of it."""
    s = np.array([score(h) for h in grid])
    i = int(np.argmax(s))
    near = grid[s >= s[i] - tol]
    return float(grid[i]), float(s[i]), (float(near.min()), float(near.max()))


def height_profile(folder, q, t_out, flare, keys=None, onset=None, base_minutes=10.0, span=EUV_SPAN_MIN, speed_km_s=None):
    """Apex heights of the shell with the shape q (kappa, alpha, lon, lat, tilt in the PyThea convention)
    in the EUV frames within `span` minutes of t_out, and the height at t_out. Returns the record written
    to `euv` of gcs_fit.json, or None when no EUV envelope was found."""
    from . import gcs
    gcs._load()
    gcs.FOV.update(FOV)
    t_out = datetime.fromisoformat(t_out) if isinstance(t_out, str) else t_out
    fr = frames(folder)
    keys = [k for k in (keys or DEFAULT_VIEWS) if k in fr]
    if not keys:
        return None
    fr = {k: fr[k] for k in keys}
    base, late = gcs.bases_of(fr, gcs.base_cutoff(onset or t_out, base_minutes))
    # the frame times of the first view (AIA 193 when present); the other views take their frame within 90 s
    times = [t for t, _ in fr[keys[0]] if abs((t - t_out).total_seconds()) <= span * 60]
    rows = []
    for t in times:
        V = gcs.views_at(fr, base, t, keys, flare)
        V = {k: v for k, v in V.items() if v[0] == t or abs((v[0] - t).total_seconds()) < 90}
        if not V:
            continue
        score = lambda h: float(np.mean([v.iou(gcs.shell([h, *q], t.isoformat())) for _, v in V.values()]))
        h, iou, (lo, hi) = best_height(score)
        r_out = max(FOV[k][1] for k in V)
        why = ("at the end of the height grid" if h <= H_GRID[0] or h >= H_GRID[-1] else
               f"apex beyond the outer edge of the field ({r_out:g} R_sun), so that the envelope constrains only the flanks"
               if h > r_out else "")
        rows.append(dict(time=t.isoformat(), views=sorted(V), height=h, iou=round(iou, 3), low=lo, high=hi,
                         used=not why, rejected=why or None))
        print(f"EUV {t:%H:%M:%S} ({', '.join(sorted(V))}): apex height {h:.2f} R_sun, IoU {iou:.3f} ({lo:.2f} to {hi:.2f})")
    for r in rows:
        if r["rejected"]:
            print(f"EUV {r['time'][11:19]} not used: {r['rejected']}")
    used = [r for r in rows if r["used"]]
    if not used:
        return None
    rec = dict(frames=rows, keys=keys, span_min=span, grid=[float(H_GRID[0]), float(H_GRID[-1]), 0.01],
               base={f"{k}/{s}": m.date.isot[:19] for (k, s), m in base.items()}, fallback_earliest=[f"{k}/{s}" for (k, s) in late])
    if len(used) >= 3:                     # a straight line through the EUV heights
        v, h0, err = gcs.extrapolate([datetime.fromisoformat(r["time"]) for r in used], [r["height"] for r in used], t_out)
        rec.update(height_t_out=round(h0, 3), height_t_out_err=round(err, 3), speed_km_s=round(v * 695700 / 60, 0), method="line")
        return rec
    # one or two frames: the frame nearest to t_out, moved to t_out with the white-light speed when known
    r = min(used, key=lambda r: abs((datetime.fromisoformat(r["time"]) - t_out).total_seconds()))
    dt = (t_out - datetime.fromisoformat(r["time"])).total_seconds()
    shift = speed_km_s * dt / 695700 if speed_km_s else 0.0
    rec.update(height_t_out=round(r["height"] + shift, 3), height_t_out_err=round((r["high"] - r["low"]) / 2, 3),
               speed_km_s=None, method=f"frame at {r['time']}" + (f", moved by {dt:+.0f} s at the white-light speed of "
                                                                   f"{speed_km_s:.0f} km/s" if speed_km_s and dt else ""))
    if len(used) == 2:
        t2 = [datetime.fromisoformat(x["time"]) for x in used]
        rec["speed_km_s_two_frames"] = round((used[1]["height"] - used[0]["height"]) * 695700 / (t2[1] - t2[0]).total_seconds(), 0)
    return rec
