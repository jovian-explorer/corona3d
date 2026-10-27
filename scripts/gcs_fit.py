"""GCS fit of the CME of 2013-10-25 from white-light coronagraph images, in two steps.

1. At the time of the reference frames (C2 15:24, COR1-A and COR1-B 15:25 UT) the six GCS
   parameters are fitted by differential evolution to the CME envelopes in the three images.
2. The shape (aspect ratio, half-angle, direction, tilt) is kept, and the apex height at the time
   of the reconstruction (15:00 UT) is extrapolated from the height-time profile of the CME front,
   as Morosan et al. (2022) scaled their 15:24 UT shell back in time. The front is measured in the
   instrument whose line of sight is most nearly perpendicular to the CME direction (the limb
   view, COR1-B here), as the largest distance of the envelope from Sun centre, deprojected with
   the fitted direction, in every frame where the envelope covers at least 1 R_sun^2 and lies
   0.2 R_sun or more inside the outer edge of the field of view, up to the first frame in which
   the front reaches that edge; a straight line through these heights gives the height at 15:00 UT.

Envelope of one frame: exposure-normalised total brightness minus a pre-event frame, smoothed
with a Gaussian of 2 px, pixels brighter than 3 robust standard deviations above the median,
the regions connected to the radial line above the flare site (S06E69) at r_seed, filled to
their convex hull, and only inside the field of view of the instrument. A shell is compared
with an envelope through the intersection over union (IoU) of its projected silhouette.

usage: python scripts/gcs_fit.py WL_DIR [--out gcs_fit.json] [--skip-step1]
WL_DIR holds the LASCO C2 (clear polariser, 1024 px), COR1-A and COR1-B (s4c1A/B) images,
e.g. data/wl of a corona3d working directory after `corona3d fetch --wl --cor1a`.
"""
import argparse
import json
from datetime import datetime
from pathlib import Path

import astropy.units as u
import numpy as np
import sunpy.map
from astropy.coordinates import SkyCoord
from astropy.io import fits
from PyThea.geometrical_models import gcs
from scipy.ndimage import binary_closing, binary_dilation, binary_fill_holes, gaussian_filter, label
from scipy.optimize import differential_evolution
from skimage.morphology import convex_hull_image
from sunpy.coordinates import frames

FLARE = (-69.0, -6.0)                                  # Stonyhurst lon, lat of the flare site
FOV = {"c2": (2.3, 6.0, 3.0), "cor1a": (1.5, 4.0, 2.2), "cor1b": (1.5, 4.0, 2.2)}   # r_in, r_out, r_seed
T_REF = "2013-10-25T15:24:30"
T_OUT = "2013-10-25T15:00:00"
KEYS = ["height", "kappa", "alpha_deg", "lon_stonyhurst", "lat_deg", "tilt_pythea_deg"]


def _t(f):
    h = fits.getheader(f)
    if h.get("TIME-OBS"):
        return datetime.strptime(h["DATE-OBS"] + " " + h["TIME-OBS"][:8], "%Y/%m/%d %H:%M:%S")
    return datetime.fromisoformat(h["DATE-OBS"][:19])


def frames_of(wl, key):
    """[(time, [files])] of total-brightness frames of one instrument, earliest first."""
    d = Path(wl)
    if key == "c2":
        fs = [f for f in sorted(d.glob("2*.fts")) if f.stem.isdigit()
              and fits.getheader(f).get("NAXIS1") == 1024 and fits.getheader(f).get("POLAR") == "Clear"]
        return [(_t(f), [f]) for f in fs]
    fs = sorted(d.glob(f"*_s[45]c1{key[-1].upper()}.fts"), key=_t)
    out, cur = [], []
    for f in fs:                                       # triplets: three exposures within 60 s
        if cur and (_t(f) - _t(cur[0])).total_seconds() > 60:
            out.append(cur); cur = []
        cur.append(f)
    out.append(cur)
    return [(_t(g[0]), g) for g in out if len(g) == 3]


def total_b(files):
    ms = [sunpy.map.Map(f) for f in files]
    d = sum(m.data.astype(float) / m.exposure_time.to_value(u.s) for m in ms)
    m = sunpy.map.Map(d * (2 / 3 if len(ms) == 3 else 1), ms[0].meta)
    return m.superpixel([2, 2] * u.pix, func=np.mean) if m.data.shape[0] > 512 else m


def envelope(dm, key, k=3.0):
    r_in, r_out, r_seed = FOV[key]
    ppr = (dm.rsun_obs / dm.scale.axis1).to_value(u.pix)
    cx, cy = (float(v) for v in dm.wcs.world_to_pixel(SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame=dm.coordinate_frame)))
    yy, xx = np.mgrid[:dm.data.shape[0], :dm.data.shape[1]]
    rr = np.hypot(xx - cx, yy - cy) / ppr
    usable = np.isfinite(dm.data) & (rr > r_in) & (rr < r_out)
    d = gaussian_filter(np.nan_to_num(dm.data), 2.0)
    med = np.median(d[usable])
    mad = np.median(np.abs(d[usable] - med)) * 1.4826
    sig = (d - med > k * mad) & usable
    s = SkyCoord(FLARE[0] * u.deg, FLARE[1] * u.deg, r_seed * u.R_sun,
                 frame=frames.HeliographicStonyhurst, obstime=dm.date).transform_to(dm.coordinate_frame)
    fx, fy = (float(v) for v in dm.wcs.world_to_pixel(s))
    seed = np.hypot(xx - fx, yy - fy) < 0.5 * ppr
    lab, _ = label(sig)
    keep = np.isin(lab, np.unique(lab[seed & sig])) & (lab > 0)
    env = binary_fill_holes(convex_hull_image(keep)) & usable if keep.any() else keep
    return env, usable


class View:
    def __init__(self, dm, key):
        self.map, self.key = dm, key
        self.mask, self.usable = envelope(dm, key)
        self.obs = dm.observer_coordinate.transform_to(frames.HeliographicStonyhurst).cartesian.xyz.to_value(u.R_sun)

    def silhouette(self, P):
        sc = SkyCoord(*(P.T * u.R_sun), frame=frames.HeliographicStonyhurst, obstime=self.map.date,
                      representation_type="cartesian").transform_to(self.map.coordinate_frame)
        px, py = (np.asarray(v) for v in self.map.wcs.world_to_pixel(sc))
        uo = self.obs / np.linalg.norm(self.obs)
        rpar = P @ uo
        rperp = np.linalg.norm(P - np.outer(rpar, uo), axis=1)
        vis = (np.linalg.norm(P, axis=1) > 1) & ~((rperp < 1) & (rpar < 0))
        px, py = px[vis], py[vis]
        m = np.zeros(self.mask.shape, bool)
        ok = (px >= 0) & (px < m.shape[1]) & (py >= 0) & (py < m.shape[0])
        m[py[ok].astype(int), px[ok].astype(int)] = True
        return binary_fill_holes(binary_closing(binary_dilation(m, iterations=2), iterations=4))

    def iou(self, P):
        s, e = self.silhouette(P) & self.usable, self.mask & self.usable
        return (s & e).sum() / max((s | e).sum(), 1)


def shell(q, obstime):
    h, kappa, alpha, lon, lat, tilt = q
    rc = gcs.rcenter_(h * u.R_sun, alpha * u.deg, kappa)
    c = SkyCoord(lon * u.deg, lat * u.deg, rc, frame=frames.HeliographicStonyhurst, obstime=obstime)
    m = gcs(c, h * u.R_sun, alpha * u.deg, kappa, tilt * u.deg, nbvertsl=10, nbvertcirc=20, nbvertcircshell=60)
    X, Y, Z = m.rotate(*(v * u.R_sun for v in m.cloud()))
    return np.c_[X.to_value(u.R_sun), Y.to_value(u.R_sun), Z.to_value(u.R_sun)]


def seq(f):
    """Observing programme of a frame (s4/s5 for COR1, one programme for C2): each frame is
    compared only with a pre-event frame of its own programme."""
    return f.name.split("_")[-1][:2] if "_" in f.name else "c2"


def bases_of(fr):
    out = {}
    for k, frs in fr.items():
        for _, files in frs:
            out.setdefault((k, seq(files[0])), total_b(files))
    return out


def views_at(fr, base, t, keys):
    """One View per instrument for the frame nearest to t (within 6 min)."""
    vs = {}
    for k in keys:
        cand = [f for f in fr[k] if abs((f[0] - t).total_seconds()) < 360]
        if not cand:
            continue
        f = min(cand, key=lambda f: abs((f[0] - t).total_seconds()))
        img = total_b(f[1])
        dm = sunpy.map.Map(img.data - base[(k, seq(f[1][0]))].data, img.meta)
        v = View(dm, k)
        if v.mask.sum() > 20:                            # an envelope was found
            vs[k] = (f[0], v)
    return vs


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("wl")
    ap.add_argument("--out", default="gcs_fit.json")
    ap.add_argument("--skip-step1", action="store_true", help="take step 1 from an existing --out file")
    a = ap.parse_args()
    keys = ("c2", "cor1a", "cor1b")
    fr = {k: frames_of(a.wl, k) for k in keys}
    base = bases_of(fr)                                  # earliest frame of each programme
    print("bases:", {f"{k}/{s}": m.date.isot[:19] for (k, s), m in base.items()})

    t_ref = datetime.fromisoformat(T_REF)
    if a.skip_step1:
        q24 = [json.loads(Path(a.out).read_text())["t1524"][k] for k in KEYS]
    else:
        V = views_at(fr, base, t_ref, keys)
        cost = lambda q: 1.0 if q[1] >= 0.95 else 1.0 - np.mean([v.iou(shell(q, T_REF)) for _, v in V.values()])
        bounds = [(2.0, 8.0), (0.1, 0.8), (0.0, 60.0), (-110.0, -40.0), (-40.0, 30.0), (-90.0, 90.0)]
        q24 = list(differential_evolution(cost, bounds, popsize=12, maxiter=60, tol=1e-4, seed=3, polish=True).x)
    V = views_at(fr, base, t_ref, keys)
    out = {"t1524": dict(zip(KEYS, [round(float(v), 3) for v in q24]),
                         iou={k: round(float(v.iou(shell(q24, T_REF))), 3) for k, (_, v) in V.items()})}
    print("step 1:", out["t1524"])

    # step 2: leading edge of the CME, in the instrument whose line of sight is most nearly
    # perpendicular to the fitted CME direction (the limb view), deprojected with that direction
    L, B = np.radians(q24[3]), np.radians(q24[4])
    d = np.array([np.cos(B) * np.cos(L), np.cos(B) * np.sin(L), np.sin(B)])
    proj = {}
    for k in keys:
        o = next(m for (kk, _), m in base.items() if kk == k).observer_coordinate.transform_to(frames.HeliographicStonyhurst).cartesian.xyz.value
        proj[k] = float(np.sqrt(1 - (d @ (o / np.linalg.norm(o))) ** 2))
    k = max(proj, key=proj.get)
    print(f"limb view: {k} (projection factor {proj[k]:.2f})")
    hts = []
    for tf, files in fr[k]:
        b = base[(k, seq(files[0]))]
        if b.date.isot[:19] == tf.isoformat():             # the base frame itself
            continue
        img = total_b(files)
        dm = sunpy.map.Map(img.data - b.data, img.meta)
        env, _ = envelope(dm, k)
        ppr = (dm.rsun_obs / dm.scale.axis1).to_value(u.pix)
        area = env.sum() / ppr**2
        if area < 1.0:                                   # no CME, or noise at the seed
            continue
        cx, cy = (float(v) for v in dm.wcs.world_to_pixel(SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame=dm.coordinate_frame)))
        yy, xx = np.nonzero(env)
        r_sky = float(np.hypot(xx - cx, yy - cy).max() / ppr)
        if r_sky > FOV[k][1] - 0.2:                      # the front has reached the edge of the
            break                                        # field: later fronts are truncated
        hts.append((tf.isoformat(), k, round(r_sky / proj[k], 3), round(float(area), 2)))
        print("  ", hts[-1])
    t0 = datetime.fromisoformat(T_OUT)
    x = np.array([(datetime.fromisoformat(t) - t0).total_seconds() / 60 for t, *_ in hts])
    y = np.array([h for _, _, h, _ in hts])
    (v, h0), cov = np.polyfit(x, y, 1, cov=True)
    out["height_time"] = dict(instrument=k, projection=round(proj[k], 3),
                              frames=[dict(time=t, height=h, area_rsun2=a) for t, _, h, a in hts])
    out["fit"] = dict(speed_rsun_per_min=round(float(v), 4), speed_km_s=round(float(v) * 695700 / 60, 0),
                      height_1500=round(float(h0), 3), height_1500_err=round(float(np.sqrt(cov[1, 1])), 3),
                      n_frames=len(hts))
    out["t1500"] = dict(out["t1524"], height=round(float(h0), 3))
    out["t1500"].pop("iou")
    print("step 2:", out["fit"])
    Path(a.out).write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
