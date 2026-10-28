"""GCS fit of a CME from white-light coronagraph images, in two steps (`corona3d gcs`).

1. At the time of the reference frames (--t-ref) the six GCS parameters are fitted by differential
   evolution to the CME envelopes in the images of the chosen instruments.
2. The shape (aspect ratio, half-angle, direction, tilt) is kept, and the apex height at the time
   of the reconstruction (--t-out) is extrapolated from the height-time profile of the CME front,
   as Morosan et al. (2022) scaled their 15:24 UT shell back in time. The front is measured in the
   instrument whose line of sight is most nearly perpendicular to the CME direction (the limb
   view), as the largest distance of the envelope from Sun centre, deprojected with the fitted
   direction, in every frame where the envelope covers at least 2 R_sun^2 and lies 0.2 R_sun or
   more inside the outer edge of the field of view, up to the first frame in which the front
   reaches that edge; a straight line through these heights gives the height at --t-out.

Envelope of one frame: exposure-normalised total brightness minus a pre-event frame, smoothed
with a Gaussian of 2 px, pixels brighter than 3 robust standard deviations above the median,
the regions connected to the radial line above the flare site at r_seed, filled to their convex
hull, and only inside the field of view of the instrument. A shell is compared with an envelope
through the intersection over union (IoU) of its projected silhouette.

The pure functions (projection_factor, front_profile, extrapolate) need only NumPy; PyThea,
SunPy, and scikit-image are imported when a fit is run.
"""
import json
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

# defaults: the event of 2013-10-25; other events set them on the command line
FLARE = (-69.0, -6.0)                                  # Stonyhurst lon, lat of the flare site (seed of the envelopes)
FOV = {"c2": (2.3, 6.0, 3.0), "cor1a": (1.5, 4.0, 2.2), "cor1b": (1.5, 4.0, 2.2),
       "cor2a": (3.0, 15.0, 4.0), "cor2b": (3.0, 15.0, 4.0)}                      # r_in, r_out, r_seed
T_REF = "2013-10-25T15:24:30"
T_OUT = "2013-10-25T15:00:00"
VIEWS = ("c2", "cor1a", "cor1b")
KEYS = ["height", "kappa", "alpha_deg", "lon_stonyhurst", "lat_deg", "tilt_pythea_deg"]

IOU_MIN = 0.5                                          # mean IoU below which a fit counts as poor (`corona3d event`)
POOR_RATIO = 0.75                                      # a view is poor below this fraction of the mean IoU of the others
IMPROVE_MIN = 0.05                                     # a refit without a view must raise the mean IoU of the kept views by this


def poor_view(iou, iou_min=IOU_MIN, ratio=POOR_RATIO):
    """Name of the view of {view: IoU} that is poor, or None. A view is poor when its IoU is below
    iou_min and below ratio times the mean IoU of the other views; of several, the lowest. A single
    view is never poor (nothing to compare with)."""
    bad = {k: v for k, v in iou.items() if len(iou) > 1 and v < iou_min
           and v < ratio * np.mean([w for j, w in iou.items() if j != k])}
    return min(bad, key=bad.get) if bad else None


def check_views(views, q, iou, solve, score, iou_min=IOU_MIN, ratio=POOR_RATIO, improve=IMPROVE_MIN, drop=True):
    """Refit without the poor views. `views` is {name: view}, `q` the parameters of the fit to all of
    them, `iou` {name: IoU} of that fit; solve(views) returns the parameters of a fit to the given
    views and score(views, q) their {name: IoU}.

    Rule (documented in docs/METHODS.md): while more than one view remains and poor_view() names a
    view, the fit is repeated without it, and the refit replaces the previous fit only when the mean
    IoU over the remaining views rises by at least `improve` compared with the previous fit on the
    same views; otherwise the previous fit is kept and the search ends. The number of views never
    falls below one. Returns (views, q, iou, record); the record lists the full fit and every refit."""
    rec = dict(rule=f"drop the view with IoU below {iou_min:g} and below {ratio:g} of the mean of the others, "
                    f"keep the refit if the mean IoU of the kept views rises by {improve:g} or more; at least one view stays",
               enabled=bool(drop), full=dict(views=list(views), iou=dict(iou), mean_iou=round(float(np.mean(list(iou.values()))), 3)),
               refits=[], dropped=[], kept="full")
    while drop and len(views) > 1:
        bad = poor_view(iou, iou_min, ratio)
        if bad is None:
            break
        rest = {k: v for k, v in views.items() if k != bad}
        q2 = solve(rest)
        iou2 = score(rest, q2)
        before = float(np.mean([iou[k] for k in rest]))
        after = float(np.mean(list(iou2.values())))
        ok = after - before >= improve
        rec["refits"].append(dict(dropped=bad, views=list(rest), iou=dict(iou2), mean_iou=round(after, 3),
                                  mean_iou_before=round(before, 3), kept=bool(ok)))
        if not ok:
            break
        views, q, iou = rest, q2, iou2
        rec["dropped"].append(bad)
        rec["kept"] = "refit"
    return views, q, iou, rec


def view_lines(rec):
    """Text lines describing a view check record (for event_report.md and the command line)."""
    if not rec:
        return []
    f = rec["full"]
    L = [f"Fit to all views ({', '.join(f['views'])}): IoU " + ", ".join(f"{k} {v:g}" for k, v in f["iou"].items())
         + f", mean {f['mean_iou']:g}."]
    for r in rec["refits"]:
        L.append(f"Refit without {r['dropped']} ({', '.join(r['views'])}): IoU " + ", ".join(f"{k} {v:g}" for k, v in r["iou"].items())
                 + f", mean {r['mean_iou']:g} against {r['mean_iou_before']:g} for the same views in the fit to all; "
                 + ("kept." if r["kept"] else "not kept (the rise is below the required gain)."))
    if not rec["refits"]:
        L.append("No view was poor." if rec["enabled"] else "Views were kept as given (--keep-all-views).")
    return L


def bounds(flare, fix_direction=False):
    """Bounds of the optimised parameters: height, kappa, alpha, [lon, lat unless fixed], tilt."""
    return [(1.5, 15.0), (0.1, 0.8), (0.0, 60.0)] + ([] if fix_direction else
            [(flare[0] - 40, flare[0] + 40), (flare[1] - 35, flare[1] + 35)]) + [(-90.0, 90.0)]


def quality(out, flare, fix_direction=False, iou_min=IOU_MIN):
    """Problems of a result of fit(): [] when it is acceptable. A fit is poor when the mean IoU
    over the instruments is below iou_min, when no instrument gave an IoU, or when a fitted
    parameter lies within 0.5 % of its range from a bound of the optimisation."""
    t, why = out["t_ref"], []
    iou = t.get("iou") or {}
    if not iou:
        why.append("no instrument had an envelope")
    elif np.mean(list(iou.values())) < iou_min:
        why.append(f"mean IoU {np.mean(list(iou.values())):.2f} is below {iou_min:g}")
    names = ["height", "kappa", "alpha_deg"] + ([] if fix_direction else ["lon_stonyhurst", "lat_deg"]) + ["tilt_pythea_deg"]
    for n, (lo, hi) in zip(names, bounds(flare, fix_direction)):
        if min(t[n] - lo, hi - t[n]) <= 0.005 * (hi - lo):
            why.append(f"{n} = {t[n]:g} is at its bound ({lo:g} to {hi:g})")
    return why


def _load():
    """Import the packages needed by a fit (not by the pure functions below)."""
    global u, sunpy, SkyCoord, fits, gcs, binary_closing, binary_dilation, binary_fill_holes
    global gaussian_filter, label, differential_evolution, convex_hull_image, frames
    import astropy.units as u
    import sunpy.map
    from astropy.coordinates import SkyCoord
    from astropy.io import fits
    from PyThea.geometrical_models import gcs
    from scipy.ndimage import binary_closing, binary_dilation, binary_fill_holes, gaussian_filter, label
    from scipy.optimize import differential_evolution
    from skimage.morphology import convex_hull_image
    from sunpy.coordinates import frames



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
    if key.startswith("cor2"):
        return [(_t(f), [f]) for f in sorted(d.glob(f"*_d4c2{key[-1].upper()}.fts"), key=_t)]
    fs = sorted(d.glob(f"*_[sn][45]c1{key[-1].upper()}.fts"), key=_t)
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
    k = m.data.shape[0] // 512                         # bin to 512 pixels across
    return m.superpixel([k, k] * u.pix, func=np.mean) if k > 1 else m


def envelope(dm, key, flare=FLARE, k=3.0):
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
    s = SkyCoord(flare[0] * u.deg, flare[1] * u.deg, r_seed * u.R_sun,
                 frame=frames.HeliographicStonyhurst, obstime=dm.date).transform_to(dm.coordinate_frame)
    fx, fy = (float(v) for v in dm.wcs.world_to_pixel(s))
    rf = np.hypot(fx - cx, fy - cy) / ppr
    if rf < r_in + 0.3:                                # seed hidden by the occulter (a CME seen
        k = (r_in + 0.5) / max(rf, 1e-3)               # near the disc): move it outward along
        fx, fy = cx + (fx - cx) * k, cy + (fy - cy) * k  # the position angle of the flare
    seed = np.hypot(xx - fx, yy - fy) < 0.5 * ppr
    lab, _ = label(sig)
    keep = np.isin(lab, np.unique(lab[seed & sig])) & (lab > 0)
    env = binary_fill_holes(convex_hull_image(keep)) & usable if keep.any() else keep
    return env, usable


class View:
    def __init__(self, dm, key, flare=FLARE):
        self.map, self.key = dm, key
        self.mask, self.usable = envelope(dm, key, flare)
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


def base_of(fr, base, k, t, files, running=None):
    """Frame subtracted from the frame at t: the earliest frame of its programme, or, with
    `running` (minutes), the frame of the same programme nearest to t minus `running`
    (a running difference); None when there is no earlier frame."""
    if not running:
        return base[(k, seq(files[0]))]
    cand = [f for f in fr[k] if seq(f[1][0]) == seq(files[0]) and f[0] < t]
    if not cand:
        return None
    t0 = t - timedelta(minutes=running)
    return total_b(min(cand, key=lambda f: abs((f[0] - t0).total_seconds()))[1])


def views_at(fr, base, t, keys, flare=FLARE, running=None):
    """One View per instrument for the frame nearest to t (within 6 min)."""
    t = datetime.fromisoformat(t) if isinstance(t, str) else t
    vs = {}
    for k in keys:
        cand = [f for f in fr[k] if abs((f[0] - t).total_seconds()) < 360]
        if not cand:
            continue
        f = min(cand, key=lambda f: abs((f[0] - t).total_seconds()))
        b = base_of(fr, base, k, f[0], f[1], running)
        if b is None:
            continue
        img = total_b(f[1])
        dm = sunpy.map.Map(img.data - b.data, img.meta)
        v = View(dm, k, flare)
        if v.mask.sum() > 20:                            # an envelope was found
            vs[k] = (f[0], v)
    return vs


def projection_factor(lon, lat, observer_xyz):
    """Fraction of a radial distance along (lon, lat) [deg, Stonyhurst] that the observer sees in
    the plane of the sky: sqrt(1 - cos^2 of the angle between the direction and the observer)."""
    L, B = np.radians(lon), np.radians(lat)
    d = np.array([np.cos(B) * np.cos(L), np.cos(B) * np.sin(L), np.sin(B)])
    o = np.asarray(observer_xyz, float)
    return float(np.sqrt(1 - (d @ (o / np.linalg.norm(o))) ** 2))


def front_profile(rows, r_out, min_area=2.0, margin=0.2):
    """rows: [(time, r_sky, area)] in time order. Keeps the uninterrupted run of frames with an
    envelope of at least min_area R_sun^2 that ends at the first frame whose front lies within
    margin of the outer edge r_out of the field (that and later fronts are truncated). A frame
    below min_area ends the run before it, so that small envelopes of earlier, unrelated
    structure (a streamer blob before the eruption) do not flatten the height-time line."""
    kept = []
    for t, r, a in rows:
        if a < min_area:
            kept = []
            continue
        if r > r_out - margin:
            break
        kept.append((t, r, a))
    return kept


def extrapolate(times, heights, t0):
    """Straight line through (time, height); returns speed [R_sun/min], height at t0, its error."""
    x = np.array([(t - t0).total_seconds() / 60 for t in times])
    (v, h0), cov = np.polyfit(x, np.asarray(heights, float), 1, cov=True)
    return float(v), float(h0), float(np.sqrt(cov[1, 1]))



def fit(wl, t_ref, t_out, flare=FLARE, keys=VIEWS, fix_direction=False, front=None, previous=None,
        running=None, height_out=None, height_note=None, drop_poor_views=True, iou_min=IOU_MIN):
    """Two-step fit of a CME from the images in the folder wl; returns the result dictionary
    (keys t_ref, t_out, height_time, fit). `previous` is an earlier result whose step 1 is reused,
    `running` the minutes of a running difference (default: difference from the first frame), and
    `height_out` an apex height at t_out given by hand (with `height_note`) in place of step 2.
    With `drop_poor_views` a view whose IoU is far below the others is left out and the fit repeated
    (see check_views); the result then has the key `view_check` with both fits."""
    _load()
    flare = tuple(flare)
    fr = {k: frames_of(wl, k) for k in keys}
    base = bases_of(fr)                                  # earliest frame of each programme
    print("bases:", {f"{k}/{s}": m.date.isot[:19] for (k, s), m in base.items()})

    tr = datetime.fromisoformat(t_ref)
    allkeys = tuple(keys)
    if previous:
        q24 = [previous["t_ref"][k] for k in KEYS]
        used = tuple(k for k in keys if k in (previous["t_ref"].get("iou") or {})) or tuple(keys)
        V = views_at(fr, base, tr, used, flare, running)
        rec = None
    else:
        V = views_at(fr, base, tr, keys, flare, running)
        if not V:
            raise SystemExit(f"no CME envelope in any of {keys} within 6 min of {t_ref}; check the times and --flare")
        full = (lambda q: [q[0], q[1], q[2], flare[0], flare[1], q[3]]) if fix_direction else list

        def solve(vs):
            cost = lambda q: 1.0 if q[1] >= 0.95 else 1.0 - np.mean([v.iou(shell(full(q), t_ref)) for _, v in vs.values()])
            return full(list(differential_evolution(cost, bounds(flare, fix_direction), popsize=12, maxiter=60, tol=1e-4,
                                                    seed=3, polish=True).x))

        score = lambda vs, q: {k: float(v.iou(shell(q, t_ref))) for k, (_, v) in vs.items()}
        q24 = solve(V)
        V, q24, _, rec = check_views(V, q24, score(V, q24), solve, score, iou_min, drop=drop_poor_views)
        for line in view_lines(rec):
            print("views:", line)
    keys = tuple(V)
    out = {"t_ref": dict(zip(KEYS, [round(float(v), 3) for v in q24]), time=t_ref,
                         iou={k: round(float(v.iou(shell(q24, t_ref))), 3) for k, (_, v) in V.items()})}
    if rec:
        out["view_check"] = rec
    print("step 1:", out["t_ref"])
    if height_out is not None:                           # step 2 replaced by an assumed height
        out["fit"] = dict(height_t_out=float(height_out), assumed=True, note=height_note or "")
        out["t_out"] = dict(out["t_ref"], height=float(height_out), time=t_out)
        out["t_out"].pop("iou")
        return out

    # step 2: leading edge of the CME, in the instrument whose line of sight is most nearly
    # perpendicular to the fitted CME direction (the limb view), deprojected with that direction
    proj = {}
    for k in allkeys:
        o = next(m for (kk, _), m in base.items() if kk == k).observer_coordinate
        proj[k] = projection_factor(q24[3], q24[4], o.transform_to(frames.HeliographicStonyhurst).cartesian.xyz.value)
    k = front or max(keys, key=proj.get)           # the limb view among the views kept in the fit
    print(f"front measured in {k} (projection factor {proj[k]:.2f})")
    rows = []
    for tf, files in fr[k]:
        b = base_of(fr, base, k, tf, files, running)
        if b is None or b.date.isot[:19] == tf.isoformat():  # the base frame itself
            continue
        img = total_b(files)
        dm = sunpy.map.Map(img.data - b.data, img.meta)
        env, _ = envelope(dm, k, flare)
        ppr = (dm.rsun_obs / dm.scale.axis1).to_value(u.pix)
        if not env.any():
            rows.append((tf, 0.0, 0.0))
            continue
        cx, cy = (float(v) for v in dm.wcs.world_to_pixel(SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame=dm.coordinate_frame)))
        yy, xx = np.nonzero(env)
        rows.append((tf, float(np.hypot(xx - cx, yy - cy).max() / ppr), float(env.sum() / ppr**2)))
    hts = [(t.isoformat(), k, round(r / proj[k], 3), round(a, 2)) for t, r, a in front_profile(rows, FOV[k][1])]
    for h in hts:
        print("  ", h)
    if len(hts) < 2:
        raise SystemExit(f"fewer than two usable fronts in {k}; the height at {t_out} cannot be extrapolated")
    v, h0, h0_err = extrapolate([datetime.fromisoformat(t) for t, *_ in hts], [h for _, _, h, _ in hts],
                                datetime.fromisoformat(t_out))
    out["height_time"] = dict(instrument=k, projection=round(proj[k], 3),
                              frames=[dict(time=t, height=h, area_rsun2=a) for t, _, h, a in hts])
    out["fit"] = dict(speed_rsun_per_min=round(v, 4), speed_km_s=round(v * 695700 / 60, 0),
                      height_t_out=round(h0, 3), height_t_out_err=round(h0_err, 3), n_frames=len(hts))
    out["t_out"] = dict(out["t_ref"], height=round(h0, 3), time=t_out)
    out["t_out"].pop("iou")
    print("step 2:", out["fit"])
    return out


def config_update(out, flare, keys, fix_direction=False, running=None):
    """The config.json keys `gcs` (at t_out), `gcs_wl` (at t_ref), and the provenance record
    `gcs_fit` for a result of fit(); the tilt of the configuration is the PyThea tilt minus 90 deg."""
    def shell_of(d, with_time):
        s = {k: d[k] for k in ("height", "kappa", "alpha_deg", "lon_stonyhurst", "lat_deg")}
        s["tilt_deg"] = round(d["tilt_pythea_deg"] - 90.0, 3)
        return dict(time=d["time"], **s) if with_time else s
    dropped = (out.get("view_check") or {}).get("dropped", [])
    rec = {"t_ref": out["t_ref"]["time"], "t_out": out["t_out"]["time"], "views": [k for k in keys if k not in dropped],
           "iou": out["t_ref"]["iou"], "flare": list(flare), "running_diff_min": running,
           "fixed": ["lon_stonyhurst", "lat_deg"] if fix_direction else [],
           "height_t_out": out["fit"]["height_t_out"]}
    if dropped:
        rec["views_dropped"] = list(dropped)
    if out["fit"].get("assumed"):
        rec.update(height_source="assumed", height_note=out["fit"]["note"])
    else:
        rec.update(height_source="measured", front=out["height_time"]["instrument"],
                   speed_km_s=out["fit"]["speed_km_s"], height_t_out_err=out["fit"]["height_t_out_err"],
                   n_frames=out["fit"]["n_frames"])
    return {"gcs": shell_of(out["t_out"], False), "gcs_wl": shell_of(out["t_ref"], True), "gcs_fit": rec}


def run(home, cfg, a):
    """The `corona3d gcs` command: fit, write gcs_fit.json and the keys of config.json."""
    from .paths import default_config, wl_dir
    if a.height_out is not None and not a.height_note:
        raise SystemExit("--height-out needs --height-note, the reason for the assumed height")
    keys = tuple(a.views.split(","))
    same = cfg["time"] == default_config()["time"]           # the event of 2013-10-25
    flare = tuple(a.flare or (cfg.get("gcs_fit") or {}).get("flare")
                  or (FLARE if same else (cfg["gcs"]["lon_stonyhurst"], cfg["gcs"]["lat_deg"])))
    t_ref = a.t_ref or (T_REF if same else cfg["gcs_wl"]["time"])
    t_out = a.t_out or (T_OUT if same else cfg["time"])
    out_file = Path(home) / "gcs_fit.json"
    previous = json.loads(out_file.read_text(encoding="utf-8")) if a.skip_step1 else None
    out = fit(a.wl or wl_dir(home), t_ref, t_out, flare, keys, a.fix_direction, a.front, previous,
              a.running_diff, a.height_out, a.height_note,
              drop_poor_views=getattr(a, "drop_poor_views", True))
    write_result(home, out, flare, keys, a.fix_direction, a.running_diff)
    if a.t_ref:
        set_wl_time(home, a.t_ref)


def write_result(home, out, flare, keys, fix_direction=False, running=None):
    """Write gcs_fit.json and the keys gcs, gcs_wl, gcs_fit of config.json for a result of fit()."""
    out_file = Path(home) / "gcs_fit.json"
    out_file.write_text(json.dumps(out, indent=1), encoding="utf-8")
    f = Path(home) / "config.json"
    user = json.loads(f.read_text(encoding="utf-8"))
    user.update(config_update(out, flare, keys, fix_direction, running))
    f.write_text(json.dumps(user, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {out_file} and the keys gcs, gcs_wl, and gcs_fit of {f}")


def set_wl_time(home, t_ref):
    """Set `wl_time` of config.json to the time given as --t-ref, so that the white-light figure uses
    frames of the time of the fit (a coronagraph can have a gap, as LASCO C2 on 2024-05-09 from 14:48
    to 18:52 UT), and record this in the `gcs_fit` provenance (`wl_time`, `wl_time_from`)."""
    f = Path(home) / "config.json"
    user = json.loads(f.read_text(encoding="utf-8"))
    user["wl_time"] = t_ref
    if "gcs_fit" in user:
        user["gcs_fit"] = dict(user["gcs_fit"], wl_time=t_ref, wl_time_from="--t-ref")
    f.write_text(json.dumps(user, indent=2) + "\n", encoding="utf-8")
    print(f"wl_time of {f} set to {t_ref} (--t-ref)")
