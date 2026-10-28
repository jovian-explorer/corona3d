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
IMPROVE_MIN = -0.01                                    # a refit without a view is kept unless the mean IoU of the kept views falls by more than this
TIE_IOU = 0.02                                         # two fits whose mean IoU differs by less than this are tied (the free fit wins)
DEV_AREA = 3.0                                         # R_sun^2: an envelope this large counts as a developed CME (t_ref choice)
DEV_MARGIN = 0.2                                       # R_sun: ... when its front is this far inside the outer edge of the field
SCAN_MIN = 150.0                                       # the automatic t_ref is searched up to this many minutes after t_out
MIN_SEP = 30.0                                         # deg: views whose observers are closer than this cannot fix the direction
CONE = 15.0                                            # deg: default half-width of the direction cone around the flare site
HEIGHT_MAX = 15.0                                      # R_sun: absolute upper bound of the apex height of a fit
HEIGHT_FOV_FACTOR = 1.5                                # the apex height of a fit is at most this factor times the largest outer field-of-view radius of the fitted views
ALPHA_MIN = 5.0                                        # deg: a half-angle below this is a degenerate (needle) shell and the fit counts as poor
BASE_MIN = 10.0                                        # base frame: the last frame at or before onset minus this many minutes


def pick_t_ref(meas, r_out, t_out, area=DEV_AREA, margin=DEV_MARGIN, scan_min=SCAN_MIN, tol_min=6.0):
    """Reference time of the six-parameter fit chosen from the data: the earliest frame time, at or after
    t_out and within scan_min minutes of it, at which the CME is well developed in the largest number
    of views. A view counts as developed when the frame of that view within tol_min minutes has an
    envelope of at least `area` R_sun^2 whose front lies at least `margin` inside the outer edge of
    the field. meas is {view: [(time, envelope area, front distance in R_sun)]}, r_out {view: outer edge}.
    Returns (time, {view: True/False}) or (None, {}) when no view is ever developed."""
    t_out = datetime.fromisoformat(t_out) if isinstance(t_out, str) else t_out
    cand = sorted({t for rows in meas.values() for t, _, _ in rows
                   if timedelta(0) <= t - t_out <= timedelta(minutes=scan_min)})
    best, best_n = (None, {}), 0
    for t in cand:
        dev = {}
        for k, rows in meas.items():
            near = [r for r in rows if abs((r[0] - t).total_seconds()) < tol_min * 60]
            _, a, f = min(near, key=lambda r: abs((r[0] - t).total_seconds())) if near else (None, 0.0, 1e9)
            dev[k] = a >= area and f <= r_out[k] - margin
        n = sum(dev.values())
        if n > best_n:
            best, best_n = (t, dev), n
    return best


def poor_view(iou, iou_min=IOU_MIN, ratio=POOR_RATIO):
    """Name of the view of {view: IoU} that is poor, or None. A view is poor when its IoU is below
    iou_min and below ratio times the mean IoU of the other views; of several, the lowest. A single
    view is never poor (nothing to compare with)."""
    bad = {k: v for k, v in iou.items() if len(iou) > 1 and v < iou_min
           and v < ratio * np.mean([w for j, w in iou.items() if j != k])}
    return min(bad, key=bad.get) if bad else None


def _tie_pick(cands, tie=TIE_IOU):
    """Index of the preferred refit among candidate dicts (keys bounds, mean_iou, views): a fit without
    parameters at a bound beats one with them; then the higher mean IoU; refits within `tie` of the best
    mean IoU go to the one that keeps more views, then to the views nearer the Sun (smaller outer edge
    of the field of view)."""
    clean = [i for i, c in enumerate(cands) if not c["bounds"]]
    pool = clean or list(range(len(cands)))
    top = max(cands[i]["mean_iou"] for i in pool)
    near = [i for i in pool if top - cands[i]["mean_iou"] < tie]
    sun = lambda i: float(np.mean([FOV[k][1] if k in FOV else 0.0 for k in cands[i]["views"]]))
    return min(near, key=lambda i: (-len(cands[i]["views"]), sun(i), -cands[i]["mean_iou"]))


def check_views(views, q, iou, solve, score, iou_min=IOU_MIN, ratio=POOR_RATIO, improve=IMPROVE_MIN, drop=True,
                limit=None, problems=None):
    """Refit without the poor views. `views` is {name: view}, `q` the parameters of the fit to all of
    them, `iou` {name: IoU} of that fit; solve(views) returns the parameters of a fit to the given
    views and score(views, q) their {name: IoU}.

    Rule (documented in docs/METHODS.md): while more than one view remains and poor_view() names a
    view, the fit is repeated without it, and the refit replaces the previous fit unless the mean
    IoU over the remaining views falls by more than -`improve` compared with the previous fit on the
    same views (a view that the shell cannot match, such as an IoU of 0, costs the others nothing and
    only lowers the mean); otherwise the previous fit is kept and the search ends. The number of views never
    falls below one.

    When `limit(views)` is not None (the observers lie within MIN_SEP of each other, so the direction is
    tied to the flare site) the joint fit must not decide which view survives: every subset obtained by
    dropping one view (each single view, for two views) is fitted separately and the refits are compared
    by _tie_pick, using `problems(q, views)`, the list of parameters of a fit that sit at a bound. The
    chosen refit then faces the same improvement test. Returns (views, q, iou, record); the record
    lists the full fit and every refit, each with the key `bounds`."""
    rec = dict(rule=f"drop the view with IoU below {iou_min:g} and below {ratio:g} of the mean of the others, "
                    f"keep the refit unless the mean IoU of the kept views falls by more than {-improve:g}; at least one view stays",
               enabled=bool(drop), full=dict(views=list(views), iou={k: round(float(v), 3) for k, v in iou.items()}, mean_iou=round(float(np.mean(list(iou.values()))), 3),
                                             bounds=list(problems(q, views)) if problems else []),
               refits=[], dropped=[], kept="full")
    while drop and len(views) > 1:
        bad = poor_view(iou, iou_min, ratio)
        if bad is None:
            break
        close = limit is not None and limit(views) is not None
        subsets = ([[k for k in views if k != d] for d in views] if close else [[k for k in views if k != bad]])
        cands = []
        for names in subsets:
            rest = {k: views[k] for k in names}
            q2 = solve(rest)
            iou2 = score(rest, q2)
            after = float(np.mean(list(iou2.values())))
            before = float(np.mean([iou[k] for k in names]))
            cands.append(dict(dropped=[k for k in views if k not in names][0], views=list(names), iou={k: round(float(v), 3) for k, v in iou2.items()},
                              mean_iou=round(after, 3), mean_iou_before=round(before, 3), gain=after - before,
                              bounds=list(problems(q2, rest)) if problems else [], q=q2, rest=rest, kept=False))
        best = _tie_pick(cands) if close else 0
        ok = cands[best]["gain"] >= improve
        cands[best]["kept"] = bool(ok)
        for i, c in enumerate(cands):
            rec["refits"].append(dict({k: v for k, v in c.items() if k not in ("q", "rest", "gain")}, candidate=close, chosen=(i == best)))
        if not ok:
            break
        c = cands[best]
        views, q, iou = c["rest"], c["q"], c["iou"]
        rec["dropped"] += [k for k in rec["full"]["views"] if k not in views and k not in rec["dropped"]]
        rec["kept"] = "refit"
    return views, q, iou, rec


def view_lines(rec):
    """Text lines describing a view check record (for event_report.md and the command line)."""
    if not rec:
        return []
    f = rec["full"]
    hit = lambda b: f"; at a bound: {', '.join(b)}" if b else ""
    L = [f"Fit to all views ({', '.join(f['views'])}): IoU " + ", ".join(f"{k} {v:g}" for k, v in f["iou"].items())
         + f", mean {f['mean_iou']:g}{hit(f.get('bounds'))}."]
    cand = [r for r in rec["refits"] if r.get("candidate")]
    for r in rec["refits"]:
        end = ("chosen." if r["chosen"] else "not chosen.") if r.get("candidate") else               ("kept." if r["kept"] else "not kept (the mean IoU of the kept views fell).")
        L.append(f"Refit without {r['dropped']} ({', '.join(r['views'])}): IoU " + ", ".join(f"{k} {v:g}" for k, v in r["iou"].items())
                 + f", mean {r['mean_iou']:g} against {r['mean_iou_before']:g} for the same views in the fit to all{hit(r.get('bounds'))}; " + end)
    if cand:
        ch = next(r for r in cand if r.get("chosen"))
        clean = [r for r in cand if not r["bounds"]]
        why = ("only this refit has no parameter at a bound" if len(clean) == 1 and not ch["bounds"] else
               "all refits have a parameter at a bound, so the higher mean IoU decides" if not clean else
               "the mean IoU decides among the refits without a parameter at a bound; a tie keeps more views and then the view nearer the Sun")
        L.append(f"The observers lie close together, so each view was fitted on its own and the choice was made among these refits, not by the fit to all views: "
                 f"{', '.join(ch['views'])} was chosen ({why}) and the refit was " + ("kept." if ch["kept"] else "not kept because the mean IoU of the kept views fell."))
    if not rec["refits"]:
        L.append("No view was poor." if rec["enabled"] else "Views were kept as given (--keep-all-views).")
    return L


def observer_separation(xyz):
    """Largest angle in degrees between any two observers, given their Cartesian positions (0 for one view)."""
    v = [np.asarray(x, float) / np.linalg.norm(x) for x in xyz]
    return max([float(np.degrees(np.arccos(np.clip(a @ b, -1, 1)))) for i, a in enumerate(v) for b in v[i + 1:]] or [0.0])


def direction_limit(xyz, cone=CONE, min_sep=MIN_SEP):
    """Whether the direction must be tied to the source region (Thernisien et al. 2011): true for one view or
    when the observers are closer than `min_sep` deg. `xyz` maps each view to its observer position. Returns None
    (free direction) or dict(cone_deg, separation_deg, views, why, text); cone 0 means the direction is fixed at
    the flare site."""
    sep = observer_separation(list(xyz.values()))
    if len(xyz) > 1 and sep >= min_sep:
        return None
    why = "single viewpoint" if len(xyz) == 1 else f"viewpoints within {sep:.0f} deg, below {min_sep:g} deg"
    return dict(cone_deg=float(cone), separation_deg=round(sep, 1), views=list(xyz), why=why,
                text=f"direction constrained by the source region ({why})")


def height_limit(views):
    """Upper bound of the apex height in R_sun for a fit to `views` (names): HEIGHT_FOV_FACTOR times the largest
    outer field-of-view radius among them, never above HEIGHT_MAX. The images cannot constrain an apex far
    outside the field, where only the legs of the shell would be matched."""
    outer = [FOV[k][1] for k in views if k in FOV]
    return min(HEIGHT_MAX, HEIGHT_FOV_FACTOR * max(outer)) if outer else HEIGHT_MAX


def bounds(flare, fix_direction=False, cone=None, height_max=HEIGHT_MAX):
    """Bounds of the optimised parameters: height, kappa, alpha, [lon, lat unless fixed], tilt. `cone` (deg)
    narrows the direction to the flare site plus or minus the cone; `height_max` is the upper height bound
    (height_limit of the fitted views)."""
    dl, db = (40.0, 35.0) if not cone else (min(cone, 40.0), min(cone, 35.0))
    return [(1.5, float(height_max)), (0.1, 0.8), (0.0, 60.0)] + ([] if fix_direction else
            [(flare[0] - dl, flare[0] + dl), (flare[1] - db, flare[1] + db)]) + [(-90.0, 90.0)]


def bound_hits(t, flare, fix_direction=False, height_max=HEIGHT_MAX):
    """Parameters of `t` (dict of KEYS, or the six values in KEYS order) within 0.5 % of their range of a bound
    of the optimisation: [(name, value, low, high)]. lon/lat are not tested when the direction is fixed or tied;
    `height_max` is the upper height bound used in the fit."""
    t = dict(zip(KEYS, t)) if not isinstance(t, dict) else t
    names = ["height", "kappa", "alpha_deg"] + ([] if fix_direction else ["lon_stonyhurst", "lat_deg"]) + ["tilt_pythea_deg"]
    return [(n, t[n], lo, hi) for n, (lo, hi) in zip(names, bounds(flare, fix_direction, None, height_max))
            if min(t[n] - lo, hi - t[n]) <= 0.005 * (hi - lo)]


def quality(out, flare, fix_direction=False, iou_min=IOU_MIN):
    """Problems of a result of fit(): [] when it is acceptable. A fit is poor when the mean IoU
    over the instruments is below iou_min, when no instrument gave an IoU, or when a fitted
    parameter lies within 0.5 % of its range from a bound of the optimisation (the upper height bound is the one
    recorded in out["height_bound"], else HEIGHT_MAX). It is also poor when the half-angle is below ALPHA_MIN (a
    degenerate shell) or when the apex height at t_ref lies outside the outer field of view of every fitted view."""
    t, why = out["t_ref"], []
    tied = bool((out.get("direction") or {}).get("constrained"))   # a prior from the source region, not a limit of the data
    fix_direction = fix_direction or tied
    iou = t.get("iou") or {}
    if not iou:
        why.append("no instrument had an envelope")
    elif np.mean(list(iou.values())) < iou_min:
        why.append(f"mean IoU {np.mean(list(iou.values())):.2f} is below {iou_min:g}")
    hmax = (out.get("height_bound") or {}).get("max", HEIGHT_MAX)
    for n, v, lo, hi in bound_hits(t, flare, fix_direction, hmax):
        why.append(f"{n} = {v:g} is at its bound ({lo:g} to {hi:g})")
    if t["alpha_deg"] < ALPHA_MIN:
        why.append(f"half-angle alpha = {t['alpha_deg']:g} deg is below {ALPHA_MIN:g} deg (a degenerate shell)")
    outer = [FOV[k][1] for k in (out.get("views_used") or iou) if k in FOV]
    if outer and t["height"] > max(outer):
        why.append(f"the apex at {t['height']:g} R_sun lies outside the field of view of every fitted view (out to {max(outer):g} R_sun)")
    return why


def mean_iou(out):
    """Mean IoU over the views of a result of fit(), or None when no view had an envelope."""
    iou = (out.get("t_ref") or {}).get("iou") or {}
    return float(np.mean(list(iou.values()))) if iou else None


def choose_fit(free, fixed, tie=TIE_IOU):
    """Better of the free and the fixed-direction fit, each a dict(status, mean_iou) with status "good",
    "poor" or "failed" (or None when the fit was not made). A good fit beats a poor one; between two fits of
    the same status the higher mean IoU of the kept views wins, and a difference below `tie` goes to the free
    fit, because the fixed fit is a special case of it. Returns "free", "fixed direction" or None."""
    ok = lambda f: f is not None and f.get("status") in ("good", "poor")
    if not ok(free) and not ok(fixed):
        return None
    if not ok(fixed):
        return "free"
    if not ok(free):
        return "fixed direction"
    if free["status"] != fixed["status"]:
        return "free" if free["status"] == "good" else "fixed direction"
    return "fixed direction" if (fixed["mean_iou"] or 0.0) - (free["mean_iou"] or 0.0) >= tie else "free"


def _load():
    """Import the packages needed by a fit (not by the pure functions below)."""
    global u, sunpy, SkyCoord, fits, gcs, binary_closing, binary_dilation, binary_fill_holes
    global gaussian_filter, label, differential_evolution, convex_hull_image, frames
    import astropy.units as u
    import sunpy.map
    from astropy.coordinates import SkyCoord
    from astropy.io import fits
    from ._extras import need
    gcs = need("gcs", "PyThea.geometrical_models", "the GCS fit").gcs
    from scipy.ndimage import binary_closing, binary_dilation, binary_fill_holes, gaussian_filter, label
    from scipy.optimize import differential_evolution
    convex_hull_image = need("gcs", "skimage.morphology", "the GCS fit").convex_hull_image
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


def base_cutoff(onset, minutes=BASE_MIN):
    """Time of the fixed pre-event base window: `onset` (the flare onset, or the reconstruction time when
    no onset is given) minus `minutes`. A datetime or ISO text is accepted."""
    onset = datetime.fromisoformat(onset) if isinstance(onset, str) else onset
    return onset - timedelta(minutes=minutes)


def bases_of(fr, cutoff=None):
    """{(instrument, programme): total-brightness map} of the frame subtracted in the envelopes. With
    `cutoff` it is the last frame of the programme at or before that time, so the base does not depend
    on when the download started; when the programme has no frame that early, the earliest frame is
    used (listed in the second returned value, {(instrument, programme): first frame time}). Without
    `cutoff` the earliest frame of each programme is used."""
    out, late = {}, {}
    for k, frs in fr.items():
        progs = {}
        for t, files in frs:
            progs.setdefault(seq(files[0]), []).append((t, files))
        for s, lst in progs.items():
            early = [f for f in lst if cutoff is None or f[0] <= cutoff]
            if not early:
                late[(k, s)] = lst[0][0].isoformat()
            out[(k, s)] = total_b((early[-1] if early and cutoff is not None else lst[0])[1])
    return out, late


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
        if img.data.shape != b.data.shape:                # a frame with another binning cannot be subtracted
            continue
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



def auto_t_ref(fr, base, keys, flare, t_out, running=None):
    """Choose t_ref from the images (see pick_t_ref): the envelope area and front distance of every
    frame of each view in the 150 minutes after t_out are measured, and the earliest time at which the
    CME is developed in most views is taken. When no view is developed, t_out plus 60 minutes is used and
    the record says so. Returns (ISO time, record for gcs_fit)."""
    t0 = datetime.fromisoformat(t_out) if isinstance(t_out, str) else t_out
    meas = {}
    for k in keys:
        meas[k] = []
        for tf, files in fr[k]:
            if not t0 <= tf <= t0 + timedelta(minutes=SCAN_MIN):
                continue
            b = base_of(fr, base, k, tf, files, running)
            if b is None or b.date.isot[:19] >= tf.isoformat():
                continue
            img = total_b(files)
            if img.data.shape != b.data.shape:            # a frame with another binning cannot be subtracted
                continue
            dm = sunpy.map.Map(img.data - b.data, img.meta)
            env, _ = envelope(dm, k, flare)
            if not env.any():
                meas[k].append((tf, 0.0, 0.0))
                continue
            ppr = (dm.rsun_obs / dm.scale.axis1).to_value(u.pix)
            cx, cy = (float(v) for v in dm.wcs.world_to_pixel(SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame=dm.coordinate_frame)))
            yy, xx = np.nonzero(env)
            meas[k].append((tf, float(env.sum() / ppr**2), float(np.hypot(xx - cx, yy - cy).max() / ppr)))
    t, dev = pick_t_ref({k: v for k, v in meas.items() if v}, {k: FOV[k][1] for k in meas}, t0)
    rule = (f"earliest frame time within {SCAN_MIN:g} min after t_out at which the envelope covers at least "
            f"{DEV_AREA:g} R_sun^2 with its front {DEV_MARGIN:g} R_sun inside the field edge in most views")
    if t is None:
        t = t0 + timedelta(minutes=60)
        print(f"t_ref: no developed CME in {keys}; t_out + 60 min = {t:%Y-%m-%dT%H:%M:%S} used")
        return t.strftime("%Y-%m-%dT%H:%M:%S"), dict(method="fallback", rule=rule, why="no developed envelope found")
    print(f"t_ref chosen: {t:%Y-%m-%dT%H:%M:%S} (developed in {', '.join(k for k, v in dev.items() if v)})")
    return t.strftime("%Y-%m-%dT%H:%M:%S"), dict(method="auto", rule=rule, developed=[k for k, v in dev.items() if v])


def fit(wl, t_ref, t_out, flare=FLARE, keys=VIEWS, fix_direction=False, front=None, previous=None,
        running=None, height_out=None, height_note=None, drop_poor_views=True, iou_min=IOU_MIN,
        onset=None, base_minutes=BASE_MIN, direction_cone=CONE):
    """Two-step fit of a CME from the images in the folder wl; returns the result dictionary
    (keys t_ref, t_out, height_time, fit). `t_ref` None: chosen from the images (auto_t_ref). `previous` is an earlier result whose step 1 is reused,
    `running` the minutes of a running difference (default: difference from the first frame), and
    `height_out` an apex height at t_out given by hand (with `height_note`) in place of step 2.
    With `drop_poor_views` a view whose IoU is far below the others is left out and the fit repeated
    (see check_views); the result then has the key `view_check` with both fits. The base frame is the
    last frame of each programme at or before `onset` (default t_out) minus `base_minutes`, recorded in
    the key `base` of the result. With one view, or views whose observers are closer than MIN_SEP degrees,
    the direction is not free: it is fixed at the flare site (`direction_cone` 0) or bounded within that many
    degrees of it (key `direction` of the result)."""
    _load()
    flare = tuple(flare)
    fr = {k: frames_of(wl, k) for k in keys}
    onset = onset or t_out
    cutoff = base_cutoff(onset, base_minutes)
    base, late = bases_of(fr, cutoff)                    # fixed pre-event window, not the first downloaded frame
    base_rec = dict(onset=str(onset), minutes=base_minutes, cutoff=cutoff.isoformat(),
                    frames={f"{k}/{s}": m.date.isot[:19] for (k, s), m in base.items()},
                    fallback_earliest=[f"{k}/{s}" for (k, s) in late])
    print("bases (last frame at or before %s):" % cutoff.isoformat(), base_rec["frames"])
    if late:
        print("  no frame that early for", base_rec["fallback_earliest"], "- earliest frame used; download an earlier start")

    choice = None
    if t_ref is None:                                    # the frame where the CME is well developed
        t_ref, choice = auto_t_ref(fr, base, keys, flare, t_out, running)
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
        obs = {k: next(m for (kk, _), m in base.items() if kk == k).observer_coordinate
                   .transform_to(frames.HeliographicStonyhurst).cartesian.xyz.value for k in keys}

        def limit(vs):
            return None if fix_direction else direction_limit({k: obs[k] for k in vs}, direction_cone)

        def solve(vs):
            lim = limit(vs)
            fixed = fix_direction or (lim is not None and lim["cone_deg"] <= 0)
            full = (lambda q: [q[0], q[1], q[2], flare[0], flare[1], q[3]]) if fixed else list
            cone = lim["cone_deg"] if lim else None
            cost = lambda q: 1.0 if q[1] >= 0.95 else 1.0 - np.mean([v.iou(shell(full(q), t_ref)) for _, v in vs.values()])
            return full(list(differential_evolution(cost, bounds(flare, fixed, cone, height_limit(vs)), popsize=12, maxiter=60, tol=1e-4,
                                                    seed=3, polish=True).x))

        score = lambda vs, q: {k: float(v.iou(shell(q, t_ref))) for k, (_, v) in vs.items()}
        problems = lambda q, vs: [h[0] for h in bound_hits(q, flare, fix_direction or limit(vs) is not None, height_limit(vs))]
        q24 = solve(V)
        V, q24, _, rec = check_views(V, q24, score(V, q24), solve, score, iou_min, drop=drop_poor_views,
                                     limit=limit, problems=problems)
        for line in view_lines(rec):
            print("views:", line)
    keys = tuple(V)
    lim = None if previous else limit(V)
    out = {"base": base_rec, "views_used": list(keys),
           "height_bound": dict(max=height_limit(keys), factor=HEIGHT_FOV_FACTOR, cap=HEIGHT_MAX, views=list(keys),
                                outer_fov=max(FOV[k][1] for k in keys if k in FOV) if any(k in FOV for k in keys) else None),
           "t_ref": dict(zip(KEYS, [round(float(v), 3) for v in q24]), time=t_ref,
                         iou={k: round(float(v.iou(shell(q24, t_ref))), 3) for k, (_, v) in V.items()})}
    if rec:
        out["view_check"] = rec
    if lim:
        out["direction"] = dict(lim, constrained=True, flare=list(flare))
        print("direction:", lim["text"], f"(cone {lim['cone_deg']:g} deg around the flare site)")
    if choice:
        out["t_ref_choice"] = choice
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
        if b is None or b.date.isot[:19] >= tf.isoformat():  # the base frame itself, or a frame before it
            continue
        img = total_b(files)
        if img.data.shape != b.data.shape:                    # a frame with another binning cannot be subtracted
            continue
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
    rec = {"t_ref": out["t_ref"]["time"], "t_out": out["t_out"]["time"], "views": list(out.get("views_used") or [k for k in keys if k not in dropped and k in (out["t_ref"].get("iou") or {k: 1})]),
           "iou": out["t_ref"]["iou"], "flare": list(flare), "running_diff_min": running,
           "fixed": ["lon_stonyhurst", "lat_deg"] if fix_direction or (out.get("direction") or {}).get("cone_deg", 1) <= 0 else [],
           "height_t_out": out["fit"]["height_t_out"]}
    if out.get("base"):
        rec["base"] = out["base"]
    if out.get("height_bound"):
        rec["height_bound"] = out["height_bound"]
    if out.get("direction"):
        rec["direction_constraint"] = {k: v for k, v in out["direction"].items() if k != "constrained"}
    if out.get("t_ref_choice"):
        rec["t_ref_choice"] = out["t_ref_choice"]
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
    t_ref = a.t_ref or (T_REF if same else None)             # None: chosen from the images
    t_out = a.t_out or (T_OUT if same else cfg["time"])
    out_file = Path(home) / "gcs_fit.json"
    previous = json.loads(out_file.read_text(encoding="utf-8")) if a.skip_step1 else None
    out = fit(a.wl or wl_dir(home), t_ref, t_out, flare, keys, a.fix_direction, a.front, previous,
              a.running_diff, a.height_out, a.height_note,
              drop_poor_views=getattr(a, "drop_poor_views", True), onset=getattr(a, "onset", None),
              base_minutes=getattr(a, "base_minutes", BASE_MIN),
              direction_cone=getattr(a, "direction_cone", CONE))
    write_result(home, out, flare, keys, a.fix_direction, a.running_diff)
    if a.t_ref:
        set_wl_time(home, a.t_ref)
    elif out.get("t_ref_choice"):
        set_wl_time(home, out["t_ref"]["time"], "t_ref chosen automatically")


def write_result(home, out, flare, keys, fix_direction=False, running=None):
    """Write gcs_fit.json and the keys gcs, gcs_wl, gcs_fit of config.json for a result of fit()."""
    out_file = Path(home) / "gcs_fit.json"
    out_file.write_text(json.dumps(out, indent=1), encoding="utf-8")
    f = Path(home) / "config.json"
    user = json.loads(f.read_text(encoding="utf-8"))
    user.update(config_update(out, flare, keys, fix_direction, running))
    f.write_text(json.dumps(user, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {out_file} and the keys gcs, gcs_wl, and gcs_fit of {f}")


def set_wl_time(home, t_ref, source="--t-ref"):
    """Set `wl_time` of config.json (to --t-ref, the automatically chosen t_ref, or the nearest time with
    coronagraph frames), so that the white-light figure uses frames of the time of the fit (a coronagraph
    can have a gap, as LASCO C2 on 2024-05-09 from 14:48 to 18:52 UT), and record this in the `gcs_fit`
    provenance (`wl_time`, `wl_time_from`)."""
    f = Path(home) / "config.json"
    user = json.loads(f.read_text(encoding="utf-8"))
    old = json.dumps(user, sort_keys=True)
    user["wl_time"] = t_ref
    if "gcs_fit" in user:
        user["gcs_fit"] = dict(user["gcs_fit"], wl_time=t_ref, wl_time_from=source)
    if json.dumps(user, sort_keys=True) == old:          # unchanged: keep the file (and the figures newer than it)
        return
    f.write_text(json.dumps(user, indent=2) + "\n", encoding="utf-8")
    print(f"wl_time of {f} set to {t_ref} ({source})")
