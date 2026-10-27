"""The three-dimensional reconstruction and the three-panel figure, every layer switchable.

MAS/MAST n_e iso-surface coloured by v_A, open and closed field lines, the closed lines threaded
by the CME shell (the streamer loops), a grey reference sphere, the GCS shell of the CME, radio
sources at their Newkirk heights, and LASCO C2 and STEREO-B COR1, COR2 images, seen from Earth,
the solar north pole, and STEREO-B (the layout of Fig. 7 of Morosan et al. 2022).
"""
from datetime import timedelta
from pathlib import Path

import astropy.units as u
import numpy as np
import pyvista as pv
from matplotlib import cm, colors
from matplotlib.figure import Figure
from matplotlib import colormaps
from scipy.spatial.transform import Rotation

from . import mast_model
from .geometry import CLOSED, OPEN, FieldTracer, deproject, fibonacci_sphere, observer_basis, pythea_shell, unit_vec
from .layers import LAYERS, ON
from .paths import config, default_config, merge, model_dir, wl_dir

WL = ("c2", "cor1", "cor2")
IMAGES_OF = {"earth": ("c2",), "stereo_b": ("cor2", "cor1")}


def newkirk_r(f_mhz, fold=4.0, harmonic=True):
    """Height [R_sun] of the plasma level of f_mhz in the fold-times Newkirk (1961) model."""
    fp = f_mhz / (2 if harmonic else 1) * 1e6
    ne = (fp / 8980.0) ** 2
    return 4.32 / np.log10(ne / (fold * 4.2e4))


def build_surface(m, p):
    PH, TH, R = np.meshgrid(m["phi"], m["theta"], m["r"], indexing="ij")
    grid = pv.StructuredGrid(R * np.sin(TH) * np.cos(PH), R * np.sin(TH) * np.sin(PH),
                             R * np.cos(TH))
    grid["ne"] = m["ne"].ravel(order="F")
    grid["va"] = m["va"].ravel(order="F")
    return grid.contour([p["ne_level"]], scalars="ne").clean()


def radio_points(p):
    """[(freq, xyz)] for the user-supplied sources, deprojected from the Earth view; arcsec are
    converted with the angular solar radius seen from Earth at the time of the reconstruction."""
    from sunpy.coordinates import sun
    rs = sun.angular_radius(p["time"]).to_value(u.arcsec)
    out = []
    for s in p["radio_sources"]:
        k = 1.0 / rs if s.get("unit") == "arcsec" else 1.0
        xyz = deproject(np.atleast_1d(s["x"]) * k, np.atleast_1d(s["y"]) * k, newkirk_r(s["freq"]), *p["earth"])
        out.append((s["freq"], xyz, np.tile(cm.cool(np.clip((s["freq"] - 150) / 300, 0, 1))[:3], (len(xyz), 1))))
    return out


def clip(lines, r_max):
    """Parts of each line inside r_max [R_sun]: lines are traced to the outer boundary of the
    model for their classification but drawn only near the Sun."""
    out = []
    for x in lines:
        inside = np.linalg.norm(x, axis=1) <= r_max
        cuts = np.flatnonzero(np.diff(inside.astype(int))) + 1
        out += [seg for seg, ok in zip(np.split(x, cuts), np.split(inside, cuts)) if ok[0] and len(seg) > 1]
    return out


def shell_of(g, p, time):
    return pythea_shell(g["height"], g["kappa"], g["alpha_deg"], g["lon_stonyhurst"], g["lat_deg"],
                        g["tilt_deg"] + 90.0, p["earth"][0], obstime=time,
                        nbvertsl=p["mesh_rings"]["nbvertsl"], nbvertcirc=p["mesh_rings"]["nbvertcirc"])


def reconstruct(p=None, on=None, data_root="."):
    """Build the layers that are switched on; p and on may hold only the keys that differ from
    the defaults. Returns a dictionary of the layers and of the times used."""
    p = merge(default_config(), p)
    on = dict(ON, **(on or {}))
    sc = {"times": {"shell": p["time"]}}
    if on["surface"] or on["open"] or on["closed"] or on["streamer"]:
        m = mast_model.load(model_dir(data_root, p), density_scale=p["density_scale"],
                            va_density=p["va_density"])
        R = Rotation.from_euler("z", p["model_dlon"], degrees=True).as_matrix()
    if on["surface"]:
        iso = build_surface(m, p)
        iso.points = iso.points @ R.T
        sc["surface"] = iso.compute_normals(auto_orient_normals=True)
    if on["open"] or on["closed"] or on["streamer"]:
        tracer = FieldTracer(m)
    if on["open"] or on["closed"]:
        lines = tracer.field_lines(fibonacci_sphere(p["n_seeds"], p["seed_r"]) @ R)
        n = {k: sum(kk == k for _, kk in lines) for k in (OPEN, CLOSED)}
        sc["open"] = clip([x @ R.T for x, k in lines if k == OPEN], p["line_r_max"])
        sc["closed"] = clip([x @ R.T for x, k in lines if k == CLOSED], p["line_r_max"])
        print(f"field lines: {n[OPEN]} open, {n[CLOSED]} closed, "
              f"{len(lines) - n[OPEN] - n[CLOSED]} undetermined (not drawn)")
    shell = shell_of(p["gcs"], p, p["time"]) if (on["gcs"] or on["streamer"]) else None
    if on["gcs"]:
        sc["gcs"] = shell
        if any(on[i] for i in WL) and p.get("gcs_wl"):
            sc["gcs"] = shell_of(p["gcs_wl"], p, p["gcs_wl"]["time"])   # shell at the image time
            sc["times"]["shell"] = p["gcs_wl"]["time"]
    if on["streamer"]:
        # closed lines threaded by the CME shell: traced from every vertex of the shell above
        # the inner boundary of the model
        pts = shell.reshape(-1, 3)
        lines = tracer.field_lines(pts[np.linalg.norm(pts, axis=1) > 1.02] @ R)
        n = {k: sum(kk == k for _, kk in lines) for k in (OPEN, CLOSED)}
        sc["streamer"] = clip([x @ R.T for x, k in lines if k == CLOSED], p["line_r_max"])
        print(f"lines through the CME shell: {n[CLOSED]} closed, {n[OPEN]} open, "
              f"{len(lines) - n[OPEN] - n[CLOSED]} undetermined, of {len(lines)}")
    if on["radio"]:
        sc["radio"], sc["radio_contours"] = radio_points(p), []
        rn = p["radio_nrh"]
        nrh_dir = Path(data_root) / rn["folder"]
        if any(nrh_dir.glob("nrh2_*")):
            from . import nrh
            from .whitelight import utc
            t0, t1 = utc(rn["start"]), utc(rn["end"])
            for j, f in enumerate(rn["freqs"]):
                h = newkirk_r(f)
                cmap = colormaps[rn["colormaps"][j % len(rn["colormaps"])]]
                blobs = []                                 # one per step between start and end
                t = t0
                while t <= t1:
                    img, tf, x, y = nrh.image_at(nrh_dir, f, t)
                    lines, (cx, cy), tb = nrh.source(img, x, y, rn["level"])
                    if tb >= rn["min_tb"] and (not blobs or blobs[-1][0] != tf):
                        blobs.append((tf, cx, cy, tb, lines))
                    t += timedelta(seconds=rn["step_s"])
                if not blobs:
                    print(f"radio: NRH {f:g} MHz has no image above {rn['min_tb']:.1e} K "
                          f"between {rn['start']} and {rn['end']}")
                    continue
                xy = np.array([(cx, cy) for _, cx, cy, _, _ in blobs])
                xyz = deproject(xy[:, 0], xy[:, 1], h, *p["earth"])
                colours = cmap(np.linspace(0, 1, len(blobs)))[:, :3]       # early to late
                sc["radio"].append((f, xyz, colours))
                if rn["contours"]:
                    sc["radio_contours"].append((f, [deproject(c[:, 0], c[:, 1], h, *p["earth"], warn=False)
                                                     for c in blobs[-1][4]]))
                sc["times"].setdefault("radio", []).append((f, blobs[0][0], blobs[-1][0], len(blobs)))
                for (tf, cx, cy, tb, _) in blobs:
                    print(f"radio: NRH {f:g} MHz {tf:%H:%M:%S} UT, peak {tb:.2e} K, centroid "
                          f"({cx:+.2f}, {cy:+.2f}) R_sun, at {h:.3f} R_sun")
        if not sc["radio"]:
            print("radio: no NRH images in", nrh_dir, "and no radio_sources given; layer empty")
    for inst in WL:
        if on[inst]:
            from .whitelight import INSTR, load_plane
            sc[inst], tf, tb = load_plane(inst, wl_dir(data_root), p["wl_time"], p[INSTR[inst]["obs"]],
                                          p["wl_window"], p["wl_mode"])
            sc["times"][inst] = (tf, tb)
    return sc


def polyline(list_of_pts):
    pts, cells, off = [], [], 0
    for a in list_of_pts:
        pts.append(a)
        cells.append(np.r_[len(a), np.arange(off, off + len(a))])
        off += len(a)
    poly = pv.PolyData(np.vstack(pts))
    poly.lines = np.concatenate(cells)
    return poly


def render(sc, on, p, cam_dir, up, frame, owner, out, name, base_half=None):
    pl = pv.Plotter(off_screen=True, window_size=p["window"])
    pl.set_background("white")
    if on["sphere"]:
        pl.add_mesh(pv.Sphere(radius=p["sphere_r"], theta_resolution=120, phi_resolution=120),
                    color="#9e9e9e", smooth_shading=True, specular=0.1)
    if on["surface"]:
        pl.add_mesh(sc["surface"], scalars="va", cmap="Spectral_r", clim=p["va_range"],
                    smooth_shading=True, specular=0.1, show_scalar_bar=False)
    for key, kw in (("open", dict(color="black", line_width=2.5, opacity=0.08)),
                    ("closed", dict(color="#8c8c8c", line_width=2.5, opacity=0.12)),
                    ("streamer", dict(color="#9a6fc0", line_width=3.0,
                                      opacity=0.35 if sc.get("radio") and on["radio"] else 0.7))):
        if on[key] and sc[key]:
            lines = sc[key]
            if key == "streamer" and len(lines) > p["streamer_max_lines"]:    # an even subset
                lines = [lines[i] for i in np.linspace(0, len(lines) - 1, p["streamer_max_lines"]).astype(int)]
            pl.add_mesh(polyline(lines), **kw)
    if on["gcs"]:
        pl.add_mesh(polyline(list(sc["gcs"])), color=p["mesh_rings"]["color"],
                    line_width=p["mesh_rings"]["width"])
    for f, lines in (sc.get("radio_contours", []) if on["radio"] else []):
        if lines:
            pl.add_mesh(polyline(lines), color=cm.cool(np.clip((f - 150) / 300, 0, 1))[:3], line_width=4.0)
    for f, xyz, rgb in (sc.get("radio", []) if on["radio"] else []):
        cloud = pv.PolyData(xyz)
        cloud["rgb"] = (np.asarray(rgb) * 255).astype(np.uint8)
        pl.add_mesh(cloud.glyph(geom=pv.Sphere(radius=p["radio_radius"]), scale=False, orient=False),
                    scalars="rgb", rgb=True, smooth_shading=True, ambient=0.35, specular=0.4)
    # coronagraph planes only in a view from their own spacecraft
    for inst in IMAGES_OF.get(owner, ()):
        if on[inst]:
            v = sc[inst]["wl"]
            lim = np.nanpercentile(np.abs(v), 98)
            pl.add_mesh(sc[inst], scalars="wl", cmap="gray", clim=(-lim, lim),
                        nan_opacity=0.0, opacity=p["wl_opacity"], show_scalar_bar=False,
                        lighting=False)
    half, dx, dy = frame
    right = np.cross(up, cam_dir)
    right /= np.linalg.norm(right)
    focal = -(dx * right + dy * np.asarray(up))      # puts Sun centre at (dx, dy) in the panel
    pl.camera.focal_point = tuple(focal)
    pl.camera.up = tuple(up)
    if p["projection"] == "perspective":
        # the camera of the published figure: distance fixed at the default framing, so that a
        # wider panel moves the camera back and keeps the perspective
        D = p["camera_distance"] * half / (base_half or half)
        pl.camera.position = tuple(focal + D * np.asarray(cam_dir))
        pl.camera.view_angle = float(np.degrees(2 * np.arctan(half / D)))
        pl.camera.clipping_range = (0.05, 10 * D)
    else:
        # orthographic projection: every depth is drawn at the scale of the plane of the sky, so
        # that the model, the shell, and the coronagraph planes share one scale in R_sun
        pl.enable_parallel_projection()
        pl.camera.position = tuple(focal + 100.0 * np.asarray(cam_dir))
        pl.camera.parallel_scale = half
        pl.camera.clipping_range = (1.0, 250.0)
    img = pl.screenshot(str(out / f"{name}.png"), return_img=True)
    pl.close()
    return img


def view_geometry(v, p):
    """Viewing direction (from Sun centre towards the viewer), screen up, and the observer that
    the view belongs to (earth, stereo_b, or None) for one entry of the configuration's views.
    An entry names an observer ("from": "earth", "stereo_b", or "north") or gives a Carrington
    longitude and latitude ("lon", "lat") of the viewing direction; "roll" turns the panel
    counter-clockwise by that angle [deg]."""
    src = v.get("from")
    if src in ("earth", "stereo_b"):
        lon, lat = p[src]
    elif src == "north":
        lon, lat = p["earth"][0], 90.0
    elif src is None and "lon" in v and "lat" in v:
        lon, lat = v["lon"], v["lat"]
    else:
        raise ValueError(f"view {v}: give 'from' (earth, stereo_b, north) or 'lon' and 'lat'")
    d = unit_vec(lon, lat)
    if abs(lat) > 89.9:                                  # from a pole: Earth's direction at the bottom
        e = unit_vec(*p["earth"])
        up = -np.sign(lat) * np.r_[e[:2], 0] / np.linalg.norm(e[:2])
    else:
        up = observer_basis(lon, lat)[1]                 # solar north up
    roll = np.radians(v.get("roll", 0.0))
    up = up * np.cos(roll) + np.cross(d, up) * np.sin(roll)
    return d, up, (src if src in ("earth", "stereo_b") else None)


def parse_view(text):
    """'b=north', 'a=earth', or 'c=300,20[,roll]' -> (panel index, view entry)."""
    k, val = text.split("=", 1)
    i = "abc".index(k.strip().lower())
    val = val.strip()
    if val in ("earth", "stereo_b", "north"):
        return i, {"from": val}
    nums = [float(x) for x in val.split(",")]
    return i, dict(lon=nums[0], lat=nums[1], roll=nums[2] if len(nums) > 2 else 0.0)


def make_figure(out="figures/corona3d.png", data_root=".", half=None, params=None, **layers):
    """Render the three-panel figure.

    data_root: working directory (its config.json, if present, is read); out: path of the
    figure, relative to data_root unless absolute; params: configuration keys that replace
    those of config.json (nested keys may be partial); layers: any of LAYERS as True/False;
    half: panel half-height [R_sun] for all panels (default: the framing of the configuration,
    or wide enough for the white-light images that are switched on).
    """
    bad = set(layers) - set(LAYERS)
    if bad:
        raise ValueError(f"unknown layers {sorted(bad)}; choose from {', '.join(LAYERS)}")
    if half is not None and half <= 0:
        raise ValueError("half must be positive")
    on = dict(ON, **layers)
    p = merge(config(data_root), params)
    out = Path(out) if Path(out).is_absolute() else Path(data_root) / out
    out.parent.mkdir(parents=True, exist_ok=True)
    sc = reconstruct(p, on, data_root)

    names = {"earth": "Earth", "north": "North Pole", "stereo_b": "STEREO B"}
    imgs, titles, owners = [], [], []
    for i, v in enumerate(p["views"]):
        d, up, owner = view_geometry(v, p)
        wide = (6.0 if on["c2"] else None) if owner == "earth" else \
            (6.0 if on["cor2"] else 4.0 if on["cor1"] else None) if owner == "stereo_b" else None
        fr = v.get("frame", [1.6, 0.0, 0.0])
        h = half if half is not None else (wide or fr[0])
        name = names.get(v.get("from")) or f"lon {v['lon']:g}°, lat {v['lat']:g}°"
        titles.append(v.get("title") or f"{'abc'[i]}. {name}")
        owners.append(owner)
        imgs.append(render(sc, on, p, d, up, (h, fr[1], fr[2]), owner, out.parent, f"view_{'abc'[i]}", fr[0]))

    t = sc["times"]
    notes = [[f"model {p['run'].split('/')[0]}, shell {str(t['shell'])[11:16]} UT"] for _ in imgs]
    for i in range(len(notes)):
        for f, ta, tb_, n in t.get("radio", []):
            notes[i].append(f"NRH {f:g} MHz, {n} images {ta:%H:%M:%S}-{tb_:%H:%M:%S} UT")
    for i, owner in enumerate(owners):
        for inst in IMAGES_OF.get(owner, ()):
            if inst in t:
                tf, tb = t[inst]
                notes[i].append(f"{inst.upper()} {tf:%H:%M} vs {tb:%H:%M} UT ({p['wl_mode']})")
    fig = Figure(figsize=(16, 4.9))           # no GUI backend: drawn off screen
    for i, (title, img) in enumerate(zip(titles, imgs)):
        ax = fig.add_axes([0.005 + i * 0.333, 0.17, 0.325, 0.81])
        ax.imshow(img)
        ax.set_xticks([]), ax.set_yticks([])
        lab, txt = title.split(" ", 1)
        ax.text(0.01, 0.985, lab, transform=ax.transAxes, fontsize=15, va="top")
        ax.text(0.06, 0.97, txt, transform=ax.transAxes, fontsize=15, va="top",
                bbox=dict(fc="white", ec="none", alpha=0.85))
        ax.text(0.99, 0.985, "\n".join(notes[i]), transform=ax.transAxes, fontsize=8, va="top",
                ha="right", bbox=dict(fc="white", ec="none", alpha=0.7))
    if on["surface"]:
        cax = fig.add_axes([0.40, 0.07, 0.24, 0.05])
        cb = fig.colorbar(cm.ScalarMappable(colors.Normalize(*p["va_range"]), "Spectral_r"),
                          cax=cax, orientation="horizontal")
        cb.ax.tick_params(labelsize=12)
        fig.text(0.395, 0.095, r"$v_A$ (km/s)", ha="right", va="center", fontsize=14)
    fig.savefig(out, dpi=150)
    print("saved", out)
    return out
