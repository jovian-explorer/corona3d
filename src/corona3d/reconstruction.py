"""The three-dimensional reconstruction and the three-panel figure, every layer switchable.

MAS/MAST n_e iso-surface coloured by v_A, open, closed, and streamer field lines, a grey
reference sphere, the GCS shell of the CME, radio sources at their Newkirk heights, and LASCO C2
and STEREO-B COR1, COR2 base-difference planes, seen from Earth, the solar north pole, and
STEREO-B (the layout of Fig. 7 of Morosan et al. 2022).
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from matplotlib import cm, colors

from . import mast_model
from .paths import default_config, model_dir, wl_dir
from .geometry import FieldTracer, deproject, fibonacci_sphere, observer_basis, pythea_shell, unit_vec

LAYERS = ("surface", "sphere", "open", "closed", "streamer", "gcs", "radio", "c2", "cor1", "cor2")
ON = dict(surface=True, sphere=True, open=True, closed=True, streamer=True, gcs=True,
          radio=False, c2=False, cor1=False, cor2=False)

# the event of 2013-10-25 (Morosan et al. 2022); a working directory's config.json overrides it
P = default_config()


def newkirk_r(f_mhz, fold=4.0, harmonic=True):
    fp = f_mhz / (2 if harmonic else 1) * 1e6
    ne = (fp / 8980.0) ** 2
    return 4.32 / np.log10(ne / (fold * 4.2e4))


def rot_z(deg):
    c, s = np.cos(np.radians(deg)), np.sin(np.radians(deg))
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])


def build_surface(m, p):
    PH, TH, R = np.meshgrid(m["phi"], m["theta"], m["r"], indexing="ij")
    grid = pv.StructuredGrid(R * np.sin(TH) * np.cos(PH), R * np.sin(TH) * np.sin(PH),
                             R * np.cos(TH))
    grid["ne"] = m["ne"].ravel(order="F")
    grid["va"] = m["va"].ravel(order="F")
    return grid.contour([p["ne_level"]], scalars="ne").clean()


def radio_points(p):
    """[(freq, xyz)] for the user-supplied sources, deprojected from the Earth view."""
    out = []
    for s in p["radio_sources"]:
        k = 1.0 / 959.6 if s.get("unit") == "arcsec" else 1.0      # arcsec -> R_sun at 1 AU
        out.append((s["freq"], deproject(np.atleast_1d(s["x"]) * k, np.atleast_1d(s["y"]) * k,
                                         newkirk_r(s["freq"]), *p["earth"])))
    return out


def reconstruct(p, on, data_root):
    """Only the layers switched on are built (field tracing is the slow part)."""
    sc = {}
    need_model = on["surface"] or on["open"] or on["closed"] or on["streamer"]
    if need_model:
        m = mast_model.load(density_scale=p["density_scale"], data_dir=model_dir(data_root, p))
        R = rot_z(p["model_dlon"])
    if on["surface"]:
        iso = build_surface(m, p)
        iso.points = iso.points @ R.T
        sc["surface"] = iso.compute_normals(auto_orient_normals=True)
    if on["open"] or on["closed"] or on["streamer"]:
        tracer = FieldTracer(m)
    if on["open"] or on["closed"]:
        lines = tracer.field_lines(fibonacci_sphere(p["n_seeds"], p["seed_r"]))
        sc["open"] = [x @ R.T for x, o in lines if o and len(x) > 2]
        sc["closed"] = [x @ R.T for x, o in lines if not o and len(x) > 2]
    g = p["gcs"]
    if on["gcs"] or on["streamer"]:
        shell = pythea_shell(g["height"], g["kappa"], g["alpha_deg"], g["lon_stonyhurst"],
                             g["lat_deg"], g["tilt_deg"] + 90.0, p["earth"][0], obstime=p["time"],
                             nbvertsl=p["mesh_rings"]["nbvertsl"],
                             nbvertcirc=p["mesh_rings"]["nbvertcirc"])
        sc["gcs"] = shell
    if on["streamer"]:
        sf = p["streamer_flank"]
        pts = shell.reshape(-1, 3)
        pts = pts[np.linalg.norm(pts, axis=1) > 1.2]
        f = pts[np.argmax(pts @ unit_vec(*p["earth"]))]            # Earth-facing flank point
        latf = np.degrees(np.arcsin(f[2] / np.linalg.norm(f)))
        lats = np.linspace(latf - sf["half_span"], latf + sf["half_span"], sf["n"])
        seeds = np.array([sf["r"] * unit_vec(p["earth"][0] + sf["lon"], la) for la in lats]) @ R
        sc["streamer"] = [x @ R.T for x, o in tracer.field_lines(seeds)
                          if not o and np.linalg.norm(x, axis=1).max() > p["streamer_min_apex"]]
    if on["radio"]:
        sc["radio"] = radio_points(p)
        if not sc["radio"]:
            print("radio: no radio_sources given, layer empty")
    for inst in ("c2", "cor1", "cor2"):
        if on[inst]:
            from .whitelight import INSTR, load_plane
            sc[inst] = load_plane(inst, wl_dir(data_root), p["wl_time"],
                                  p[INSTR[inst]["obs"]])
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


def render(sc, on, p, cam_dir, up, frame, view, out):
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
                    ("streamer", dict(color="#9a6fc0", line_width=4.0, opacity=0.85))):
        if on[key] and sc[key]:
            pl.add_mesh(polyline(sc[key]), **kw)
    if on["gcs"]:
        pl.add_mesh(polyline(list(sc["gcs"])), color=p["mesh_rings"]["color"],
                    line_width=p["mesh_rings"]["width"])
    for f, xyz in sc.get("radio", []):
        c = cm.cool(np.clip((f - 150) / 300, 0, 1))[:3]
        pl.add_mesh(pv.PolyData(xyz).glyph(geom=pv.Sphere(radius=p["radio_radius"]),
                                           scale=False, orient=False),
                    color=c, smooth_shading=True)
    # white-light planes: C2 in the Earth view, COR1/COR2-B in the STEREO-B view
    for inst in {"earth": ("c2",), "stb": ("cor2", "cor1")}.get(view, ()):
        if inst in sc:
            v = sc[inst]["wl"]
            lim = np.nanpercentile(np.abs(v), 98)
            pl.add_mesh(sc[inst], scalars="wl", cmap="gray", clim=(-lim, lim),
                        nan_opacity=0.0, opacity=p["wl_opacity"], show_scalar_bar=False,
                        lighting=False)
    half, dx, dy = frame
    right = np.cross(up, cam_dir)
    right /= np.linalg.norm(right)
    focal = -(dx * right + dy * np.asarray(up))
    D = p["camera_distance"] * half / p["frame"][view][0]     # wider FOV: camera pulled back
    pl.camera.focal_point = tuple(focal)
    pl.camera.up = tuple(up)
    pl.camera.position = tuple(focal + D * cam_dir)
    pl.camera.view_angle = float(np.degrees(2 * np.arctan(half / D)))
    pl.camera.clipping_range = (0.1, 10 * D)
    img = pl.screenshot(str(out / f"view_{view}.png"), return_img=True)
    pl.close()
    return img


def make_figure(out="figures/corona3d.png", data_root=".", half=None, params=None, **layers):
    """Render the three-panel figure.

    data_root: working directory holding data/ (see `corona3d fetch`); params: configuration
    keys that replace those of P; layers: any of LAYERS as True/False; half: panel half-height
    [R_sun] for all panels (default: the framing of the paper, or wide enough for the
    white-light images that are switched on).
    """
    on = dict(ON, **layers)
    bad = set(layers) - set(LAYERS)
    if bad:
        raise ValueError(f"unknown layers {bad}; choose from {LAYERS}")
    p = dict(P, **(params or {}))
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    sc = reconstruct(p, on, data_root)

    e = unit_vec(*p["earth"])
    _, up_e, _ = observer_basis(*p["earth"])
    _, up_b, _ = observer_basis(*p["stereo_b"])
    away = -np.r_[e[:2], 0] / np.linalg.norm(e[:2])
    views = [("a. Earth", e, up_e, "earth"),
             ("b. North Pole", np.array([0, 0, 1.0]), away, "north"),
             ("c. STEREO B", unit_vec(*p["stereo_b"]), up_b, "stb")]
    wide = {"earth": 6.0 if on["c2"] else None,
            "stb": 6.0 if on["cor2"] else 4.0 if on["cor1"] else None}
    imgs = []
    for _, d, u, n in views:
        fr = p["frame"][n]
        h = half or wide.get(n) or fr[0]
        imgs.append(render(sc, on, p, d, u, (h, fr[1], fr[2]), n, out.parent))

    fig = plt.figure(figsize=(16, 5.6))
    for i, ((title, *_), img) in enumerate(zip(views, imgs)):
        ax = fig.add_axes([0.005 + i * 0.333, 0.17, 0.325, 0.81])
        ax.imshow(img)
        ax.set_xticks([]), ax.set_yticks([])
        lab, txt = title.split(" ", 1)
        ax.text(0.01, 0.985, lab, transform=ax.transAxes, fontsize=15, va="top")
        ax.text(0.06, 0.97, txt, transform=ax.transAxes, fontsize=15, va="top",
                bbox=dict(fc="white", ec="none", alpha=0.85))
    if on["surface"]:
        cax = fig.add_axes([0.40, 0.07, 0.24, 0.05])
        cb = fig.colorbar(cm.ScalarMappable(colors.Normalize(*p["va_range"]), "Spectral_r"),
                          cax=cax, orientation="horizontal")
        cb.set_ticks(np.arange(100, 701, 100))
        cb.ax.tick_params(labelsize=12)
        fig.text(0.395, 0.095, r"$v_A$ (km/s)", ha="right", va="center", fontsize=14)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("saved", out)
    return out
