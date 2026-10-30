"""The three-dimensional reconstruction and the three-panel figure, every layer switchable.

MAS/MAST n_e iso-surface coloured by v_A, open and closed field lines, the closed lines threaded
by the CME shell (the streamer loops), a grey reference sphere, the GCS shell of the CME, radio
sources at their Newkirk heights as blobs and contours, and LASCO C2, STEREO COR1, COR2 and Solar
Orbiter/Metis images. Each panel is seen from an observer (Earth, the solar north pole, the STEREO
spacecraft of the event, or Solar Orbiter) or from any Carrington direction. For the validation event
of 2013 October 25 the default panels share the viewpoints of Fig. 7 of Morosan et al. (2022), so that
the two reconstructions can be compared panel by panel.
"""
from datetime import timedelta
from pathlib import Path

import astropy.units as u
import numpy as np
from matplotlib import cm, colors
from matplotlib.figure import Figure
from matplotlib import colormaps
from scipy.spatial.transform import Rotation

from . import mast_model
from .geometry import CLOSED, OPEN, FieldTracer, deproject, fibonacci_sphere, observer_basis, pythea_shell, unit_vec
from ._extras import need
from .layers import LAYERS, ON
from .paths import config, default_config, merge, model_dir, wl_dir

pv = need("render", "pyvista", "the drawing of the figure")

WL = ("c2", "cor1", "cor2", "metis")
IMAGES_OF = {"earth": ("c2",), "stereo": ("cor2", "cor1"), "solo": ("metis",)}
OBSERVERS = ("earth", "stereo", "solo")               # observers with a position in the configuration


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


def chain_colours(base, frac):
    """Colours of the blobs of one frequency: its base colour, light at the start of the radio time
    range (frac 0), the base colour at the middle, and dark at the end (frac 1)."""
    rgb = np.asarray(colors.to_rgb(base))
    f = np.clip(np.asarray(frac, float), 0, 1)[:, None]
    light, dark = rgb + (1 - rgb) * 0.55, rgb * 0.5
    return np.where(f < 0.5, light + (rgb - light) * 2 * f, rgb + (dark - rgb) * (2 * f - 1))


def clip(lines, r_max):
    """Parts of each line inside r_max [R_sun]: lines are traced to the outer boundary of the
    model for their classification but drawn only near the Sun."""
    out = []
    for x in lines:
        inside = np.linalg.norm(x, axis=1) <= r_max
        cuts = np.flatnonzero(np.diff(inside.astype(int))) + 1
        out += [seg for seg, ok in zip(np.split(x, cuts), np.split(inside, cuts)) if ok[0] and len(seg) > 1]
    return out


SHELL_COLOR = {"mach_a": ("M_A", "viridis", (0.0, 4.0)), "mach_f": ("M_f", "viridis", (0.0, 4.0)),
               "theta_bn": (r"$	heta_{Bn}$ (deg)", "coolwarm", (0.0, 90.0)), "vn": ("$V_n$ (km/s)", "magma", None),
               "va": (r"$v_A$ (km/s)", "Spectral_r", None)}


def _model(p, data_root, m):
    return m if m is not None else mast_model.load(model_dir(data_root, p), density_scale=p["density_scale"],
                                                   va_density=p["va_density"])


def shell_scalars(p, rings, g, data_root, m=None):
    """(label, values per vertex of the drawn shell, colour map, limits) for `shell_color`: the Alfven or
    fast-magnetosonic Mach number, theta_Bn, the normal speed, or the Alfven speed of the upstream corona
    (corona3d.shock). Vertices that do not advance, or lie below shock.r_min, are NaN."""
    from . import shock
    key = p["shell_color"]
    if key not in SHELL_COLOR:
        raise ValueError(f"shell_color {key!r}: choose from {', '.join(SHELL_COLOR)}")
    s = shock.settings(p)
    m = _model(p, data_root, m)
    mast_model.load_plasma(model_dir(data_root, p), m)
    v, _, _ = shock.apex_speed(data_root, p, p.get("shock_speed_km_s"), "shock_speed_km_s")
    pts = rings.reshape(-1, 3)
    ok = np.linalg.norm(pts, axis=1) > s["r_min"]                  # the legs reach Sun centre
    d = shock.diagnostics(pts[ok], shock.shell_normals(rings).reshape(-1, 3)[ok],
                          shock.Plasma(m, p, s["temperature_k"], s["flow"]), v, float(g["height"]), s["gamma"])
    vals = np.full(len(pts), np.nan)
    vals[ok] = np.where(d["vn"] > 0, d[key], np.nan) if key == "theta_bn" else d[key]
    label, cmap, lim = SHELL_COLOR[key]
    lim = tuple(p.get("shell_clim") or lim or ((0.0, v) if key == "vn" else tuple(p["va_range"])))
    return label, vals, cmap, lim


def height_segments(p, sc, data_root, m=None):
    """Segments along the line of sight from Earth between the lowest and highest height of every radio
    position over the density models of shock.height_models (corona3d.heights)."""
    from . import heights, shock
    models = shock.settings(p)["height_models"]
    rows = [(s_["freq"], x, y) for s_ in p["radio_sources"] for x, y in zip(np.atleast_1d(s_["x"]), np.atleast_1d(s_["y"]))
            if s_.get("unit") != "arcsec"]
    rows += [(b["freq"], x, y) for b in sc.get("radio_blobs", []) for x, y in zip(b["x"], b["y"])]
    if not rows:
        return []
    f, x, y = (np.array(c, float) for c in zip(*rows))
    r = heights.ranges(dict(freq_mhz=f, x_rsun=x, y_rsun=y), p, _model(p, data_root, m) if "mas" in models else None, models)
    return [np.array([[a, b, c], [d, e, g]]) for a, b, c, d, e, g in
            zip(r["x_min"], r["y_min"], r["z_min"], r["x_max"], r["y_max"], r["z_max"]) if np.isfinite([a, d]).all()]


def shell_of(g, p, time):
    return pythea_shell(g["height"], g["kappa"], g["alpha_deg"], g["lon_stonyhurst"], g["lat_deg"],
                        g["tilt_deg"] + 90.0, p["earth"][0], obstime=time,
                        nbvertsl=p["mesh_rings"]["nbvertsl"], nbvertcirc=p["mesh_rings"]["nbvertcirc"])


def burst_shell(p):
    """Configuration for a coronagraph figure at the burst time (`wl_at` "burst"): the shell is the one
    at t_out (`gcs` of a fitted configuration, `gcs_wl` of one without `gcs_fit`) placed at the time of
    the event, and wl_time is the time of the event unless it was given. A height given with `wl_height`
    replaces the fitted one; it is called assumed, with `wl_height_note`, as is any height that the fit
    itself marks as assumed. Returns a copy of p with `gcs_wl`, `wl_time` and `wl_shell` set."""
    p = dict(p)
    fit = p.get("gcs_fit")
    h_given = p.get("wl_height")
    if h_given is not None and not p.get("wl_height_note"):
        raise ValueError("--wl-height needs --wl-height-note, the reason for the assumed height")
    if fit:
        g, assumed, note = dict(p["gcs"]), fit.get("height_source") == "assumed", fit.get("height_note", "")
    elif h_given is not None:
        g, assumed, note = dict(p["gcs_wl"]), False, ""
    else:
        raise ValueError("the configuration has no gcs_fit, so the height of the shell at the burst time is not known: "
                         "give --wl-height and --wl-height-note")
    if h_given is not None:
        g["height"], assumed, note = float(h_given), True, p["wl_height_note"]
    g["time"] = p["time"]
    p["gcs_wl"] = g
    p["wl_time"] = p.get("wl_time_given") or p["time"]
    p["wl_shell"] = dict(height=g["height"], assumed=assumed, note=note, time=p["time"])
    return p


def reconstruct(p=None, on=None, data_root="."):
    """Build the layers that are switched on; p and on may hold only the keys that differ from
    the defaults. Returns a dictionary of the layers and of the times used."""
    from .paths import colormap_rule
    p = merge(default_config(), colormap_rule(p or {}))    # radio_colormaps without radio_colors: as for a config.json
    from .pb import resolve
    p = resolve(p)                                         # density_scale "pb": the scale fitted by corona3d pb
    if p.get("wl_at") == "burst":
        p = burst_shell(p)
    on = dict(ON, **(on or {}))
    sc = {"times": {"shell": p["time"]}}
    m = None
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
            if p.get("wl_shell"):
                s = p["wl_shell"]
                sc["times"]["shell_height"] = f"height {s['height']:g} R_sun" + (", assumed" if s["assumed"] else "")
    if on["gcs"] and p.get("shell_color"):
        g = p["gcs_wl"] if sc["times"]["shell"] != p["time"] else p["gcs"]
        sc["gcs_scalars"] = shell_scalars(p, sc["gcs"], g, data_root, m)
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
        sc["radio"], sc["radio_contours"], sc["radio_blobs"] = radio_points(p), [], []
        from . import nrh, ovro_lwa
        from .whitelight import utc
        readers = {"nrh": nrh, "ovro_lwa": ovro_lwa}
        from . import heights
        ne_cache = []
        ne_loader = lambda: ne_cache[0] if ne_cache else ne_cache.append(heights.ne_sampler(_model(p, data_root, m), p)) or ne_cache[0]
        names = {"nrh": "NRH", "ovro_lwa": "OVRO-LWA"}
        j = 0                                              # colours run over all frequencies
        # one time scale for every chain, fixed by radio_time_range (default: the start and end of the
        # entries), so that a blob keeps its colour when the interval is extended, as in an animation
        ri = p["radio_images"]
        tr = p.get("radio_time_range") or ([min((rs["start"] for rs in ri), key=utc), max((rs["end"] for rs in ri), key=utc)]
                                           if ri else [p["time"], p["time"]])
        T0, T1 = utc(tr[0]), utc(tr[1])
        span = max((T1 - T0).total_seconds(), 1.0)
        sc["radio_time"], sc["radio_legend"], sc["radio_contour_colours"] = (T0, T1), [], []
        base_of = {}                                       # one colour per (instrument, frequency)
        for rs in p["radio_images"]:
            mod, folder = readers[rs["instrument"]], Path(data_root) / rs["folder"]
            if not folder.exists() or not any(folder.iterdir()):
                print(f"radio: no {names[rs['instrument']]} files in {folder}")
                j += len(rs["freqs"])                      # the other frequencies keep their colours
                continue
            t0, t1 = utc(rs["start"]), utc(rs["end"])
            for f in rs["freqs"]:
                h = newkirk_r(f)
                label = f"{names[rs['instrument']]} {f:g} MHz"
                if p.get("radio_colors"):
                    base = base_of.setdefault((rs["instrument"], f), p["radio_colors"][len(base_of) % len(p["radio_colors"])])
                    cmap = None
                else:
                    base, cmap = None, colormaps[p["radio_colormaps"][j % len(p["radio_colormaps"])]]
                j += 1
                blobs, last = [], None                     # one per step between start and end
                t = t0
                while t <= t1:
                    img, tf, x, y = mod.image_at(folder, f, t)
                    lines, (cx, cy), tb = nrh.source(img, x, y, rs["level"])
                    if tb >= rs["min_tb"] and (not blobs or blobs[-1][0] != tf):
                        blobs.append((tf, cx, cy, tb, lines))
                        last = (img, x, y)
                    t += timedelta(seconds=rs["step_s"])
                if not blobs:
                    print(f"radio: {label} has no image above {rs['min_tb']:.1e} K between {rs['start']} and {rs['end']}")
                    continue
                xy = np.array([(cx, cy) for _, cx, cy, _, _ in blobs])
                xyz = deproject(xy[:, 0], xy[:, 1], h, *p["earth"])
                if heights.custom(p, rs):                  # model surface, F emission, scattering, or offset (heights.place_blobs)
                    xyz, r_each, _ = heights.place_blobs(rs, f, xy, [tf for tf, *_ in blobs], p, data_root, ne_loader)
                    h = float(np.median(r_each))
                frac = [(tf - T0).total_seconds() / span for tf, *_ in blobs]     # early to late
                colours = chain_colours(base, frac) if base else cmap(np.clip(frac, 0, 1))[:, :3]
                if base and (label, base) not in sc["radio_legend"]:
                    sc["radio_legend"].append((label, base))
                sc["radio"].append((f, xyz, colours))
                sc["radio_blobs"].append(dict(instrument=rs["instrument"], freq=f, r=h, xyz=xyz,
                                              time=[tf for tf, *_ in blobs], x=xy[:, 0], y=xy[:, 1],
                                              tb=np.array([b[3] for b in blobs])))
                if rs["contours"]:                         # of the last image, at each level of its peak
                    sc["radio_contour_colours"].append(colors.to_rgb(base) if base else cm.cool(np.clip((f - 150) / 300, 0, 1))[:3])
                    levels = sorted(rs.get("contour_levels") or [rs["level"]])
                    sc["radio_contours"].append((f, [(lev, [deproject(c[:, 0], c[:, 1], h, *p["earth"], warn=False)
                                                            for c in nrh.source(*last, lev)[0]])
                                                     for lev in levels]))
                sc["times"].setdefault("radio", []).append((label, blobs[0][0], blobs[-1][0], len(blobs)))
                for (tf, cx, cy, tb, _) in blobs:
                    print(f"radio: {label} {tf:%H:%M:%S} UT, peak {tb:.2e} K, centroid "
                          f"({cx:+.2f}, {cy:+.2f}) R_sun, at {h:.3f} R_sun")
        if not sc["radio"]:
            print("radio: no radio images and no radio_sources; layer empty")
        if p.get("radio_height_range") and sc["radio"]:
            sc["radio_ranges"] = height_segments(p, sc, data_root, m)
    for inst in WL:
        if on[inst]:
            from .whitelight import INSTR, WL_MAX_GAP_MIN, load_plane, utc
            if INSTR[inst]["obs"] not in p:
                raise ValueError(f"{inst}: the configuration has no position of {INSTR[inst]['obs']}; run corona3d fetch --solo")
            burst = p.get("wl_at") == "burst"
            sc[inst], tf, tb = load_plane(inst, wl_dir(data_root), p["wl_time"], p[INSTR[inst]["obs"]],
                                          p["wl_window"], p["wl_mode"], p["stereo_sc"], nearest=burst)
            sc["times"][inst] = (tf, tb)
            gap = (tf - utc(p["wl_time"])).total_seconds() / 60.0
            if burst and abs(gap) > WL_MAX_GAP_MIN:        # a frame far from the burst: say so on the figure
                sc["times"].setdefault("wl_gap", {})[inst] = gap
                print(f"{inst}: nearest frame, {abs(gap):.0f} min {'after' if gap > 0 else 'before'} the burst")
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


VARIABLE = ("surface", "wl")                         # opacity may vary with the value of the pixel
PER_ITEM = ("radio", "radio_contours")               # opacity may vary along the chain or over the levels


def check_opacity(op):
    """Raise ValueError unless every opacity is a number from 0 to 1, a list of such numbers (for
    the layers in VARIABLE and PER_ITEM), or the name of a PyVista transfer function such as
    "linear" or "sigmoid" (for the layers in VARIABLE)."""
    for k, v in op.items():
        if isinstance(v, str) and k in VARIABLE:
            continue
        vals = v if isinstance(v, (list, tuple)) and k in VARIABLE + PER_ITEM else [v]
        if not len(vals) or not all(isinstance(x, (int, float)) and not isinstance(x, bool) and 0 <= x <= 1
                                    for x in vals):
            kinds = ("a number from 0 to 1" + (", a list of them" if k in VARIABLE + PER_ITEM else "")
                     + (", or a transfer function such as linear" if k in VARIABLE else ""))
            raise ValueError(f"opacity of {k}: {v!r}; give {kinds}")


def opacities(value, n):
    """n opacities from a number (all equal) or a list, which is stretched over the n items when its
    length differs from n (the first value to the first item, the last value to the last)."""
    v = np.atleast_1d(np.asarray(value, float))
    if len(v) == n or len(v) == 1:
        return np.broadcast_to(v, (n,)).copy()
    return np.interp(np.linspace(0, 1, n), np.linspace(0, 1, len(v)), v)


def render(sc, on, p, cam_dir, up, frame, owner, out, name, base_half=None):
    op = p["opacity"]
    pl = pv.Plotter(off_screen=True, window_size=p["window"])
    pl.set_background("white")
    if on["sphere"]:
        pl.add_mesh(pv.Sphere(radius=p["sphere_r"], theta_resolution=120, phi_resolution=120),
                    color="#9e9e9e", smooth_shading=True, specular=0.1, opacity=op["sphere"])
    if on["surface"]:
        pl.add_mesh(sc["surface"], scalars="va", cmap="Spectral_r", clim=p["va_range"],
                    smooth_shading=True, specular=0.1, show_scalar_bar=False, opacity=op["surface"])
    for key, kw in (("open", dict(color="black", line_width=2.5, opacity=op["open"])),
                    ("closed", dict(color="#8c8c8c", line_width=2.5, opacity=op["closed"])),
                    ("streamer", dict(color="#9a6fc0", line_width=3.0,
                                      opacity=op["streamer_with_radio"] if sc.get("radio") and on["radio"]
                                      else op["streamer"]))):
        if on[key] and sc[key]:
            lines = sc[key]
            if key == "streamer" and len(lines) > p["streamer_max_lines"]:    # an even subset
                lines = [lines[i] for i in np.linspace(0, len(lines) - 1, p["streamer_max_lines"]).astype(int)]
            pl.add_mesh(polyline(lines), **kw)
    if on["gcs"] and sc.get("gcs_scalars"):                # the shell coloured by a shock quantity
        label, vals, cmap, lim = sc["gcs_scalars"]
        poly = polyline(list(sc["gcs"]))
        poly["shock"] = vals
        shell = pl.add_mesh(poly, scalars="shock", cmap=cmap, clim=lim, nan_color="#b0b0b0", show_scalar_bar=False,
                            line_width=p["mesh_rings"]["width"], opacity=op["gcs"])
    elif on["gcs"]:
        shell = pl.add_mesh(polyline(list(sc["gcs"])), color=p["mesh_rings"]["color"],
                            line_width=p["mesh_rings"]["width"], opacity=op["gcs"])
    for (f, levels), rgb in zip(sc.get("radio_contours", []) if on["radio"] else [], sc.get("radio_contour_colours", [])):
        # the levels run from the lowest to the highest, and so do the opacities
        for (lev, lines), a in zip(levels, opacities(op["radio_contours"], len(levels))):
            if lines:
                pl.add_mesh(polyline(lines), color=rgb, line_width=4.0, opacity=float(a))
    if on["radio"] and sc.get("radio_ranges"):             # height range of each source over the density models
        pl.add_mesh(polyline(sc["radio_ranges"]), color="#303030", line_width=3.0, opacity=0.9)
    blobs = []
    for f, xyz, rgb in (sc.get("radio", []) if on["radio"] else []):
        cloud = pv.PolyData(xyz)
        alpha = opacities(op["radio"], len(xyz))                # from the first image to the last
        cloud["rgb"] = (np.asarray(rgb) * 255).astype(np.uint8)
        if (alpha < 1).any():                                   # one alpha per blob, carried by its colour
            cloud["rgb"] = np.c_[cloud["rgb"], (alpha * 255).round().astype(np.uint8)]
        glyphs = cloud.glyph(geom=pv.Sphere(radius=p["radio_radius"]), scale=False, orient=False)
        # flatter lighting keeps the colour (frequency) and its lightness (time) readable, and a black
        # outline separates the blobs from the model surface beneath them
        blobs += [pl.add_mesh(glyphs, scalars="rgb", rgb=True, smooth_shading=True, ambient=0.5, diffuse=0.6,
                              specular=0.1)]
        if (alpha < 1).any():                                   # each outline fades with its blob
            blobs += [pl.add_silhouette(pv.Sphere(radius=p["radio_radius"], center=c), color="black", line_width=1.5,
                                        opacity=float(a)) for c, a in zip(xyz, alpha) if a > 0]
        else:
            blobs += [pl.add_silhouette(glyphs, color="black", line_width=1.5)]
    if on["gcs"] and p.get("shell_on_top"):
        # a second layer that shares the camera, so that the shell stays visible where the
        # surface of the model would hide it (a low shell at the time of the burst)
        from vtkmodules.vtkRenderingCore import vtkRenderer
        top = vtkRenderer()
        top.SetLayer(1)
        top.SetActiveCamera(pl.renderer.GetActiveCamera())
        top.InteractiveOff()
        pl.ren_win.SetNumberOfLayers(2)
        pl.ren_win.AddRenderer(top)
        for actor in [shell] + blobs:                 # the radio sources keep their depth order with the shell
            pl.renderer.RemoveActor(actor)
            top.AddActor(actor)
    # coronagraph planes only in a view from their own spacecraft
    for inst in IMAGES_OF.get(owner, ()):
        if on[inst]:
            v = sc[inst]["wl"]
            lim = np.nanpercentile(np.abs(v), 98)
            pl.add_mesh(sc[inst], scalars="wl", cmap="gray", clim=(-lim, lim),
                        nan_opacity=0.0, opacity=p.get("wl_opacity", op["wl"]), show_scalar_bar=False,
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
    the view belongs to (earth, stereo, solo, or None) for one entry of the configuration's views.
    An entry names an observer ("from": "earth", "stereo", "solo", or "north") or gives a Carrington
    longitude and latitude ("lon", "lat") of the viewing direction; "roll" turns the panel
    counter-clockwise by that angle [deg]."""
    src = v.get("from")
    if src in OBSERVERS:
        if src not in p:
            raise ValueError(f"view {v}: the configuration has no position of {src}; for Solar Orbiter run corona3d fetch --solo")
        lon, lat = p[src]
    elif src == "north":
        lon, lat = p["earth"][0], 90.0
    elif src is None and "lon" in v and "lat" in v:
        lon, lat = v["lon"], v["lat"]
    else:
        raise ValueError(f"view {v}: give 'from' (earth, stereo, solo, north) or 'lon' and 'lat'")
    d = unit_vec(lon, lat)
    if abs(lat) > 89.9:                                  # from a pole: Earth's direction at the bottom
        e = unit_vec(*p["earth"])
        up = -np.sign(lat) * np.r_[e[:2], 0] / np.linalg.norm(e[:2])
    else:
        up = observer_basis(lon, lat)[1]                 # solar north up
    roll = np.radians(v.get("roll", 0.0))
    up = up * np.cos(roll) + np.cross(d, up) * np.sin(roll)
    return d, up, (src if src in OBSERVERS else None)


def parse_view(text):
    """'b=north', 'a=earth', 'c=stereo', 'c=solo' (also the legacy 'c=stereo_a' and 'c=stereo_b'), or 'c=300,20[,roll]' -> (panel index, view entry)."""
    k, val = text.split("=", 1)
    i = "abc".index(k.strip().lower())
    val = val.strip()
    if val in ("stereo_a", "stereo_b"):                # legacy names: the STEREO spacecraft
        return i, {"from": "stereo", "sc": val[-1].upper()}
    if val in ("earth", "stereo", "north", "solo"):
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
    check_opacity(dict(p["opacity"], **({"wl": p["wl_opacity"]} if "wl_opacity" in p else {})))
    out = Path(out) if Path(out).is_absolute() else Path(data_root) / out
    out.parent.mkdir(parents=True, exist_ok=True)
    sc = reconstruct(p, on, data_root)

    names = {"earth": "Earth", "north": "North Pole", "stereo": f"STEREO {p['stereo_sc']}", "solo": "Solar Orbiter"}
    imgs, titles, owners = [], [], []
    for i, v in enumerate(p["views"]):
        d, up, owner = view_geometry(v, p)
        wide = (6.0 if on["c2"] else None) if owner == "earth" else \
            (6.0 if on["cor2"] else 4.0 if on["cor1"] else None) if owner == "stereo" else None
        if owner == "solo" and on["metis"]:              # the outer edge of the Metis image, which scales with distance
            wide = float(np.linalg.norm(sc["metis"].points[np.isfinite(sc["metis"]["wl"])], axis=1).max())
        fr = v.get("frame", [1.6, 0.0, 0.0])
        h = half if half is not None else (wide or fr[0])
        name = names.get(v.get("from")) or f"lon {v['lon']:g}°, lat {v['lat']:g}°"
        titles.append(v.get("title") or f"{'abc'[i]}. {name}")
        owners.append(owner)
        imgs.append(render(sc, on, p, d, up, (h, fr[1], fr[2]), owner, out.parent, f"view_{'abc'[i]}", fr[0]))

    t = sc["times"]
    height = f", {t['shell_height']}" if t.get("shell_height") else ""
    notes = [[f"model {p['run'].split('/')[0]}, shell {str(t['shell'])[11:16]} UT{height}"] for _ in imgs]
    for i in range(len(notes)):
        for label, ta, tb_, n in t.get("radio", []):
            notes[i].append(f"{label}, {n} images {ta:%H:%M:%S}-{tb_:%H:%M:%S} UT")
    for i, owner in enumerate(owners):
        for inst in IMAGES_OF.get(owner, ()):
            if inst in t:
                tf, tb = t[inst]
                notes[i].append(f"{inst.upper()} {tf:%H:%M} vs {tb:%H:%M} UT ({p['wl_mode']})")
                gap = (t.get("wl_gap") or {}).get(inst)
                if gap is not None:
                    notes[i][-1] += f", nearest frame, {abs(gap):.0f} min {'after' if gap > 0 else 'before'} the burst"
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
    if on["radio"] and sc.get("radio_legend"):
        # which colour is which frequency, and how the lightness of a blob gives the time of its image
        from matplotlib.lines import Line2D
        hs = [Line2D([], [], ls="", marker="o", ms=9, mfc=c, mec="black", label=lab) for lab, c in sc["radio_legend"]]
        fig.legend(handles=hs, loc="lower left", bbox_to_anchor=(0.003, 0.0), fontsize=8, frameon=False,
                   handletextpad=0.1, labelspacing=0.25, borderaxespad=0.3)
        T0, T1 = sc["radio_time"]
        cax = fig.add_axes([0.155, 0.07, 0.1, 0.05])
        grey = colors.ListedColormap(chain_colours("#808080", np.linspace(0, 1, 256)))
        cb = fig.colorbar(cm.ScalarMappable(colors.Normalize(0, 1), grey), cax=cax, orientation="horizontal",
                          ticks=[0, 1])
        cb.ax.set_xticklabels([f"{T0:%H:%M:%S}", f"{T1:%H:%M:%S}"])
        cb.ax.tick_params(labelsize=10)
        fig.text(0.205, 0.135, "radio image time (UT)", ha="center", fontsize=10)
    if on["gcs"] and sc.get("gcs_scalars"):
        label, _, cmap, lim = sc["gcs_scalars"]
        cax = fig.add_axes([0.73, 0.07, 0.18, 0.05])
        cb = fig.colorbar(cm.ScalarMappable(colors.Normalize(*lim), cmap), cax=cax, orientation="horizontal")
        cb.ax.tick_params(labelsize=12)
        fig.text(0.725, 0.095, f"shell {label}", ha="right", va="center", fontsize=14)
    fig.savefig(out, dpi=150)
    print("saved", out)
    return out
