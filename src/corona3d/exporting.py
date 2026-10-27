"""Write the layers of a reconstruction to data files, with a JSON sidecar of their provenance.

Formats: ``npz`` (one NumPy archive of every layer), ``csv`` (one table per layer), and ``vtk``
(PyVista/VTK files for ParaView: ``.vtp`` for the surface, lines, shells, and radio sources, and
``.vts`` for the coronagraph planes). Every format is accompanied by ``metadata.json``. Positions
are in R_sun, Heliographic Carrington Cartesian (see ``corona3d.api``).
"""
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from . import __version__
from .api import FRAME, RADIO_COLUMNS, UNITS, Reconstruction, build, carrington_angles

__all__ = ["export", "metadata", "FORMATS", "LINE_CLASSES"]

FORMATS = ("npz", "csv", "vtk")
LINE_CLASSES = {"open": 0, "closed": 1, "streamer": 2}      # class code in the vtk files


def _inputs(rec, home):
    """Sources of the inputs that the enabled layers used."""
    from . import fetch, ovro_lwa
    from .paths import model_dir, wl_dir
    p, on, out = rec.config, rec.layers, {}
    if on["surface"] or on["open"] or on["closed"] or on["streamer"]:
        d = model_dir(home, p)
        out["mhd_model"] = dict(
            name="PSI MAS/MAST " + p["run"], url=fetch.PSI.format(run=p["run"], var="{var}"),
            files={f.name: f.stat().st_size for f in sorted(d.glob("*.hdf"))},
            density_scale=p["density_scale"], va_density=p["va_density"], model_dlon_deg=p["model_dlon"])
    if on["gcs"]:
        out["gcs"] = dict(shell=p["gcs"], shell_image_time=p.get("gcs_wl"), note=p.get("gcs_note"),
                          fit=p.get("gcs_fit"), mesh=p["mesh_rings"])
    if on["radio"]:
        sets = []
        for r in p["radio_images"]:
            folder = Path(home) / r["folder"]
            files = sorted(f.name for f in folder.iterdir()) if folder.exists() else []
            sets.append(dict(r, files=files, source=ovro_lwa.S3 if r["instrument"] == "ovro_lwa"
                             else "Nancay Radioheliograph, level-2 images (obtained from the observatory)"))
        out["radio"] = dict(image_sets=sets, user_sources=p["radio_sources"],
                            height_model="four-fold Newkirk (1961), harmonic emission")
    for inst in rec.whitelight:
        out[inst] = dict(folder=str(wl_dir(home)), source="Virtual Solar Observatory (corona3d fetch --wl)",
                         mode=p["wl_mode"], window=p["wl_window"], **{k: rec.whitelight[inst][k]
                                                                      for k in ("time", "base_time")})
    return out


def metadata(rec, home="."):
    """The provenance record of a reconstruction, as a dictionary that is written as JSON.

    Parameters
    ----------
    rec : Reconstruction
        The layers from ``corona3d.api.build``.
    home : str or Path
        Working directory (to list the input files).

    Returns
    -------
    dict
        corona3d version, the SHA-256 of the configuration used and the configuration itself,
        model run, event time, coordinate frame and units, the observers, the layers, and the
        sources of every input.
    """
    p = rec.config
    text = json.dumps(p, sort_keys=True)
    return dict(
        corona3d_version=__version__, created_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        config_sha256=hashlib.sha256(text.encode()).hexdigest(), config=p,
        event=p["event"], event_time_utc=p["time"], model_run=p["run"],
        frame=FRAME, units=UNITS,
        radio_columns=dict(
            x_rsun_y_rsun="plane-of-sky position seen from Earth, x west and y north [R_sun]",
            x_arcsec_y_arcsec="the same in arcsec, with the angular solar radius at the event time",
            r_rsun="distance from Sun centre at the Newkirk height of the frequency [R_sun]",
            height_rsun="r_rsun - 1 [R_sun]", x_y_z="Carrington Cartesian [R_sun]",
            lon_carrington_deg_lat_deg="Carrington longitude [0, 360) and latitude [deg]",
            lon_stonyhurst_deg="Stonyhurst longitude [-180, 180) [deg] (Carrington longitude minus Earth's)"),
        line_classes=LINE_CLASSES, layers=[k for k, v in rec.layers.items() if v],
        observers={k: dict(lon_carrington_deg=v["lon_deg"], lat_deg=v["lat_deg"],
                           direction=v["direction"].tolist()) for k, v in rec.observers.items()},
        shell_times_utc={k: v["time"] for k, v in rec.gcs.items()}, inputs=_inputs(rec, home))


def _lines_table(rec):
    """All field lines as points with line id and class: (n, 3), (n,), (n,) of str."""
    pts, ids, cls, i = [], [], [], 0
    for k, lines in rec.lines.items():
        for x in lines:
            pts.append(x), ids.append(np.full(len(x), i)), cls.append(np.full(len(x), k))
            i += 1
    if not pts:
        return np.empty((0, 3)), np.empty(0, int), np.empty(0, str)
    return np.vstack(pts), np.concatenate(ids), np.concatenate(cls)


def _arrays(rec):
    """Every layer as a flat dictionary of arrays (the content of the npz file)."""
    a = {}
    if rec.surface:
        a.update({f"surface_{k}": v for k, v in rec.surface.items()})
    pts, ids, cls = _lines_table(rec)
    if len(pts):
        a.update(lines_xyz=pts, lines_id=ids, lines_class=cls)
    for k, g in rec.gcs.items():
        a[f"gcs_{k}_rings"] = g["rings"]
    for c, v in rec.radio.items():
        a[f"radio_{c}"] = v
    for k, g in rec.whitelight.items():
        a.update({f"{k}_{n}": g[n] for n in ("points", "wl", "dims")})
    for k, o in rec.observers.items():
        a[f"observer_{k}"] = np.r_[o["lon_deg"], o["lat_deg"], o["direction"]]
    return a


def _csv(path, header, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def _write_csv(rec, d):
    if rec.surface:
        s = rec.surface
        _csv(d / "surface_vertices.csv", ("x", "y", "z", "va_kms", "ne_cm3"),
             np.c_[s["vertices"], s["va_kms"], s["ne_cm3"]].tolist())
        _csv(d / "surface_faces.csv", ("i", "j", "k"), s["faces"].tolist())
    pts, ids, cls = _lines_table(rec)
    if len(pts):
        _csv(d / "field_lines.csv", ("line_id", "class", "x", "y", "z"),
             ([i, c, *q] for i, c, q in zip(ids.tolist(), cls.tolist(), pts.tolist())))
    for k, g in rec.gcs.items():
        n, m, _ = g["rings"].shape
        ii, jj = np.meshgrid(range(n), range(m), indexing="ij")
        _csv(d / f"gcs_{k}.csv", ("axis_index", "circle_index", "x", "y", "z"),
             np.c_[ii.ravel(), jj.ravel(), g["rings"].reshape(-1, 3)].tolist())
    if rec.radio:
        _csv(d / "radio_sources.csv", RADIO_COLUMNS, zip(*(rec.radio[c].tolist() for c in RADIO_COLUMNS)))
    for k, g in rec.whitelight.items():
        _csv(d / f"{k}_plane.csv", ("x", "y", "z", "wl"), np.c_[g["points"], g["wl"]].tolist())
    _csv(d / "observers.csv", ("observer", "lon_carrington_deg", "lat_deg", "dir_x", "dir_y", "dir_z"),
         [[k, o["lon_deg"], o["lat_deg"], *o["direction"]] for k, o in rec.observers.items()])


def _write_vtk(rec, d):
    import pyvista as pv
    from .reconstruction import polyline
    if rec.surface:
        s = rec.surface
        m = pv.PolyData(s["vertices"], np.c_[np.full(len(s["faces"]), 3), s["faces"]].ravel())
        m["va_kms"], m["ne_cm3"] = s["va_kms"], s["ne_cm3"]
        m.save(d / "surface.vtp")
    for k, lines in rec.lines.items():
        if lines:
            m = polyline(lines)
            m["line_class"] = np.full(m.n_points, LINE_CLASSES[k])
            m.save(d / f"lines_{k}.vtp")
    for k, g in rec.gcs.items():
        n, m_, _ = g["rings"].shape
        x, y, z = (g["rings"][..., i].reshape(n, m_, 1) for i in range(3))
        pv.StructuredGrid(x, y, z).save(d / f"gcs_{k}.vts")
    if rec.radio and len(rec.radio["x"]):
        m = pv.PolyData(np.c_[rec.radio["x"], rec.radio["y"], rec.radio["z"]])
        for c in ("freq_mhz", "tb_k", "x_rsun", "y_rsun", "r_rsun", "lon_carrington_deg", "lat_deg",
                  "lon_stonyhurst_deg"):
            m[c] = rec.radio[c]
        m.save(d / "radio_sources.vtp")
    for k, g in rec.whitelight.items():
        m = pv.StructuredGrid()
        m.points, m.dimensions = g["points"], tuple(int(i) for i in g["dims"])
        m["wl"] = g["wl"]
        m.save(d / f"{k}_plane.vts")


def export(home=".", out="export", fmt="npz", rec=None, params=None, **layers):
    """Write every enabled layer of a working directory to data files and a JSON sidecar.

    Parameters
    ----------
    home : str or Path
        Working directory (see ``corona3d.api.build``).
    out : str or Path
        Output folder (created); relative paths are relative to ``home``.
    fmt : {"npz", "csv", "vtk"}
        ``npz``: ``layers.npz`` with every array. ``csv``: one table per layer. ``vtk``: PyVista
        files (``surface.vtp``, ``lines_<class>.vtp`` (point array ``line_class``), ``gcs_<name>.vts``, ``radio_sources.vtp``,
        ``<instrument>_plane.vts``). The sidecar ``metadata.json`` is written in every case.
    rec : Reconstruction, optional
        A reconstruction already built; otherwise it is built from ``home``, ``params``, ``layers``.
    params, **layers
        As for ``corona3d.api.build``.

    Returns
    -------
    Path
        The output folder.
    """
    if fmt not in FORMATS:
        raise ValueError(f"format {fmt!r}: choose from {', '.join(FORMATS)}")
    rec = rec or build(home, params, **layers)
    d = Path(out) if Path(out).is_absolute() else Path(home) / out
    d.mkdir(parents=True, exist_ok=True)
    if fmt == "npz":
        np.savez_compressed(d / "layers.npz", **_arrays(rec))
    elif fmt == "csv":
        _write_csv(rec, d)
    else:
        _write_vtk(rec, d)
    (d / "metadata.json").write_text(json.dumps(metadata(rec, home), indent=2, default=str) + "\n",
                                     encoding="utf-8")
    return d
