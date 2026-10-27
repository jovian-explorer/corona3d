"""The Python API (layers as arrays) and `corona3d export`, on the synthetic run (uniform
density, dipole field) and synthetic NRH files."""
import json

import numpy as np
import pytest

import corona3d
from corona3d.api import carrington_angles
from corona3d.cli import main
from corona3d.geometry import unit_vec
from corona3d.reconstruction import newkirk_r
from test_nrh import T0, write_nrh

RADIO = [{"freq": 173, "x": [-1.19, -1.09], "y": [-0.15, -0.08]},
         {"freq": 228, "x": -1100, "y": -190, "unit": "arcsec"}]
LOW = dict(ne_level=1e8, n_seeds=40)                       # a few lines keep the tests quick


@pytest.fixture
def rec(synthetic_home):
    from conftest import write_mas
    from corona3d.paths import config, model_dir
    run = model_dir(synthetic_home, config(synthetic_home))
    for f in run.glob("*.hdf"):
        f.unlink()
    write_mas(run, n_model=2.0, radial=True)                 # n_e = 4e8 / R^2
    d = synthetic_home / "data" / "nrh"
    d.mkdir()
    write_nrh(d / "nrh2_1509_h60_20131025_145800c00_i.fts", 150.9, T0, 40,
              lambda k: (-1.0 - 0.01 * k, -0.1 - 0.01 * k))
    entry = dict(instrument="nrh", folder="data/nrh", freqs=[150.9], start="2013-10-25T14:59:00",
                 end="2013-10-25T15:03:00", step_s=60, level=0.5, min_tb=1e8, contours=False)
    return corona3d.build(synthetic_home, dict(LOW, radio_sources=RADIO, radio_images=[entry]), radio=True)


def test_top_level_names():
    for k in ("load_config", "build", "export", "Reconstruction", "reconstruct", "make_figure"):
        assert k in corona3d.__all__ and getattr(corona3d, k)
    assert callable(corona3d.export)                       # the function, not the module


def test_unknown_layer_is_refused(synthetic_home):
    with pytest.raises(ValueError, match="unknown layers"):
        corona3d.build(synthetic_home, flares=True)


def test_surface_lines_and_shell_have_the_documented_shapes(rec):
    s = rec.surface
    n = len(s["vertices"])
    assert s["vertices"].shape == (n, 3) and s["faces"].shape[1] == 3 and s["faces"].max() < n
    assert s["va_kms"].shape == s["ne_cm3"].shape == (n,)
    assert np.allclose(s["ne_cm3"], 1e8, rtol=1e-2)                    # the iso-level of ne_level
    assert set(rec.lines) == {"open", "closed", "streamer"}
    assert rec.lines["closed"] and rec.lines["streamer"]               # a dipole closes the low lines
    assert all(x.ndim == 2 and x.shape[1] == 3 for k in rec.lines for x in rec.lines[k])
    g = rec.gcs["time"]["rings"]
    assert g.ndim == 3 and g.shape[2] == 3 and set(rec.gcs) == {"time", "image"}
    assert rec.gcs["image"]["time"] == "2013-10-25T15:24"
    assert np.linalg.norm(rec.gcs["image"]["rings"], axis=2).max() > np.linalg.norm(g, axis=2).max()
    assert rec.observers["earth"]["lon_deg"] == pytest.approx(357.372)
    assert np.allclose(rec.observers["earth"]["direction"], unit_vec(357.372, 4.993))


def test_radio_table_positions_and_frames(rec):
    r = rec.radio
    n = len(r["freq_mhz"])
    assert n == 3 + 5 and all(len(v) == n for v in r.values())          # 3 user positions, 5 NRH blobs
    assert list(r["instrument"][:3]) == ["user"] * 3 and set(r["instrument"][3:]) == {"nrh"}
    assert np.isnan(r["tb_k"][:3]).all() and (r["tb_k"][3:] > 1e8).all()
    assert r["time_utc"][3].startswith("2013-10-25T14:5") and r["time_utc"][0] == ""
    q = np.c_[r["x"], r["y"], r["z"]]
    assert np.allclose(r["r_rsun"], np.linalg.norm(q, axis=1))
    assert np.allclose(r["r_rsun"][[0, 1]], newkirk_r(173)) and np.allclose(r["r_rsun"][2], newkirk_r(228))
    assert np.allclose(r["height_rsun"], r["r_rsun"] - 1)
    lon, lat, _ = carrington_angles(q)
    assert np.allclose(lon, r["lon_carrington_deg"]) and np.allclose(lat, r["lat_deg"])
    assert np.allclose(r["lon_stonyhurst_deg"], (lon - 357.372 + 180) % 360 - 180)
    assert r["x_arcsec"][2] == pytest.approx(-1100, rel=1e-6) and r["y_arcsec"][2] == pytest.approx(-190, rel=1e-6)
    assert -1.2 < r["x_rsun"][2] < -1.1                                  # about 960 arcsec per R_sun
    assert (r["lon_stonyhurst_deg"] < 0).all() and (np.abs(r["lon_stonyhurst_deg"]) < 90).all()


def test_export_npz_sidecar_and_round_trip(rec, synthetic_home):
    d = corona3d.export(synthetic_home, "out_npz", "npz", rec=rec)
    z = np.load(d / "layers.npz")
    assert np.array_equal(z["surface_vertices"], rec.surface["vertices"])
    assert np.array_equal(z["surface_va_kms"], rec.surface["va_kms"])
    assert len(np.unique(z["lines_id"])) == sum(len(v) for v in rec.lines.values())
    assert set(np.unique(z["lines_class"])) == {k for k, v in rec.lines.items() if v}
    assert np.array_equal(z["gcs_time_rings"], rec.gcs["time"]["rings"])
    assert np.array_equal(z["radio_lon_stonyhurst_deg"], rec.radio["lon_stonyhurst_deg"])
    m = json.loads((d / "metadata.json").read_text())
    assert m["corona3d_version"] == corona3d.__version__ and len(m["config_sha256"]) == 64
    assert m["event_time_utc"] == "2013-10-25T15:00" and m["model_run"] == rec.config["run"]
    assert "Carrington" in m["frame"] and m["units"]["position"] == "R_sun"
    assert m["config"]["density_scale"] == 2.0 and m["observers"]["earth"]["lat_deg"] == 4.993
    assert {"mhd_model", "radio", "gcs"} <= set(m["inputs"])
    assert m["inputs"]["radio"]["image_sets"][0]["files"] == ["nrh2_1509_h60_20131025_145800c00_i.fts"]


def test_export_csv_round_trip(rec, synthetic_home):
    d = corona3d.export(synthetic_home, "out_csv", "csv", rec=rec)
    v = np.loadtxt(d / "surface_vertices.csv", delimiter=",", skiprows=1)
    assert np.allclose(v[:, :3], rec.surface["vertices"]) and np.allclose(v[:, 3], rec.surface["va_kms"])
    f = np.loadtxt(d / "surface_faces.csv", delimiter=",", skiprows=1, dtype=int)
    assert np.array_equal(f, rec.surface["faces"])
    first = (d / "field_lines.csv").read_text().splitlines()
    assert first[0] == "line_id,class,x,y,z" and first[1].split(",")[1] in ("open", "closed", "streamer")
    rows = (d / "radio_sources.csv").read_text().splitlines()
    assert len(rows) == 1 + len(rec.radio["freq_mhz"]) and "lon_stonyhurst_deg" in rows[0]


def test_export_vtk_round_trip(rec, synthetic_home):
    import pyvista as pv
    d = corona3d.export(synthetic_home, "out_vtk", "vtk", rec=rec)
    m = pv.read(d / "surface.vtp")
    assert np.allclose(m.points, rec.surface["vertices"], atol=1e-5)
    assert np.allclose(m["va_kms"], rec.surface["va_kms"], rtol=1e-5)
    assert m.n_cells == len(rec.surface["faces"])
    c = pv.read(d / "lines_closed.vtp")
    assert c.n_lines == len(rec.lines["closed"]) and set(c["line_class"]) == {1}
    assert pv.read(d / "gcs_time.vts").dimensions[:2] == rec.gcs["time"]["rings"].shape[:2]
    assert pv.read(d / "radio_sources.vtp").n_points == len(rec.radio["freq_mhz"])


def test_cli_export_writes_enabled_layers_only(synthetic_home):
    main(["--home", str(synthetic_home), "export", "--no-surface", "--no-open", "--no-closed",
          "--no-streamer", "--format", "csv", "--out", "e"])
    names = {f.name for f in (synthetic_home / "e").iterdir()}
    assert names == {"gcs_time.csv", "gcs_image.csv", "observers.csv", "metadata.json"}
    with pytest.raises(ValueError, match="format"):
        corona3d.export(synthetic_home, "x", "hdf5")


def test_gcs_note_does_not_break_the_config(synthetic_home):
    from corona3d.paths import config, event
    assert event("2024-05-09")["gcs_note"]
    (synthetic_home / "config.json").write_text(json.dumps(dict(config(synthetic_home), gcs_note="assumed height")))
    assert corona3d.load_config(synthetic_home)["gcs_note"] == "assumed height"
    rec = corona3d.build(synthetic_home, surface=False, open=False, closed=False, streamer=False)
    assert rec.gcs["time"]
