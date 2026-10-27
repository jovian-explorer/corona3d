"""The working directory, the layer switches, and a render without any downloaded data."""
import json

import pytest

from corona3d import LAYERS, make_figure
from corona3d.cli import main


def test_render_outside_a_working_directory_stops_with_advice(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("CORONA3D_HOME", raising=False)
    with pytest.raises(SystemExit) as e:
        main(["render"])
    assert "corona3d init" in str(e.value)


def test_init_writes_the_event_and_keeps_an_edited_config(tmp_path):
    main(["--home", str(tmp_path), "init"])
    cfg = json.loads((tmp_path / "config.json").read_text())
    assert cfg["time"] == "2013-10-25T15:00" and cfg["ne_level"] == 9e7
    (tmp_path / "config.json").write_text(json.dumps({"density_scale": 2.5}))
    main(["--home", str(tmp_path), "init"])
    assert json.loads((tmp_path / "config.json").read_text()) == {"density_scale": 2.5}


def test_unknown_layer_is_refused(tmp_path):
    with pytest.raises(ValueError, match="unknown layers"):
        make_figure(tmp_path / "x.png", tmp_path, flares=True)


def test_render_of_the_layers_that_need_no_data(tmp_path):
    """Sphere, GCS shell, and two radio sources: no model or image is read."""
    main(["--home", str(tmp_path), "init"])
    off = [f"--no-{k}" for k in ("surface", "open", "closed", "streamer")]
    src = json.dumps([{"freq": 173, "x": [-1.19, -1.09], "y": [-0.15, -0.08]},
                      {"freq": 228, "x": -1100, "y": -190, "unit": "arcsec"}])
    main(["--home", str(tmp_path), "render", *off, "--radio", "--radio-sources", src,
          "--out", "figures/test.png"])
    for f in ("test.png", "view_a.png", "view_b.png", "view_c.png"):
        assert (tmp_path / "figures" / f).stat().st_size > 10_000
    assert set(LAYERS) == {"surface", "sphere", "open", "closed", "streamer", "gcs", "radio",
                           "c2", "cor1", "cor2"}


def test_partial_nested_keys_are_merged(tmp_path):
    from corona3d.paths import config
    main(["--home", str(tmp_path), "init"])
    (tmp_path / "config.json").write_text(json.dumps({"gcs": {"height": 2.0}}))
    g = config(tmp_path)["gcs"]
    assert g["height"] == 2.0                                      # the other keys are kept
    assert set(g) >= {"height", "kappa", "alpha_deg", "lon_stonyhurst", "lat_deg", "tilt_deg"}


def _disc_radius(png):
    """Radius in pixels of the grey sphere, from the width of its row through the centre."""
    from PIL import Image
    import numpy as np
    a = np.asarray(Image.open(png).convert("L"), float)
    row = a[a.shape[0] // 2]
    return (row < 250).sum() / 2


def test_half_sets_the_scale_of_the_panels(tmp_path):
    """Orthographic panels: halving the half-height doubles the radius of the sphere on screen."""
    main(["--home", str(tmp_path), "init"])
    off = [f"--no-{k}" for k in ("surface", "open", "closed", "streamer", "gcs")]
    r = {}
    for h in (3.0, 1.5):
        main(["--home", str(tmp_path), "render", *off, "--half", str(h), "--out", f"h{h}/f.png"])
        r[h] = _disc_radius(tmp_path / f"h{h}" / "view_b.png")
    assert r[1.5] / r[3.0] == pytest.approx(2.0, rel=0.03)
    with pytest.raises(ValueError, match="positive"):
        make_figure("x.png", tmp_path, half=0, surface=False, open=False, closed=False, streamer=False)


def test_model_fetch_skips_present_files_and_cleans_failures(tmp_path, monkeypatch):
    import io
    import urllib.request
    from corona3d import fetch
    from corona3d.paths import config, model_dir
    main(["--home", str(tmp_path), "init"])
    cfg = config(tmp_path)
    d = model_dir(tmp_path, cfg)
    d.mkdir(parents=True)
    (d / "rho002.hdf").write_bytes(b"present")
    calls = []

    def fake(url, timeout):
        calls.append(url)
        if "bp002" in url:
            raise OSError("connection reset")
        return io.BytesIO(b"data")
    monkeypatch.setattr(urllib.request, "urlopen", fake)
    assert fetch.model(tmp_path, cfg) == 1                        # bp002 failed
    assert not any("rho002" in u for u in calls)                  # present file skipped
    assert (d / "br002.hdf").read_bytes() == b"data"
    assert not (d / "bp002.hdf").exists() and not list(d.glob("*.part"))
    assert calls[0] == "https://www.predsci.com/data/runs/cr2143-medium/hmi_mast_mas_std_0201/corona/br002.hdf"


def test_damaged_images_are_deleted_and_counted(tmp_path, monkeypatch):
    import numpy as np
    from astropy.io import fits
    from corona3d import fetch
    from corona3d.paths import config
    main(["--home", str(tmp_path), "init"])
    wl = tmp_path / "data" / "wl"
    wl.mkdir(parents=True)
    fits.PrimaryHDU(np.zeros((64, 64), np.float32)).writeto(wl / "good.fts")
    raw = (wl / "good.fts").read_bytes()
    (wl / "cut.fts").write_bytes(raw[:4000])                     # a truncated download

    class Res(list):
        errors = []
    import sunpy.net
    monkeypatch.setattr(sunpy.net.Fido, "search", lambda *a, **k: [])
    monkeypatch.setattr(sunpy.net.Fido, "fetch", lambda *a, **k: Res())
    assert fetch.whitelight(tmp_path, config(tmp_path)) == 1
    assert (wl / "good.fts").exists() and not (wl / "cut.fts").exists()


def test_views_can_be_set_for_every_panel():
    import numpy as np
    from corona3d.geometry import unit_vec
    from corona3d.paths import default_config
    from corona3d.reconstruction import parse_view, view_geometry
    p = default_config()
    assert parse_view("b=north") == (1, {"from": "north"})
    assert parse_view("c=300,20,45") == (2, {"lon": 300.0, "lat": 20.0, "roll": 45.0})
    d, up, owner = view_geometry({"lon": 300.0, "lat": 0.0}, p)
    assert np.allclose(d, unit_vec(300, 0)) and np.allclose(up, [0, 0, 1]) and owner is None
    d, up, _ = view_geometry({"lon": 300.0, "lat": 0.0, "roll": 90.0}, p)
    assert np.allclose(up, np.cross(d, [0, 0, 1]))            # turned by 90 deg about the view axis
    d, up, owner = view_geometry({"from": "stereo_b"}, p)
    assert owner == "stereo_b" and np.allclose(d, unit_vec(*p["stereo_b"]))
    d, up, _ = view_geometry({"from": "north"}, p)
    assert np.allclose(d, [0, 0, 1], atol=1e-9) and up @ unit_vec(*p["earth"]) < -0.99   # Earth below
    with pytest.raises(ValueError):
        view_geometry({"from": "mars"}, p)


def test_perspective_projection_magnifies_the_near_side(tmp_path):
    """The camera of the published figure (perspective, 6.65 R_sun) draws the sphere larger than
    the orthographic projection at the same half-height, since its near side is closer."""
    main(["--home", str(tmp_path), "init"])
    off = dict(surface=False, open=False, closed=False, streamer=False, gcs=False)
    r = {}
    for proj in ("orthographic", "perspective"):
        make_figure(f"{proj}/f.png", tmp_path, params={"projection": proj}, **off)
        r[proj] = _disc_radius(tmp_path / proj / "view_b.png")
    assert 1.005 < r["perspective"] / r["orthographic"] < 1.4      # the limb lies slightly nearer
