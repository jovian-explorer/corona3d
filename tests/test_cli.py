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
    for f in ("test.png", "view_earth.png", "view_north.png", "view_stb.png"):
        assert (tmp_path / "figures" / f).stat().st_size > 10_000
    assert set(LAYERS) == {"surface", "sphere", "open", "closed", "streamer", "gcs", "radio",
                           "c2", "cor1", "cor2"}
