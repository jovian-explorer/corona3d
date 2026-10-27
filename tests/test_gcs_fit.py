"""`corona3d gcs`: the height-time step (the limb view, the selection of the front, and the
extrapolation to the time of the reconstruction) on synthetic numbers with a known answer, and the
keys written into config.json by a mocked fit."""
import json
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pytest

from corona3d import gcs as G
from corona3d.cli import main

FIT_JSON = Path(__file__).resolve().parent.parent / "scripts" / "gcs_fit_2013-10-25.json"
needs_tree = pytest.mark.skipif(not FIT_JSON.exists(), reason="source tree not present")


def test_projection_factor_is_one_at_the_limb_and_zero_along_the_line_of_sight():
    assert G.projection_factor(-53.0, 0.0, [np.cos(np.radians(37)), np.sin(np.radians(37)), 0]) == \
        pytest.approx(1.0, abs=1e-9)                          # 90 deg from the observer
    assert G.projection_factor(0.0, 0.0, [1, 0, 0]) == pytest.approx(0.0, abs=1e-9)


def test_front_profile_and_extrapolation_recover_a_known_rise():
    """Front at 1.9 R_sun at 15:00 rising at 0.09 R_sun/min; noise frames before the eruption
    (small envelopes) and truncated fronts at the edge of a 4 R_sun field must be left out."""
    t0 = datetime(2013, 10, 25, 15, 0)
    rows = [(t0 - timedelta(minutes=10), 2.4, 0.3),          # pre-event noise at the seed
            (t0 - timedelta(minutes=2), 2.1, 0.6)]
    for m in (5, 8, 13, 18):
        rows.append((t0 + timedelta(minutes=m), 1.9 + 0.09 * m, 2.0 + m / 4))
    rows += [(t0 + timedelta(minutes=23), 3.95, 9.0),           # front at the edge: truncated
             (t0 + timedelta(minutes=28), 3.70, 9.5)]           # what remains after it left
    kept = G.front_profile(rows, r_out=4.0)
    assert [r for _, r, _ in kept] == pytest.approx([1.9 + 0.09 * m for m in (5, 8, 13, 18)])
    v, h0, err = G.extrapolate([t for t, _, _ in kept], [r for _, r, _ in kept], t0)
    assert v == pytest.approx(0.09, abs=1e-9) and h0 == pytest.approx(1.9, abs=1e-9)
    assert err == pytest.approx(0.0, abs=1e-6)


def fake_fit(*args, **kw):
    return json.loads(FIT_JSON.read_text())


@needs_tree
def test_gcs_command_writes_the_shell_and_its_provenance_into_config_json(tmp_path, monkeypatch):
    main(["--home", str(tmp_path), "init"])
    monkeypatch.setattr(G, "fit", fake_fit)
    main(["--home", str(tmp_path), "gcs", "--views", "c2,cor1a,cor1b", "--fix-direction"])
    cfg = json.loads((tmp_path / "config.json").read_text())
    assert cfg["gcs"] == dict(height=1.931, kappa=0.396, alpha_deg=23.46, lon_stonyhurst=-53.262,
                              lat_deg=3.146, tilt_deg=-146.654)
    assert cfg["gcs_wl"]["height"] == 5.538 and cfg["gcs_wl"]["tilt_deg"] == -146.654
    rec = cfg["gcs_fit"]
    assert rec["views"] == ["c2", "cor1a", "cor1b"] and rec["iou"] == dict(c2=0.761, cor1a=0.854, cor1b=0.744)
    assert rec["speed_km_s"] == 1039.0 and rec["height_t_out_err"] == 0.158
    assert rec["fixed"] == ["lon_stonyhurst", "lat_deg"] and rec["flare"] == [-69.0, -6.0]
    assert (tmp_path / "gcs_fit.json").exists()


@needs_tree
def test_the_fit_file_of_2013_10_25_uses_the_renamed_keys():
    d = json.loads(FIT_JSON.read_text())
    assert {"t_ref", "t_out"} <= set(d) and not {"t1524", "t1500"} & set(d)
    assert d["t_out"]["height"] == d["fit"]["height_t_out"] == 1.931


@needs_tree
def test_assumed_height_is_recorded_as_assumed_with_its_note(tmp_path, monkeypatch):
    main(["--home", str(tmp_path), "init"])
    seen = {}

    def fit(*args):
        seen["args"] = args
        d = json.loads(FIT_JSON.read_text())
        d["t_out"]["height"] = args[-2]
        d["fit"] = dict(height_t_out=args[-2], assumed=True, note=args[-1])
        return d
    monkeypatch.setattr(G, "fit", fit)
    with pytest.raises(SystemExit):                              # a height without a reason
        main(["--home", str(tmp_path), "gcs", "--height-out", "1.77"])
    main(["--home", str(tmp_path), "gcs", "--height-out", "1.77", "--height-note", "kinematic estimate",
          "--running-diff", "10"])
    cfg = json.loads((tmp_path / "config.json").read_text())
    assert cfg["gcs"]["height"] == 1.77
    assert cfg["gcs_fit"]["height_source"] == "assumed" and cfg["gcs_fit"]["height_note"] == "kinematic estimate"
    assert cfg["gcs_fit"]["running_diff_min"] == 10.0 and "speed_km_s" not in cfg["gcs_fit"]
    assert seen["args"][-3] == 10.0


def test_running_difference_takes_the_frame_nearest_to_t_minus_the_interval(monkeypatch):
    t0 = datetime(2024, 5, 9, 17, 0)
    fr = {"cor1a": [(t0 + timedelta(minutes=5 * i), [Path(f"x_s4c1A_{i}")]) for i in range(6)]}
    monkeypatch.setattr(G, "total_b", lambda files: files[0].name)
    monkeypatch.setattr(G, "seq", lambda f: "s4")
    t = t0 + timedelta(minutes=20)
    assert G.base_of(fr, {("cor1a", "s4"): "first"}, "cor1a", t, fr["cor1a"][4][1], None) == "first"
    assert G.base_of(fr, {}, "cor1a", t, fr["cor1a"][4][1], 10) == "x_s4c1A_2"
    assert G.base_of(fr, {}, "cor1a", t0, fr["cor1a"][0][1], 10) is None
