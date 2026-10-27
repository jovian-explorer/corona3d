"""Configurations written by version 0.1.0 and by the commits before 0.2 still load."""
import json

from corona3d.paths import config, legacy

V010 = {"time": "2013-10-25T15:00", "earth": [357.372, 4.993], "stereo_b": [215.329, -0.984],
        "views": [{"from": "earth", "frame": [1.513, 0.087, 0.158]}, {"from": "north"},
                  {"from": "stereo_b", "frame": [1.455, -0.181, 0.093]}]}
NRH = {"folder": "data/nrh", "start": "2013-10-25T14:59:00", "end": "2013-10-25T15:03:00", "step_s": 60,
       "freqs": [150.9, 432.0], "colormaps": ["spring", "winter"], "level": 0.5, "min_tb": 1e8, "contours": False}


def test_stereo_b_keys_become_stereo(tmp_path):
    (tmp_path / "config.json").write_text(json.dumps(V010))
    p = config(tmp_path)
    assert p["stereo"] == [215.329, -0.984] and p["stereo_sc"] == "B" and "stereo_b" not in p
    assert [v["from"] for v in p["views"]] == ["earth", "north", "stereo"]
    assert p["views"][2]["frame"] == [1.455, -0.181, 0.093]


def test_stereo_a_view_sets_the_spacecraft_and_new_keys_win():
    assert legacy({"views": [{"from": "stereo_a"}]})["stereo_sc"] == "A"
    out = legacy({"stereo": [1, 2], "stereo_b": [3, 4], "stereo_sc": "A"})
    assert out["stereo"] == [1, 2] and out["stereo_sc"] == "A" and "stereo_b" not in out


def test_radio_nrh_becomes_a_radio_images_entry(tmp_path):
    (tmp_path / "config.json").write_text(json.dumps(dict(V010, radio_nrh=NRH)))
    p = config(tmp_path)
    (e,) = p["radio_images"]
    assert e["instrument"] == "nrh" and e["freqs"] == [150.9, 432.0] and e["step_s"] == 60
    assert "colormaps" not in e and p["radio_colormaps"] == ["spring", "winter"] and "radio_nrh" not in p
    assert legacy({"radio_nrh": NRH, "radio_images": []})["radio_images"] == []
