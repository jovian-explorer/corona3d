"""Small offline checks: a pending `density_scale: "pb"`, bursts that cross midnight, the colour-map rule,
the C2 axes from the header, the OVRO-LWA position check that fails, and a field line that touches a surface."""
import json

import numpy as np
import pytest

from corona3d import bursts, paths, pb, sources
from corona3d import shockphys as sp


def test_pb_scale_without_a_fit_does_not_stop_the_commands(tmp_path):
    (tmp_path / "config.json").write_text(json.dumps({"density_scale": "pb"}))
    with pytest.warns(UserWarning, match="corona3d pb"):
        cfg = paths.config(tmp_path)
    assert cfg["density_scale"] == paths.default_config()["density_scale"]
    with pytest.raises(ValueError, match="corona3d pb"):               # the places that use the scale stay strict
        pb.resolve(cfg)
    with pytest.raises(ValueError, match="corona3d pb"):
        sources.analyse(tmp_path)
    (tmp_path / "config.json").write_text(json.dumps({"density_scale": "pb", "pb_fit": {"scale": 3.5}}))
    assert paths.config(tmp_path)["density_scale"] == 3.5


def test_burst_that_crosses_midnight_ends_the_next_day():
    text = " 700       2350   ////      0010  CUL  C   RSP  023-090   II/1   900\n"
    (r,) = bursts.parse(text, "2013-10-25")
    assert r["begin"].day == 25 and r["end"].day == 26 and (r["end"] - r["begin"]).total_seconds() == 1200


def test_colormaps_of_the_user_replace_the_default_colours():
    assert paths.colormap_rule({"radio_colormaps": ["hot", "cool"]})["radio_colors"] is None
    assert "radio_colors" not in paths.colormap_rule({"radio_colormaps": ["spring", "winter"]})
    assert paths.colormap_rule({"radio_colormaps": ["hot"], "radio_colors": ["red"]})["radio_colors"] == ["red"]


def test_c2_axes_centre_on_the_sun_and_have_north_up():
    from astropy.io import fits
    h = fits.Header()
    h.update({"R1COL": 20, "R2COL": 1043, "CRPIX1": 512.5, "CRPIX2": 512.5, "CROTA1": 0.0, "DATE-OBS": "2013/10/25",
              "TIME-OBS": "15:12:05.000", "TELESCOP": "SOHO", "INSTRUME": "LASCO", "DETECTOR": "C2", "CTYPE1": "SOLAR-X",
              "CTYPE2": "SOLAR-Y", "CDELT1": 11.9, "CDELT2": 11.9, "CRVAL1": 0, "CRVAL2": 0, "CUNIT1": "arcsec", "CUNIT2": "arcsec"})
    x, y = pb.c2_axes(h, (1024, 1024))
    assert abs(x[511, 511]) < 0.01 and abs(y[511, 511]) < 0.01
    assert y[900, 511] > y[100, 511] and x[511, 900] > x[511, 100]


def test_failed_ovro_check_is_recorded_not_raised(tmp_path, monkeypatch):
    (tmp_path / "config.json").write_text("{}")
    monkeypatch.setattr(sources, "ovro_check", lambda *a, **k: (_ for _ in ()).throw(OSError("offline")))
    rows, summ = sources.analyse(tmp_path, ovro_offset=True, log=lambda *a: None)
    assert summ["ovro_lwa_offset"] == {"failed": "offline"}


def test_crossing_from_a_point_on_the_surface_counts():
    class S:
        def nearest(self, line):
            return np.zeros(len(line)), np.array([0.0, -1.0, 0.5, 1.0]), np.arange(len(line))
    c = sp.crossings(np.arange(12.0).reshape(4, 3), S(), 1.0)
    assert [x["j"] for x in c] == [0, 1] and not c[0]["outward"] and c[1]["outward"]
