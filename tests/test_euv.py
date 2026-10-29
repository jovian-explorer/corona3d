"""EUV images: the recognition of AIA and EUVI files, the height of best agreement, the apex height at
t_out from a series of EUV frames, and the record of an EUV height in the fit provenance and the report.
The default fit, without --euv, does not read EUV images at all."""
from datetime import datetime, timedelta
from types import SimpleNamespace

import numpy as np
import pytest
from astropy.io import fits

from corona3d import euv, gcs


def test_aia_and_euvi_headers_are_recognised():
    assert euv.key_of({"TELESCOP": "SDO/AIA", "WAVELNTH": 193}) == "aia193"
    assert euv.key_of({"INSTRUME": "AIA_3", "WAVELNTH": 211}) == "aia211"
    assert euv.key_of({"TELESCOP": "SDO/AIA", "WAVELNTH": 171}) is None
    assert euv.key_of({"DETECTOR": "EUVI", "WAVELNTH": 195, "OBSRVTRY": "STEREO_B"}) == "euvi195b"
    assert euv.key_of({"DETECTOR": "EUVI", "WAVELNTH": 195, "OBSRVTRY": "STEREO_A"}) == "euvi195a"
    assert euv.key_of({"DETECTOR": "COR1", "WAVELNTH": 0}) is None


def test_frames_are_grouped_by_view_and_time(tmp_path):
    for name, hd in (("20131025_150122_n4euB.fts", dict(DETECTOR="EUVI", WAVELNTH=195, OBSRVTRY="STEREO_B",
                                                          **{"DATE-OBS": "2013-10-25T15:01:22.000"})),
                     ("20131025_145622_n4euB.fts", dict(DETECTOR="EUVI", WAVELNTH=195, OBSRVTRY="STEREO_B",
                                                          **{"DATE-OBS": "2013-10-25T14:56:22.000"})),
                     ("20131025_145822_n4c1B.fts", dict(DETECTOR="COR1", WAVELNTH=0, **{"DATE-OBS": "2013-10-25T14:58:22"}))):
        fits.PrimaryHDU(np.zeros((4, 4)), header=fits.Header(hd)).writeto(tmp_path / name)
    fr = euv.frames(tmp_path)
    assert list(fr) == ["euvi195b"]
    assert [t for t, _ in fr["euvi195b"]] == [datetime(2013, 10, 25, 14, 56, 22), datetime(2013, 10, 25, 15, 1, 22)]


def test_best_height_and_its_span():
    h, s, (lo, hi) = euv.best_height(lambda h: 0.8 - abs(h - 1.4), tol=0.02)
    assert h == pytest.approx(1.4) and s == pytest.approx(0.8) and (lo, hi) == pytest.approx((1.38, 1.42))


def _profile(monkeypatch, tmp_path, true, minutes, speed=None):
    """height_profile on synthetic EUVI views whose IoU peaks at true(t)."""
    t0 = datetime(2013, 10, 25, 15, 0)
    times = [t0 + timedelta(minutes=m) for m in minutes]
    monkeypatch.setattr(euv, "frames", lambda folder: {"euvi195b": [(t, ["f"]) for t in times]})
    monkeypatch.setattr(gcs, "_load", lambda: None)
    monkeypatch.setattr(gcs, "bases_of", lambda fr, cut: ({("euvi195b", "n4"): SimpleNamespace(date=SimpleNamespace(isot="2013-10-25T14:45:00"))}, {}))
    monkeypatch.setattr(gcs, "shell", lambda q, t: (q[0], datetime.fromisoformat(t)))

    class View:
        def iou(self, shell):
            h, t = shell
            return 1.0 - abs(h - true(t))
    monkeypatch.setattr(gcs, "views_at", lambda fr, base, t, keys, flare: {"euvi195b": (t, View())})
    return t0, times, euv.height_profile(tmp_path, [0.4, 25, -60, 0, -30], t0, (-69, -6), speed_km_s=speed)


def test_height_profile_fits_a_line_through_the_euv_heights(monkeypatch, tmp_path):
    """Synthetic views whose IoU peaks at the true height of a shell rising at 300 km/s."""
    t0 = datetime(2013, 10, 25, 15, 0)
    true = lambda t: 1.3 + 300.0 / 695700 * (t - t0).total_seconds()
    _, times, rec = _profile(monkeypatch, tmp_path, true, (-5, 0, 5))
    assert [r["height"] for r in rec["frames"]] == pytest.approx([true(t) for t in times], abs=0.006)
    assert rec["method"] == "line" and rec["height_t_out"] == pytest.approx(1.3, abs=0.005)
    assert rec["speed_km_s"] == pytest.approx(300.0, rel=0.05)


def test_an_apex_beyond_the_euv_field_is_not_used_and_one_frame_is_moved_to_t_out(monkeypatch, tmp_path):
    t0 = datetime(2013, 10, 25, 15, 0)
    true = lambda t: 1.5 + 1000.0 / 695700 * (t - t0).total_seconds()     # 1.5 at 15:00, 1.93 at 15:05
    _, times, rec = _profile(monkeypatch, tmp_path, true, (1, 5), speed=1000.0)
    assert [r["used"] for r in rec["frames"]] == [True, False]           # 1.93 R_sun lies beyond 1.7 R_sun
    assert "beyond the outer edge" in rec["frames"][1]["rejected"]
    assert rec["height_t_out"] == pytest.approx(1.5, abs=0.006)          # 15:01 frame moved back by 60 s
    assert "moved by -60 s" in rec["method"]


def test_an_euv_height_is_recorded_with_the_extrapolated_one():
    out = {"t_ref": {"height": 5.5, "kappa": .4, "alpha_deg": 23, "lon_stonyhurst": -53, "lat_deg": 3,
                     "tilt_pythea_deg": -56, "time": "2013-10-25T15:24:30", "iou": {"c2": .8}},
           "t_out": {"height": 1.4, "kappa": .4, "alpha_deg": 23, "lon_stonyhurst": -53, "lat_deg": 3,
                     "tilt_pythea_deg": -56, "time": "2013-10-25T15:00:00"},
           "height_time": {"instrument": "cor1b", "projection": 1.0, "frames": []},
           "fit": {"speed_rsun_per_min": 0.09, "speed_km_s": 1039.0, "height_t_out": 1.4, "height_t_out_err": 0.05,
                   "n_frames": 3, "height_from": "euv", "height_t_out_wl": 1.931},
           "euv": {"keys": ["aia193", "euvi195b"], "method": "line", "speed_km_s": 700.0}}
    rec = gcs.config_update(out, (-69, -6), ("c2",))["gcs_fit"]
    assert rec["height_from"] == "euv" and rec["height_t_out"] == 1.4 and rec["height_t_out_wl"] == 1.931
    assert rec["euv_views"] == ["aia193", "euvi195b"] and rec["speed_km_s"] == 1039.0


def test_euv_folder_is_only_read_with_the_option(tmp_path):
    assert gcs.euv_folder(tmp_path, SimpleNamespace()) is None
    with pytest.raises(SystemExit, match="corona3d fetch --euv"):
        gcs.euv_folder(tmp_path, SimpleNamespace(euv=True))
