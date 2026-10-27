"""Synthetic coronagraph frames: the plane must put a feature at its true place in three
dimensions at any roll angle, and the frame selection must refuse inputs that would draw a
wrong or blank image."""
import astropy.units as u
import numpy as np
import pytest
import sunpy.map
from astropy.coordinates import SkyCoord
from astropy.io import fits
from sunpy.coordinates import frames

from corona3d.geometry import observer_basis
from corona3d.whitelight import _frames, _image, load_plane

OBS = (215.329, -0.984)                  # STEREO-B, Carrington lon/lat on 2013-10-25
WIN = ["2013-10-25T14:35", "2013-10-25T15:40"]


def _write(d, name, time, blob=None, roll=0.0, n=512, scale=30.0, level=1.0):
    obs = SkyCoord(-142.0 * u.deg, OBS[1] * u.deg, 1 * u.AU, frame=frames.HeliographicStonyhurst,
                   obstime=time)
    ref = SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame=frames.Helioprojective(observer=obs, obstime=time))
    hdr = sunpy.map.make_fitswcs_header((n, n), ref, scale=[scale, scale] * u.arcsec / u.pix,
                                        rotation_angle=roll * u.deg)
    hdr["exptime"] = 6.0
    m = sunpy.map.Map(np.zeros((n, n)), hdr)
    data = np.full((n, n), level)
    if blob is not None:                 # a 4-pixel square at (x, y) R_sun in solar north-up axes
        rs = m.rsun_obs.to_value(u.arcsec)
        c = SkyCoord(blob[0] * rs * u.arcsec, blob[1] * rs * u.arcsec, frame=m.coordinate_frame)
        i, j = (int(round(float(v))) for v in m.wcs.world_to_pixel(c))
        k = max(2, n // 256)
        data[j - k:j + k, i - k:i + k] += 100.0
    sunpy.map.Map(data, m.meta).save(str(d / name), overwrite=True)


def _pair(d, roll=0.0, n_base=512, n=512):
    _write(d, "20131025_143900_d4c2B.fts", "2013-10-25T14:39:51", roll=roll, n=n_base, scale=30.0 * 512 / n_base)
    _write(d, "20131025_152400_d4c2B.fts", "2013-10-25T15:24:51", blob=(5.0, 2.0), roll=roll, n=n,
           scale=30.0 * 512 / n)


@pytest.mark.parametrize("roll", [0.0, 6.6, 90.0, 175.3, 271.0])
@pytest.mark.parametrize("mode", ["ratio", "difference"])
def test_plane_places_a_feature_in_3d(tmp_path, roll, mode):
    _pair(tmp_path, roll)
    g, tf, tb = load_plane("cor2", tmp_path, "2013-10-25T15:24", OBS, WIN, mode)
    v = np.asarray(g["wl"])
    p = g.points[np.nanargmax(v)]
    ex, ey, ez = observer_basis(*OBS)
    assert p @ ex == pytest.approx(5.0, abs=0.15)          # west
    assert p @ ey == pytest.approx(2.0, abs=0.15)          # north
    assert abs(p @ ez) < 1e-9                              # in the plane of the sky
    r = np.hypot(g.points @ ex, g.points @ ey)
    assert np.isnan(v[r < 3.9]).all()                      # COR2 occulter left transparent
    assert (tf.hour, tf.minute, tb.hour, tb.minute) == (15, 24, 14, 39)


def test_frames_of_different_size_are_compared_on_one_grid(tmp_path):
    _pair(tmp_path, n_base=1024, n=2048)                   # the COR2 base of 2013-10-25 is 1024 px
    g, _, _ = load_plane("cor2", tmp_path, "2013-10-25T15:24", OBS, WIN)
    p = g.points[np.nanargmax(np.asarray(g["wl"]))]
    assert p @ observer_basis(*OBS)[0] == pytest.approx(5.0, abs=0.15)


def test_small_images_are_not_binned_to_nothing(tmp_path):
    _write(tmp_path, "x.fts", "2013-10-25T15:00:00", n=128, scale=120.0)
    assert _image([tmp_path / "x.fts"]).data.shape == (128, 128)


def test_inputs_that_would_draw_a_wrong_or_blank_image_are_refused(tmp_path):
    _pair(tmp_path)
    with pytest.raises(ValueError, match="base frame"):     # nearest frame is the base itself
        load_plane("cor2", tmp_path, "2013-10-25T14:40", OBS, WIN)
    with pytest.raises(ValueError, match="outside wl_window"):
        load_plane("cor2", tmp_path, "2013-10-25T16:30", OBS, WIN)
    with pytest.raises(FileNotFoundError, match="no cor2 frames"):
        load_plane("cor2", tmp_path, "2013-10-26T00:10", OBS, ["2013-10-26T00:00", "2013-10-26T01:00"])
    g, tf, _ = load_plane("cor2", tmp_path, "2013-10-25T15:24Z", OBS, WIN)   # time zone accepted
    assert tf.minute == 24


def test_cor1_triplets_are_grouped_by_time_not_by_file_name(tmp_path):
    for s in ("144500", "144512", "144524",          # one triplet
              "145550", "145602", "145614",          # a triplet across a minute boundary
              "152500", "152512"):                   # an incomplete triplet
        h = fits.Header({"DATE-OBS": f"2013-10-25T{s[:2]}:{s[2:4]}:{s[4:]}.000"})
        fits.PrimaryHDU(np.zeros((2, 2)), h).writeto(tmp_path / f"20131025_{s}_s4c1B.fts")
    fr = _frames(tmp_path, "cor1")
    assert [f[0].strftime("%H:%M:%S") for f in fr] == ["14:45:00", "14:55:50"]
    assert all(len(f[1]) == 3 for f in fr)


def test_c2_frames_are_selected_from_their_headers(tmp_path):
    def c2(name, nx=1024, **kw):                           # LASCO level 0.5: DATE-OBS and TIME-OBS
        h = fits.Header({"DATE-OBS": "2013/10/25", "TIME-OBS": kw.pop("t"), **kw})
        fits.PrimaryHDU(np.zeros((1, nx), np.float32), h).writeto(tmp_path / name)
    c2("22478907.fts", t="14:48:05.407", DETECTOR="C2", POLAR="Clear")
    c2("22478908.fts", t="14:58:00.000", DETECTOR="C2", POLAR="+60 Deg")
    c2("22478910.fts", t="14:59:00.000", DETECTOR="C2", POLAR="Clear", nx=512)
    c2("22478909.fts", t="15:00:00.000")                   # no DETECTOR or POLAR keyword
    fr = _frames(tmp_path, "c2")
    assert [f[1][0].name for f in fr] == ["22478907.fts"]
    assert fr[0][0].strftime("%H:%M:%S") == "14:48:05"


def test_cor1_triplets_of_the_normal_sequence_are_read(tmp_path):
    """STEREO-A COR1 in 2024 writes its triplets as n4c1A (0, 120, 240 deg, 18 s apart)."""
    for s in ("144100", "144118", "144136"):
        h = fits.Header({"DATE-OBS": f"2024-05-09T{s[:2]}:{s[2:4]}:{s[4:]}.006"})
        fits.PrimaryHDU(np.zeros((2, 2)), h).writeto(tmp_path / f"20240509_{s}_n4c1A.fts")
    fr = _frames(tmp_path, "cor1", "A")
    assert len(fr) == 1 and len(fr[0][1]) == 3 and _frames(tmp_path, "cor1", "B") == []


def test_a_gap_in_the_frames_is_refused_and_the_nearest_frame_times_are_named(tmp_path):
    _pair(tmp_path)                                         # frames at 14:39:51 and 15:24:51
    wide = ["2013-10-25T14:35", "2013-10-25T18:00"]
    with pytest.raises(ValueError, match=r"within 30 min.*2013-10-25T15:24:51 \(before\)"):
        load_plane("cor2", tmp_path, "2013-10-25T17:55", OBS, wide)
    with pytest.raises(FileNotFoundError, match=r"nearest available frames: 2013-10-25T15:24:51 \(before\)"):
        load_plane("cor2", tmp_path, "2013-10-25T17:55", OBS, ["2013-10-25T17:00", "2013-10-25T18:00"])
    _write(tmp_path, "20131025_175900_d4c2B.fts", "2013-10-25T17:59:51", roll=0.0)
    with pytest.raises(ValueError, match=r"15:24:51 \(before\) and 2013-10-25T17:59:51 \(after\)"):
        load_plane("cor2", tmp_path, "2013-10-25T16:40", OBS, wide)
