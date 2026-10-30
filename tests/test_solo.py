"""Solar Orbiter as an observer, offline with synthetic headers: the time of a frame seen from 1 AU, the
Metis field from the header, the position of the spacecraft in config.json, the recognition of EUI/FSI
images, the filter as the observing programme, the Metis plane in three dimensions, and the panel seen
from Solar Orbiter."""
import json
from datetime import datetime, timedelta

import astropy.units as u
import numpy as np
import pytest
import sunpy.map
from astropy.coordinates import SkyCoord
from astropy.io import fits
from sunpy.coordinates import frames, sun

from corona3d import euv, gcs, solo
from corona3d.geometry import observer_basis, unit_vec

DSUN = 0.688 * solo.AU_M
WIN = ["2024-05-09T16:32:00", "2024-05-09T20:02:00"]


def _metis(d, name, t_obs, t_avg, blob=None, n=256, scale=81.1):
    """A Metis tB file seen from Stonyhurst lon 167.3, lat 7.5 at 0.688 AU: zero outside the field (as the
    level-2 files), a uniform corona of 1 MSB inside it, and a bright square at blob (x, y) R_sun."""
    obs = SkyCoord(167.3 * u.deg, 7.5 * u.deg, DSUN * u.m, frame=frames.HeliographicStonyhurst, obstime=t_obs)
    ref = SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame=frames.Helioprojective(observer=obs, obstime=t_obs))
    hdr = sunpy.map.make_fitswcs_header((n, n), ref, scale=[scale, scale] * u.arcsec / u.pix, rotation_angle=-13.8 * u.deg)
    m = sunpy.map.Map(np.zeros((n, n)), hdr)
    rs = m.rsun_obs.to_value(u.arcsec)
    hdr.update({"INSTRUME": "Metis", "OBSRVTRY": "Solar Orbiter", "TELESCOP": "SOLO/Metis/VLD", "DATE-AVG": t_avg,
                "INN_FOV": 1.6, "OUT_FOV": 3.4, "RSUN_ARC": rs, "XPOSURE": 420.0, "BUNIT": "MSB"})
    yy, xx = np.mgrid[:n, :n]
    cx, cy = (float(v) for v in m.wcs.world_to_pixel(ref))
    r = np.hypot(xx - cx, yy - cy) * scale / rs
    data = np.where((r > 1.6 * 3600 / rs) & (r < 3.4 * 3600 / rs), 1.0, 0.0)
    if blob is not None:
        c = SkyCoord(blob[0] * rs * u.arcsec, blob[1] * rs * u.arcsec, frame=m.coordinate_frame)
        i, j = (int(round(float(v))) for v in m.wcs.world_to_pixel(c))
        data[j - 2:j + 2, i - 2:i + 2] += 100.0
    fits.PrimaryHDU(data, header=fits.Header(hdr)).writeto(d / name)
    return d / name


def test_frame_time_is_the_middle_of_the_sequence_seen_from_earth():
    from sunpy.coordinates import get_earth
    h = {"DATE-OBS": "2024-05-09T18:29:55.426", "DATE-AVG": "2024-05-09T18:44:42.349", "DSUN_OBS": 102961478291.0}
    d_earth = float(get_earth(datetime(2024, 5, 9, 18, 44, 42)).radius.to_value("m"))   # 1.0108 AU in May
    assert d_earth == pytest.approx(solo.AU_M, rel=0.02) and d_earth != solo.AU_M
    shift = (d_earth - 102961478291.0) / solo.C_M_S                    # 161 s, not the 155.6 s of 1 AU
    expect = datetime(2024, 5, 9, 18, 44, 42, 349000) + timedelta(seconds=shift) + timedelta(microseconds=500000)
    assert solo.frame_time(h) == expect.replace(microsecond=0)
    near = solo.frame_time({"DATE-OBS": "2024-05-09T17:40:45.211", "DSUN_OBS": solo.AU_M})
    assert 0 < (near - datetime(2024, 5, 9, 17, 40, 45)).total_seconds() < 10     # Earth is 1.0108 AU away


def test_metis_field_scales_with_the_distance_of_the_spacecraft():
    rs = 1393.68855646                                                  # arcsec, 0.688 AU on 2024-05-09
    r_in, r_out, r_seed = solo.metis_fov({"RSUN_ARC": rs, "INN_FOV": 1.6, "OUT_FOV": 3.4})
    assert (r_in, r_out) == pytest.approx((1.6 * 3600 / rs, 3.4 * 3600 / rs), abs=1e-3)
    assert r_in == pytest.approx(4.13, abs=0.01) and r_in < r_seed < r_out
    near = solo.metis_fov({"DSUN_OBS": 0.3 * solo.AU_M})                 # nominal 1.6 to 2.9 deg, closer to the Sun
    assert near[0] == pytest.approx(1.80, abs=0.01) and near[1] == pytest.approx(3.27, abs=0.01)


def test_observer_from_the_header_is_the_carrington_position_at_the_event_time():
    h = {"HGLN_OBS": 167.3365, "HGLT_OBS": 7.527, "DSUN_OBS": DSUN, "DATE-OBS": "2024-05-09T16:29:55"}
    lon, lat, au = solo.observer(h, "2024-05-09T17:32:00")
    earth = sun.L0("2024-05-09T17:32:00").to_value(u.deg)              # Carrington longitude of Earth
    assert (lon - earth - 167.3365 + 180) % 360 - 180 == pytest.approx(0.0, abs=0.05)
    assert lat == pytest.approx(7.527, abs=1e-3) and au == pytest.approx(0.688, abs=1e-4)


def test_fetched_files_set_the_observer_in_config_json(tmp_path):
    (tmp_path / "config.json").write_text(json.dumps({"time": "2024-05-09T17:32:00"}), encoding="utf-8")
    (tmp_path / "data" / "wl").mkdir(parents=True)
    _metis(tmp_path / "data" / "wl", "solo_L2_metis-vl-tb_20240509T162955_V01.fits", "2024-05-09T16:29:55", "2024-05-09T16:44:42")
    lon, lat, au = solo.set_observer(tmp_path)
    cfg = json.loads((tmp_path / "config.json").read_text(encoding="utf-8"))
    assert cfg["solo"] == [lon, lat] and cfg["solo_dsun_au"] == pytest.approx(0.688)
    assert lat == pytest.approx(7.5, abs=1e-3)
    assert solo.set_observer(tmp_path / "data") is None                  # a folder without Solar Orbiter files


def test_eui_fsi_headers_are_recognised_and_timed(tmp_path):
    assert euv.key_of({"TELESCOP": "SOLO/EUI/FSI", "DETECTOR": "FSI", "WAVELNTH": 174}) == "fsi174"
    assert euv.key_of({"TELESCOP": "SOLO/EUI/FSI", "DETECTOR": "FSI", "WAVELNTH": 304}) is None
    hdr = fits.Header({"TELESCOP": "SOLO/EUI/FSI", "DETECTOR": "FSI", "INSTRUME": "EUI", "WAVELNTH": 174, "FILTER": "Zirconium_174_n25",
                       "DATE-OBS": "2024-05-09T17:40:45.211", "DATE-AVG": "2024-05-09T17:40:50.211", "DSUN_OBS": DSUN})
    f = tmp_path / "solo_L2_eui-fsi174-image_20240509T174045211_V02.fits"
    fits.HDUList([fits.PrimaryHDU(), fits.CompImageHDU(np.ones((8, 8), np.float32), header=hdr)]).writeto(f)
    fr = euv.frames(tmp_path)
    assert list(fr) == ["fsi174"] and fr["fsi174"][0][0] == datetime(2024, 5, 9, 17, 43, 31)
    assert gcs.seq(f) == "Zirconium_174_n25"                           # the two FSI filters are separate programmes


def test_metis_frames_set_the_field_and_are_not_divided_by_the_exposure(tmp_path):
    f = _metis(tmp_path, "solo_L2_metis-vl-tb_20240509T182955_V01.fits", "2024-05-09T18:29:55", "2024-05-09T18:44:42")
    fr = gcs.frames_of(tmp_path, "metis")
    assert [t for t, _ in fr] == [datetime(2024, 5, 9, 18, 47, 23)]
    assert gcs.FOV["metis"][0] == pytest.approx(4.13, abs=0.01) and gcs.height_limit(["metis"]) == pytest.approx(1.5 * gcs.FOV["metis"][1])
    m = solo.load_map(f)
    assert np.nanmax(m.data) == pytest.approx(1.0) and np.isnan(m.data).any()    # MSB kept (XPOSURE 420 s), zero pixels outside the field
    assert gcs._t(f) == fr[0][0]


@pytest.mark.parametrize("mode", ["ratio", "difference"])
def test_metis_plane_places_a_feature_in_3d(tmp_path, mode):
    pytest.importorskip("pyvista")
    from corona3d.whitelight import load_plane
    _metis(tmp_path, "solo_L2_metis-vl-tb_20240509T162955_V01.fits", "2024-05-09T16:29:55", "2024-05-09T16:44:42")
    f = _metis(tmp_path, "solo_L2_metis-vl-tb_20240509T182955_V01.fits", "2024-05-09T18:29:55", "2024-05-09T18:44:42", blob=(-5.0, 2.0))
    obs = solo.observer(solo.header_of(f), "2024-05-09T18:44:42")[:2]
    g, tf, tb = load_plane("metis", tmp_path, "2024-05-09T18:47", obs, WIN, mode)
    v = np.asarray(g["wl"])
    p = g.points[np.nanargmax(v)]
    ex, ey, ez = observer_basis(*obs)
    assert p @ ex == pytest.approx(-5.0, abs=0.2)                      # east of Sun centre seen from Solar Orbiter
    assert p @ ey == pytest.approx(2.0, abs=0.2)
    assert abs(p @ ez) < 1e-9
    r = np.hypot(g.points @ ex, g.points @ ey)
    assert np.isnan(v[r < 4.0]).all()                                  # Metis occulter left transparent
    assert (tf.hour, tf.minute, tb.hour, tb.minute) == (18, 47, 16, 47)


def test_panel_seen_from_solar_orbiter():
    pytest.importorskip("pyvista")
    from corona3d.reconstruction import parse_view, view_geometry
    assert parse_view("c=solo") == (2, {"from": "solo"})
    p = {"earth": [317.145, -3.229], "solo": [124.4, 7.53]}
    d, up, owner = view_geometry({"from": "solo"}, p)
    assert owner == "solo" and d == pytest.approx(unit_vec(124.4, 7.53))
    assert up == pytest.approx(observer_basis(124.4, 7.53)[1])          # solar north up
    with pytest.raises(ValueError, match="fetch --solo"):
        view_geometry({"from": "solo"}, {"earth": [317.145, -3.229]})


def test_cli_offers_the_solar_orbiter_download_and_layer(capsys):
    from corona3d.cli import main
    from corona3d.layers import LAYERS, ON
    with pytest.raises(SystemExit):
        main(["fetch", "--help"])
    out = capsys.readouterr().out
    assert "--solo" in out and "--solohi" in out
    assert "metis" in LAYERS and ON["metis"] is False                   # off by default: default figures unchanged
