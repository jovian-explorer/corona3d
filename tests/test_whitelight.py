"""A synthetic coronagraph pair: the base-difference plane must put a feature at its true place
in three dimensions, whatever the roll of the image, and the COR1 triplets must be grouped."""
import astropy.units as u
import numpy as np
import pytest
import sunpy.map
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames

from corona3d.geometry import observer_basis
from corona3d.whitelight import _frames, load_plane

OBS = (215.329, -0.984)                  # STEREO-B, Carrington lon/lat on 2013-10-25


def _write(d, name, time, blob=None, roll=0.0, n=512, scale=30.0):
    obs = SkyCoord(-142.0 * u.deg, OBS[1] * u.deg, 1 * u.AU, frame=frames.HeliographicStonyhurst,
                   obstime=time)
    ref = SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame=frames.Helioprojective(observer=obs, obstime=time))
    hdr = sunpy.map.make_fitswcs_header((n, n), ref, scale=[scale, scale] * u.arcsec / u.pix,
                                        rotation_angle=roll * u.deg)
    hdr["exptime"] = 6.0
    m = sunpy.map.Map(np.zeros((n, n)), hdr)
    data = np.ones((n, n))
    if blob is not None:                 # a 2-pixel square at (x, y) R_sun in solar north-up axes
        rs = m.rsun_obs.to_value(u.arcsec)
        c = SkyCoord(blob[0] * rs * u.arcsec, blob[1] * rs * u.arcsec, frame=m.coordinate_frame)
        i, j = (int(round(float(v))) for v in m.wcs.world_to_pixel(c))
        data[j - 2:j + 2, i - 2:i + 2] += 100.0
    sunpy.map.Map(data, m.meta).save(str(d / name), overwrite=True)


@pytest.mark.parametrize("roll", [0.0, 6.6, 175.3])
def test_plane_places_a_feature_in_3d(tmp_path, roll):
    _write(tmp_path, "20131025_143900_d4c2B.fts", "2013-10-25T14:39:51", roll=roll)
    _write(tmp_path, "20131025_152400_d4c2B.fts", "2013-10-25T15:24:51", blob=(5.0, 2.0), roll=roll)
    g = load_plane("cor2", tmp_path, "2013-10-25T15:24", OBS)
    v = np.asarray(g["wl"])
    p = g.points[np.nanargmax(v)]
    ex, ey, ez = observer_basis(*OBS)
    assert p @ ex == pytest.approx(5.0, abs=0.15)          # west
    assert p @ ey == pytest.approx(2.0, abs=0.15)          # north
    assert abs(p @ ez) < 1e-9                              # in the plane of the sky
    r = np.hypot(g.points @ ex, g.points @ ey)
    assert np.isnan(v[r < 3.9]).all()                      # COR2 occulter left transparent


def test_cor1_triplets_are_grouped(tmp_path):
    from astropy.io import fits
    for s in ("144500", "144512", "144524", "152500", "152512"):   # second triplet incomplete
        h = fits.Header({"DATE-OBS": f"2013-10-25T{s[:2]}:{s[2:4]}:{s[4:]}.000"})
        fits.PrimaryHDU(np.zeros((2, 2)), h).writeto(tmp_path / f"20131025_{s}_s4c1B.fts")
    fr = _frames(tmp_path, "cor1")
    assert len(fr) == 1 and len(fr[0][1]) == 3
    assert fr[0][0].strftime("%H:%M:%S") == "14:45:00"
