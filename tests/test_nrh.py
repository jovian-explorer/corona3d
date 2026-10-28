"""Synthetic Nancay Radioheliograph files: the reader must place a source at its true sky
position, the radio layer must turn a burst into a chain of blobs at the Newkirk height of each
frequency, and times outside the files must be refused."""
from datetime import datetime

import numpy as np
import pytest
from astropy.io import fits

from corona3d import nrh
from corona3d.geometry import observer_basis
from corona3d.reconstruction import newkirk_r, reconstruct

EARTH = (357.372, 4.993)


def write_nrh(path, freq, start_ms, n, src=None, tb=1e9):
    """64 x 64 images every 10 s from start_ms (ms after midnight); src(k) -> (x, y) [R_sun] of a
    Gaussian source in image k, or None for a quiet image."""
    yy, xx = np.mgrid[:64, :64]
    imgs = np.full((n, 64, 64), 1e5, np.float32)
    for k in range(n):
        if src is not None and src(k) is not None:
            x, y = src(k)
            imgs[k] += tb * np.exp(-((xx - 31 - 16 * x) ** 2 + (yy - 31 - 16 * y) ** 2) / (2 * 1.5**2))
    cols = [fits.Column("TIME", "1J", array=start_ms + 9900 * np.arange(n)),
            fits.Column("NUMIM", "1J", array=np.arange(n)),
            fits.Column("STOKESI", "4096E", dim="(64,64)", array=imgs),
            fits.Column("STOKESV", "4096E", dim="(64,64)", array=np.zeros_like(imgs))]
    t = fits.BinTableHDU.from_columns(cols)
    t.header.update(FREQ=freq, CRPIX1=32, CRPIX2=32, SOLAR_R=16.0, EXPTIME=9.9)
    fits.HDUList([fits.PrimaryHDU(header=fits.Header({"DATE-OBS": "2013-10-25"})), t]).writeto(path)


T0 = (14 * 3600 + 58 * 60) * 1000                     # 14:58:00 UT


def test_reader_places_a_source_at_its_sky_position(tmp_path):
    write_nrh(tmp_path / "nrh2_1509_h60_20131025_145800c00_i.fts", 150.9, T0, 30, lambda k: (-1.1, -0.2))
    img, tf, x, y = nrh.image_at(tmp_path, 150.9, datetime(2013, 10, 25, 15, 0, 1))
    assert tf == datetime(2013, 10, 25, 14, 59, 58, 800000)  # nearest image: 14:58:00 + 12 x 9.9 s
    lines, (cx, cy), tb = nrh.source(img, x, y, 0.5)
    assert (cx, cy) == pytest.approx((-1.1, -0.2), abs=0.02)  # x west, y north, R_sun
    assert 0.9e9 < tb <= 1e9 + 1e5                            # peak pixel, source between pixels
    assert lines and all(np.hypot(c[:, 0] + 1.1, c[:, 1] + 0.2).max() < 0.2 for c in lines)
    with pytest.raises(ValueError, match="files cover 14:58:00"):
        nrh.image_at(tmp_path, 150.9, datetime(2013, 10, 25, 15, 10))
    with pytest.raises(FileNotFoundError):
        nrh.image_at(tmp_path, 228.0, datetime(2013, 10, 25, 15, 0))


def test_radio_layer_draws_a_chain_of_blobs_per_frequency(synthetic_home):
    d = synthetic_home / "data" / "nrh"
    d.mkdir()
    # 150.9 MHz: a source moving south-east over the burst; 432 MHz: quiet until 15:01, then a source
    write_nrh(d / "nrh2_1509_h60_20131025_145800c00_i.fts", 150.9, T0, 40, lambda k: (-1.0 - 0.01 * k, -0.1 - 0.01 * k))
    write_nrh(d / "nrh2_4320_h60_20131025_145800c00_i.fts", 432.0, T0, 40,
              lambda k: (-0.85, -0.2) if k >= 18 else None)
    on = dict(surface=False, sphere=False, open=False, closed=False, streamer=False, gcs=False, radio=True)
    entry = dict(instrument="nrh", folder="data/nrh", freqs=[150.9, 432.0], start="2013-10-25T14:59:00",
                 end="2013-10-25T15:03:00", step_s=60, level=0.5, min_tb=1e8, contours=True)
    sc = reconstruct({"radio_images": [entry]}, on, synthetic_home)
    (f1, xyz1, rgb1), (f2, xyz2, _) = sc["radio"]
    assert (f1, f2) == (150.9, 432.0)
    assert len(xyz1) == 5 and len(xyz2) == 3                 # quiet 432 MHz images left out
    assert np.allclose(np.linalg.norm(xyz1, axis=1), newkirk_r(150.9))
    assert np.allclose(np.linalg.norm(xyz2, axis=1), newkirk_r(432.0))
    ex, _, ez = observer_basis(*EARTH)
    assert xyz2[0] @ ex == pytest.approx(-0.85, abs=0.02) and (xyz2 @ ez > 0).all()
    assert not np.allclose(rgb1[0], rgb1[-1])                # coloured from early to late
    assert len(sc["radio_contours"]) == 2


def test_contours_are_drawn_at_every_level_and_nest_inwards(synthetic_home):
    d = synthetic_home / "data" / "nrh"
    d.mkdir()
    write_nrh(d / "nrh2_1509_h60_20131025_145800c00_i.fts", 150.9, T0, 40, lambda k: (-1.0, -0.1))
    on = dict(surface=False, sphere=False, open=False, closed=False, streamer=False, gcs=False, radio=True)
    entry = dict(instrument="nrh", folder="data/nrh", freqs=[150.9], start="2013-10-25T14:59:00",
                 end="2013-10-25T15:01:00", step_s=60, level=0.5, min_tb=1e8, contours=True,
                 contour_levels=[0.9, 0.3, 0.6])
    (f, levels), = reconstruct({"radio_images": [entry]}, on, synthetic_home)["radio_contours"]
    assert [lev for lev, _ in levels] == [0.3, 0.6, 0.9]          # from the lowest to the highest
    ex, ey, _ = observer_basis(*EARTH)
    centre = np.array([-1.0, -0.1])
    size = [max(np.hypot(c @ ex - centre[0], c @ ey - centre[1]).max() for c in lines) for _, lines in levels]
    assert size[0] > size[1] > size[2]                            # higher levels lie inside lower ones
    entry.pop("contour_levels")
    (_, levels), = reconstruct({"radio_images": [entry]}, on, synthetic_home)["radio_contours"]
    assert [lev for lev, _ in levels] == [0.5]                    # without levels: the level of the centroid
