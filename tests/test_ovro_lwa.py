"""Synthetic OVRO-LWA spectral-image files: the reader must restore the band images to the full
grid and place a source at its sky position, choose the nearest band and file, refuse times and
frequencies that the files do not cover, and the download must select files by time."""
import io
from datetime import datetime

import numpy as np
import pytest

from corona3d import nrh, ovro_lwa

FREQS = np.array([34.1e6, 52.4e6, 84.6e6])


def write_lwa(path, src=(1.1, 0.05), n=64, small=32, tb=2e7):
    """One file: n x n grid at 90 arcsec per pixel, Sun centre at the reference pixel, each band
    stored at small x small; a Gaussian source at src [R_sun] in every band."""
    import h5py
    rsun = 960.0
    yy, xx = np.mgrid[:small, :small]
    k = n / small                                       # stored pixel = k grid pixels
    cx, cy = (n / 2 + src[0] * rsun / 90.0) / k, (n / 2 + src[1] * rsun / 90.0) / k
    img = 1e5 + tb * np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * 1.0**2))
    with h5py.File(path, "w") as f:
        for i in range(len(FREQS)):
            f[f"FITS_pol0ch{i:04d}"] = img
        cv = f.create_dataset("ch_vals", data=np.vstack([FREQS, np.full(3, 4.6e6), np.ones((3, 3))]).astype("f4"))
        cv.attrs.update(NAXIS1=n, NAXIS2=n, CRPIX1=n / 2 + 1, CRPIX2=n / 2 + 1, CDELT1=90.0, CDELT2=90.0,
                        CRVAL1=0.0, CRVAL2=0.0, RSUN_OBS=rsun, CTYPE1="HPLN-TAN", CTYPE2="HPLT-TAN",
                        original_shape=np.array([1, 3, n, n]))
        cv.attrs["arr_name"] = np.array(["cfreqs", "cdelts", "bmaj", "bmin", "bpa"], dtype=object)


def name(t):
    return f"ovro-lwa-352.lev1_mfs_10s.{t:%Y-%m-%dT%H%M%S}Z.image_I.hdf"


def test_reader_restores_the_grid_and_places_the_source(tmp_path):
    for t in (datetime(2024, 5, 14, 14, 10, 2), datetime(2024, 5, 14, 14, 11, 2)):
        write_lwa(tmp_path / name(t))
    img, tf, x, y = ovro_lwa.image_at(tmp_path, 52.0, datetime(2024, 5, 14, 14, 10, 50))
    assert tf == datetime(2024, 5, 14, 14, 11, 2) and img.shape == (64, 64)
    _, (cx, cy), tb = nrh.source(img, x, y, 0.5)
    assert (cx, cy) == pytest.approx((1.1, 0.05), abs=0.1)     # x west, y north, R_sun
    with pytest.raises(ValueError, match="bands are 34.1, 52.4, 84.6"):
        ovro_lwa.image_at(tmp_path, 70.0, datetime(2024, 5, 14, 14, 10, 2))
    with pytest.raises(ValueError, match="within 2 min"):
        ovro_lwa.image_at(tmp_path, 52.4, datetime(2024, 5, 14, 15, 0))


def test_download_selects_files_by_time(tmp_path, monkeypatch):
    import urllib.request

    from corona3d import fetch
    keys = [f"image_hdf/2024/05/14/{name(datetime(2024, 5, 14, 14, m, 2))}" for m in (5, 10, 11, 30)]
    monkeypatch.setattr(ovro_lwa, "list_remote", lambda day, product="mfs": keys)
    got = []
    monkeypatch.setattr(urllib.request, "urlopen", lambda url, timeout: got.append(url) or io.BytesIO(b"x"))
    cfg = {"radio_images": [dict(instrument="ovro_lwa", folder="data/ovro_lwa", start="2024-05-14T14:10:00",
                                 end="2024-05-14T14:12:00")]}
    assert fetch.ovro_lwa(tmp_path, cfg) == 0
    assert sorted(p.name for p in (tmp_path / "data" / "ovro_lwa").iterdir()) == [k.split("/")[-1] for k in keys[1:3]]
    assert all(u.startswith(ovro_lwa.S3) for u in got)
