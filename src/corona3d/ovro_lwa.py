"""Spectral images of the Owens Valley Radio Observatory Long Wavelength Array (OVRO-LWA).

Each file `ovro-lwa-352.lev1[.5]_{mfs,fch}_10s.YYYY-MM-DDTHHMMSSZ.image_I.hdf` holds one 10 s
integration: the FITS header as attributes of `ch_vals`, the centre frequency of every band in
`ch_vals[0]`, and one Stokes I brightness-temperature image [K] per band (`FITS_pol0chNNNN`),
stored at reduced size and restored to the full grid (NAXIS1 x NAXIS2) by interpolation, as in
lwa-solar-util (recover_fits_from_h5). The images are helioprojective (HPLN-TAN, HPLT-TAN), with
solar north up, x towards the west, and y towards the north.
"""
import re
from datetime import datetime
from pathlib import Path

import numpy as np

S3 = "https://ovro-lwa-solar.s3.us-west-2.amazonaws.com"
NAME = re.compile(r"ovro-lwa-352\.lev1(?:\.5)?_(mfs|fch)_10s\.(\d{4}-\d{2}-\d{2}T\d{6})Z\.image_I\.hdf$")


def file_time(path):
    return datetime.strptime(NAME.search(Path(path).name).group(2), "%Y-%m-%dT%H%M%S")


def files(folder):
    """[(time, path)] of the OVRO-LWA image files in a folder, in time order."""
    return sorted((file_time(f), f) for f in Path(folder).glob("ovro-lwa-352.*image_I.hdf") if NAME.search(f.name))


def read(path):
    """Header (dict), band centre frequencies [MHz], and the full-grid images (n_band, ny, nx) [K]."""
    from ._extras import need
    h5py = need("radio", "h5py", "reading OVRO-LWA images")
    from scipy.ndimage import zoom
    with h5py.File(path, "r") as f:
        hdr = dict(f["ch_vals"].attrs)
        names = list(hdr["arr_name"])
        freqs = np.asarray(f["ch_vals"][names.index("cfreqs")], float) / 1e6
        ny, nx = int(hdr["NAXIS2"]), int(hdr["NAXIS1"])
        imgs = np.empty((len(freqs), ny, nx))
        for k in range(len(freqs)):
            small = f[f"FITS_pol0ch{k:04d}"][:]
            imgs[k] = small[0, 0] if small.shape[0] == 1 else zoom(small, nx / small.shape[-1], order=3,
                                                                     prefilter=False)[:ny, :nx]
    return hdr, freqs, imgs


def axes(hdr):
    """x and y of the pixel centres [R_sun], from the linear part of the TAN projection (the field
    of a few degrees makes the projection correction negligible)."""
    k = 1.0 / float(hdr["RSUN_OBS"])
    x = (np.arange(int(hdr["NAXIS1"])) + 1 - hdr["CRPIX1"]) * hdr["CDELT1"] + hdr["CRVAL1"]
    y = (np.arange(int(hdr["NAXIS2"])) + 1 - hdr["CRPIX2"]) * hdr["CDELT2"] + hdr["CRVAL2"]
    return x * k, y * k


def image_at(folder, freq, time, tol_mhz=3.0):
    """The image of the band nearest to `freq` [MHz] in the file nearest to `time`; returns
    (image [K], file time, x and y in R_sun) like nrh.image_at."""
    fs = files(folder)
    if not fs:
        raise FileNotFoundError(f"no OVRO-LWA image file in {folder}")
    tf, path = min(fs, key=lambda t: abs((t[0] - time).total_seconds()))
    if abs((tf - time).total_seconds()) > 120:
        raise ValueError(f"no OVRO-LWA image within 2 min of {time:%Y-%m-%d %H:%M:%S} UT; the files cover "
                         f"{fs[0][0]:%H:%M:%S}-{fs[-1][0]:%H:%M:%S} UT")
    hdr, freqs, imgs = read(path)
    k = int(np.argmin(abs(freqs - freq)))
    if abs(freqs[k] - freq) > tol_mhz:
        raise ValueError(f"no OVRO-LWA band within {tol_mhz} MHz of {freq} MHz; the bands are "
                         + ", ".join(f"{f:.1f}" for f in freqs) + " MHz")
    x, y = axes(hdr)
    return imgs[k], tf, x, y


def list_remote(day, product="mfs"):
    """Keys of the image files of one day (datetime.date) in the public OVRO-LWA bucket."""
    import urllib.request
    keys, token = [], None
    prefix = f"image_hdf/{day:%Y/%m/%d}/ovro-lwa-352.lev1_{product}_10s"
    while True:
        url = f"{S3}/?list-type=2&prefix={prefix}" + (f"&continuation-token={urllib.request.quote(token)}" if token else "")
        with urllib.request.urlopen(url, timeout=60) as r:
            text = r.read().decode()
        keys += re.findall(r"<Key>([^<]+)</Key>", text)
        m = re.search(r"<NextContinuationToken>([^<]+)</NextContinuationToken>", text)
        if not m:
            return keys
        token = m.group(1)
