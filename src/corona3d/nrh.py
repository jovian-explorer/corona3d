"""Images of the Nancay Radioheliograph (NRH2 2-D files, nrh2_FFFF_..._i.fts).

Each file holds the Stokes I brightness temperature [K] of one frequency as a binary table of
64 x 64 images, one row per integration (about 10 s). The images are rotated to solar north up
(ROT_P = 1), with x towards the west and y towards the north; Sun centre lies at pixel CRPIX
(FITS convention, counted from 1) and one solar radius spans SOLAR_R pixels.
"""
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
from astropy.io import fits


def read(path):
    """Times (datetime, UTC), frequency [MHz], images (n, 64, 64) [K], Sun centre (x, y) as
    0-based pixels, and pixels per solar radius."""
    with fits.open(path) as h:
        hd, d = h[1].header, h[1].data
        day = datetime.fromisoformat(h[0].header["DATE-OBS"][:10])
        times = [day + timedelta(milliseconds=int(ms)) for ms in d["TIME"]]
        return dict(times=times, freq=float(hd["FREQ"]), images=np.asarray(d["STOKESI"], float),
                    centre=(hd["CRPIX1"] - 1.0, hd["CRPIX2"] - 1.0), pix_rsun=float(hd["SOLAR_R"]),
                    exptime=float(hd["EXPTIME"]))


def files(folder):
    """{frequency [MHz]: [paths]} of the NRH2 image files in a folder."""
    out = {}
    for f in sorted(Path(folder).glob("nrh2_*_i.f*ts")):
        out.setdefault(int(f.name.split("_")[1]) / 10.0, []).append(f)
    return out


def coverage(folder, freq):
    """[(first, last)] frame times (datetime, UTC) of each NRH2 file of `freq` [MHz] in a folder."""
    out = []
    for p in files(folder).get(round(freq, 1), []):
        with fits.open(p) as h:
            day = datetime.fromisoformat(h[0].header["DATE-OBS"][:10])
            ms = np.asarray(h[1].data["TIME"], float)
        out.append((day + timedelta(milliseconds=float(ms.min())), day + timedelta(milliseconds=float(ms.max()))))
    return out


def image_at(folder, freq, time):
    """The image of `freq` [MHz] nearest to `time` (datetime), from the files that cover it;
    returns (image [K], frame time, x and y of every pixel in R_sun, x west and y north)."""
    paths = files(folder).get(round(freq, 1))
    if not paths:
        raise FileNotFoundError(f"no NRH file at {freq} MHz in {folder}")
    for p in paths:
        r = read(p)
        if r["times"][0] - timedelta(seconds=r["exptime"]) <= time <= r["times"][-1] + timedelta(seconds=r["exptime"]):
            i = int(np.argmin([abs((t - time).total_seconds()) for t in r["times"]]))
            ny, nx = r["images"].shape[1:]
            x = (np.arange(nx) - r["centre"][0]) / r["pix_rsun"]
            y = (np.arange(ny) - r["centre"][1]) / r["pix_rsun"]
            return r["images"][i], r["times"][i], x, y
    spans = ", ".join(f"{read(p)['times'][0]:%H:%M:%S}-{read(p)['times'][-1]:%H:%M:%S}" for p in paths)
    raise ValueError(f"no NRH image at {freq} MHz near {time:%H:%M:%S} UT; the files cover {spans} UT")


def source(img, x, y, level=0.5):
    """Contours at `level` times the peak brightness temperature, as lists of (x, y) in R_sun,
    and the brightness-weighted centroid of the pixels above that level."""
    import contourpy
    lev = level * np.nanmax(img)
    gen = contourpy.contour_generator(x, y, img)
    lines = [np.asarray(c) for c in gen.lines(lev) if len(c) > 2]
    w = np.where(img >= lev, img, 0.0)
    X, Y = np.meshgrid(x, y)
    centroid = (float((w * X).sum() / w.sum()), float((w * Y).sum() / w.sum()))
    return lines, centroid, float(np.nanmax(img))
