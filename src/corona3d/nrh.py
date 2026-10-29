"""Images of the Nancay Radioheliograph (NRH2 2-D files, nrh2_FFFF_..._i.fts).

Each file holds the Stokes I brightness temperature [K] of one frequency as a binary table of
64 x 64 images, one row per integration (about 10 s). The images are rotated to solar north up
(ROT_P = 1), with x towards the west and y towards the north; Sun centre lies at pixel CRPIX
(FITS convention, counted from 1) and one solar radius spans SOLAR_R pixels. The files of the Radio
Solar Database also hold the Stokes V images in the column STOKESV, from which the degree of circular
polarisation of a source is measured (stokes_v_at, polarisation).
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


# ------------------------------------------------------------------ circular polarisation (Stokes V)
def v_file(path):
    """The Stokes V file that belongs to the Stokes I file `path` (nrh2_..._v.fts beside nrh2_..._i.fts),
    or None."""
    p = Path(path)
    for ext in (".fts", ".fits"):
        q = p.with_name(p.name.rsplit("_i.", 1)[0] + "_v" + ext)
        if q.exists():
            return q
    return None


def stokes_v_at(folder, freq, time):
    """(Stokes V image [K], its source) of the frame of `freq` [MHz] nearest to `time`, from the column
    STOKESV of the Stokes I file (the NRH2 files of the Radio Solar Database hold both) or from a separate
    _v file; (None, reason) when no Stokes V is available."""
    for p in files(folder).get(round(freq, 1), []):
        with fits.open(p) as h:
            d, day = h[1].data, datetime.fromisoformat(h[0].header["DATE-OBS"][:10])
            times = [day + timedelta(milliseconds=int(ms)) for ms in d["TIME"]]
            exp = float(h[1].header["EXPTIME"])
            if not times[0] - timedelta(seconds=exp) <= time <= times[-1] + timedelta(seconds=exp):
                continue
            i = int(np.argmin([abs((t - time).total_seconds()) for t in times]))
            if "STOKESV" in d.columns.names:
                return np.asarray(d["STOKESV"][i], float), f"{p.name} (column STOKESV)"
        q = v_file(p)
        if q is None:
            return None, f"{p.name} has no STOKESV column and no _v file"
        with fits.open(q) as h:
            d = h[1].data
            col = "STOKESV" if "STOKESV" in d.columns.names else d.columns.names[-1]
            return np.asarray(d[col][i], float), f"{q.name} (column {col})"
    return None, f"no NRH file at {freq} MHz covering {time:%H:%M:%S} UT"


def polarisation(img_i, img_v, level=0.5):
    """Degree of circular polarisation V/I of a source: (sum V / sum I over the pixels with I at or above
    `level` times the peak, V / I at the peak pixel). The sign is that of the NRH files."""
    img_i, img_v = np.asarray(img_i, float), np.asarray(img_v, float)
    w = img_i >= level * np.nanmax(img_i)
    k = np.unravel_index(np.nanargmax(img_i), img_i.shape)
    return float(img_v[w].sum() / img_i[w].sum()), float(img_v[k] / img_i[k])


def polarisation_rows(home, p, rt):
    """Columns dcp_source, dcp_peak, stokes_v (the file it came from), and dcp_valid for the rows of a
    radio table; NaN for rows that are not NRH images or have no Stokes V. A degree of polarisation can
    not exceed 1 in magnitude, so that a value |V/I| > 1 (a Stokes V map whose artefacts exceed the
    source, as at 432 MHz on 2013-10-25) is kept in the table with dcp_valid false."""
    from .whitelight import utc
    n = len(rt["x_rsun"])
    src, pk, how = np.full(n, np.nan), np.full(n, np.nan), [""] * n
    for j, (inst, f, t) in enumerate(zip(rt["instrument"], rt["freq_mhz"], rt["time_utc"])):
        if inst != "nrh" or not t:
            continue
        e = next((e for e in p["radio_images"] if e["instrument"] == "nrh" and any(abs(f - x) < 0.05 for x in e["freqs"])), None)
        if e is None:
            continue
        folder = Path(home) / e["folder"]
        tt = utc(t)
        img, _, _, _ = image_at(folder, f, tt)
        v, how[j] = stokes_v_at(folder, f, tt)
        if v is not None:
            src[j], pk[j] = polarisation(img, v, e["level"])
    valid = np.isfinite(src) & (np.abs(src) <= 1) & (np.abs(pk) <= 1)
    return dict(dcp_source=src, dcp_peak=pk, dcp_valid=valid, stokes_v=np.array(how, dtype=str))
