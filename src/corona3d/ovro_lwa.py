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


# ------------------------------------------------------------------ position check on the quiet Sun
def disc(img, x, y, frac=0.5, r_search=3.0, iterations=5):
    """Centroid (x, y) [R_sun] and equivalent radius [R_sun] of the radio disc of one image: the pixels within
    r_search of the current centre brighter than frac times the median brightness inside 1 R_sun of it, without
    weighting, so that a compact bright source moves the centroid only by its area. Starts at Sun centre and
    iterates; returns (cx, cy, radius, disc brightness [K], peak / disc brightness)."""
    X, Y = np.meshgrid(x, y)
    cx = cy = 0.0
    tb = np.nan
    m = np.zeros(img.shape, bool)
    for _ in range(iterations):
        R = np.hypot(X - cx, Y - cy)
        tb = float(np.nanmedian(img[R < 1.0]))
        m = (img > frac * tb) & (R < r_search)
        if tb <= 0 or not m.any():
            return np.nan, np.nan, np.nan, tb, np.nan
        cx, cy = float(X[m].mean()), float(Y[m].mean())
    area = m.sum() * abs((x[1] - x[0]) * (y[1] - y[0]))
    return cx, cy, float(np.sqrt(area / np.pi)), tb, float(np.nanmax(img) / tb)


QUIET = dict(radius=(0.9, 2.5), contrast=5.0, max_tb=5e6)
# a quiet disc: equivalent radius [R_sun], largest peak / disc brightness, and largest disc brightness [K] (the quiet
# Sun at metre wavelengths is near 1e6 K, so that a brighter "disc" is a noise storm or a burst)


def quiet_offsets(folder, freqs, t0=None, t1=None, tol_mhz=3.0, quiet=QUIET):
    """Offset of the radio disc from Sun centre per band, from the images of a folder between t0 and t1 whose disc is
    quiet (equivalent radius within quiet['radius'], peak below quiet['contrast'] times the disc brightness, and disc
    brightness below quiet['max_tb']).
    Returns {freq: dict(dx, dy, err_x, err_y, n, n_frames, frames)} with dx, dy [R_sun] the median offset and err the
    standard error of the median (1.2533 sigma / sqrt n); n is 0 when no image is quiet."""
    out = {f: dict(rows=[]) for f in freqs}
    fs = [(t, p) for t, p in files(folder) if (t0 is None or t >= t0) and (t1 is None or t <= t1)]
    for t, p in fs:
        hdr, bands, imgs = read(p)
        x, y = axes(hdr)
        for f in freqs:
            k = int(np.argmin(abs(bands - f)))
            if abs(bands[k] - f) > tol_mhz:
                continue
            cx, cy, rad, tb, con = disc(imgs[k], x, y)
            ok = np.isfinite(cx) and quiet["radius"][0] <= rad <= quiet["radius"][1] and con <= quiet["contrast"]                 and tb <= quiet["max_tb"]
            out[f]["rows"].append(dict(time=t.isoformat(), dx=cx, dy=cy, radius=rad, tb_disc=tb, contrast=con, quiet=bool(ok)))
    res = {}
    for f, d in out.items():
        q = [r for r in d["rows"] if r["quiet"]]
        if q:
            dx, dy = np.array([r["dx"] for r in q]), np.array([r["dy"] for r in q])
            n = len(q)
            res[f] = dict(dx=float(np.median(dx)), dy=float(np.median(dy)),
                          err_x=float(1.2533 * np.std(dx) / np.sqrt(n)) if n > 1 else np.nan,
                          err_y=float(1.2533 * np.std(dy) / np.sqrt(n)) if n > 1 else np.nan, n=n,
                          n_frames=len(d["rows"]), frames=d["rows"])
        else:
            res[f] = dict(dx=np.nan, dy=np.nan, err_x=np.nan, err_y=np.nan, n=0, n_frames=len(d["rows"]), frames=d["rows"])
    return res


def refraction_fit(res):
    """Least-squares fit of the offsets with a + b / f^2 (ionospheric refraction scales as f^-2), per axis;
    returns {axis: (a, b)} [R_sun, R_sun MHz^2] from the bands with quiet images, or {} for fewer than three bands."""
    fs = [f for f, r in res.items() if r["n"]]
    if len(fs) < 3:
        return {}
    A = np.c_[np.ones(len(fs)), 1.0 / np.array(fs) ** 2]
    out = {}
    for ax in ("dx", "dy"):
        coef, *_ = np.linalg.lstsq(A, np.array([res[f][ax] for f in fs]), rcond=None)
        out[ax] = (float(coef[0]), float(coef[1]))
    return out


OFFSET_FILE = "analysis/ovro_lwa_offset.json"


def offset_for(home, entry, freq, tol_mhz=3.0):
    """(dx, dy) [R_sun] to subtract from the OVRO-LWA positions of `freq` for a radio_images entry: its
    `position_offset` as [dx, dy] given by hand, or "quiet_sun" for the offsets measured with `corona3d sources
    --ovro-offset` (analysis/ovro_lwa_offset.json); (0, 0) when the key is absent or the band has no quiet image."""
    import json
    po = entry.get("position_offset")
    if not po:
        return 0.0, 0.0
    if isinstance(po, (list, tuple)):
        return float(po[0]), float(po[1])
    f = Path(home) / OFFSET_FILE
    if po != "quiet_sun" or not f.exists():
        raise ValueError(f"position_offset {po!r}: give [dx, dy] in R_sun, or \"quiet_sun\" after `corona3d sources "
                         "--ovro-offset` has written " + OFFSET_FILE)
    rec = json.loads(f.read_text(encoding="utf-8"))["bands"]
    best = min(rec, key=lambda k: abs(float(k) - freq))
    r = rec[best]
    if abs(float(best) - freq) > tol_mhz or not r["n"]:
        return 0.0, 0.0
    return float(r["dx"]), float(r["dy"])


def quiet_dir(h, entry):
    return Path(h) / (entry["folder"].rstrip("/") + "_quiet")


def fetch_quiet(h, entry, start, end, step_s=600):
    """Download one image file every step_s seconds between start and end (datetimes) from the public bucket
    into the folder of the entry with the suffix _quiet (kept apart from the burst images), for the position
    check; returns the number of files that failed."""
    import urllib.request
    d = quiet_dir(h, entry)
    d.mkdir(parents=True, exist_ok=True)
    keys = [k for day in sorted({start.date(), end.date()}) for k in list_remote(day)]
    keys = sorted((file_time(k), k) for k in keys if start <= file_time(k) <= end)
    chosen, last = [], None
    for t, k in keys:
        if last is None or (t - last).total_seconds() >= step_s:
            chosen.append(k)
            last = t
    failed = 0
    for k in chosen:
        f = d / k.split("/")[-1]
        if f.exists():
            continue
        tmp = f.with_suffix(".part")
        try:
            with urllib.request.urlopen(f"{S3}/{k}", timeout=120) as r, open(tmp, "wb") as o:
                while chunk := r.read(1 << 20):
                    o.write(chunk)
            tmp.replace(f)
        except OSError as e:
            tmp.unlink(missing_ok=True)
            print(f"{f.name}: failed ({e})")
            failed += 1
    return failed
