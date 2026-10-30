"""Calibration of the model density with the polarised brightness of the corona (`corona3d pb`).

The electrons of the corona scatter the photospheric light by Thomson scattering, and the scattered
light is polarised tangentially to the limb. The polarised brightness pB of a line of sight is an integral
of the electron density along it, weighted by the scattering geometry of van de Hulst (1950) with the limb
darkening of the photosphere, in the form of Billings (1966) and Hayes et al. (2001):

    pB(rho) = K / (1 - u/3) * Int n_e(r) [(1 - u) A(r) + u B(r)] (rho / r)^2 dz
    tB(rho) = K / (1 - u/3) * Int n_e(r) {2 [(1 - u) C(r) + u D(r)] - [(1 - u) A(r) + u B(r)] (rho / r)^2} dz

with K = pi r_e^2 R_sun / 2 (r_e the classical electron radius, dz in R_sun), pB and tB in units of the mean
brightness of the solar disc, rho the distance of the line of sight from Sun centre, r the distance of the
scattering point, u the limb-darkening coefficient, and the coefficients of sin(Omega) = 1/r:

    A = cos(Omega) sin^2(Omega)
    B = -1/8 [1 - 3 sin^2 - cos^2 / sin (1 + 3 sin^2) ln((1 + sin) / cos)]
    C = 4/3 - cos - cos^3 / 3
    D = 1/8 [5 + sin^2 - cos^2 / sin (5 - sin^2) ln((1 + sin) / cos)]

Far from the Sun every coefficient tends to sin^2(Omega), and the integrands to the Thomson scattering of a
point source. pB is linear in n_e, so that the ratio of the observed pB to the pB of the model density
(density_scale 1) is the factor by which the model density must be multiplied.

The observed pB is formed from polarisation sequences of three images: the STEREO/SECCHI COR1 triplets at
0, 120 and 240 deg (level 0.5, calibrated as in secchi_prep: SEB onboard-processing factors, bias, exposure,
the background of the instrumental stray light of the SECCHI archive, the flat field, and the calibration
factor), and the LASCO C2 sequences at +60, 0 and -60 deg (level 0.5, calibrated as in c2_calibrate: bias,
exposure, calibration factor of the polariser, vignetting). With the Mueller matrix of an ideal polariser,
the images I_k at the polariser angles theta_k give I, Q and U, and pB = sqrt(Q^2 + U^2).
"""
import json
import urllib.request
from datetime import datetime
from pathlib import Path

import numpy as np

R_E = 2.8179403262e-13            # cm, classical electron radius
R_SUN_CM = 6.957e10
K_THOMSON = np.pi * R_E**2 / 2 * R_SUN_CM        # cm^3 per R_sun of path: pB = K * Int n_e (...) dz

SSW = "https://sohoftp.nascom.nasa.gov/solarsoft"
SSC_BKG = "https://stereo-ssc.nascom.nasa.gov/data/ins_data/secchi_backgrounds"
DEFAULTS = dict(
    instruments=["cor1", "c2"],
    u=0.63,                       # limb-darkening coefficient of the visible continuum
    r_range=dict(cor1=[1.8, 3.6], c2=[2.6, 5.6]),   # COR1 below 1.8 R_sun: the vignetting edge of the occulter
    r_step=0.1,                   # R_sun, radial bins of the comparison
    pa_step=10.0,                 # deg, position-angle bins
    background="monthly_roll",    # COR1 stray-light background: monthly_roll, monthly_min, or none
    time=None,                    # time of the frames (default: the start of wl_window; of the sequences before the event time, the nearest)
    radial_ranges=[[1.8, 2.4], [2.4, 3.0], [3.0, 4.0], [4.0, 5.6]],
    los_step=0.02, los_r_out=12.0,
    fetch=True,                   # download the C2 polarisation sequences of the day when the folder has none
)


def settings(p):
    s = dict(DEFAULTS, **(p.get("pb") or {}))
    s["r_range"] = dict(DEFAULTS["r_range"], **((p.get("pb") or {}).get("r_range") or {}))
    return s


# ------------------------------------------------------------------ Thomson scattering
def vdh_coefficients(r):
    """van de Hulst coefficients A, B, C, D at the distance r [R_sun] (r > 1)."""
    s = 1.0 / np.clip(np.asarray(r, float), 1.0 + 1e-9, None)
    s2 = s * s
    c = np.sqrt(1.0 - s2)
    with np.errstate(divide="ignore", invalid="ignore"):
        lg = np.where(c > 0, np.log((1.0 + s) / np.where(c > 0, c, 1.0)), 0.0)
        k = np.where(c > 0, c * c / s * lg, 0.0)
    A = c * s2
    B = -(1.0 - 3.0 * s2 - k * (1.0 + 3.0 * s2)) / 8.0
    C = 4.0 / 3.0 - c - c**3 / 3.0
    D = (5.0 + s2 - k * (5.0 - s2)) / 8.0
    return A, B, C, D


def thomson(ne_at, x, y, observer, u=0.63, step=0.02, r_out=12.0, chunk=200):
    """(pB, tB) [mean solar brightness] of the lines of sight through the sky-plane points (x, y) [R_sun]
    of an observer (Carrington longitude, latitude), for the density ne_at(points) [cm^-3]. The line of
    sight runs from -sqrt(r_out^2 - rho^2) to +sqrt(...) through the plane of the sky; points inside the Sun
    and, for rho < 1, behind it are left out."""
    from .geometry import observer_basis
    ex, ey, ez = observer_basis(*observer)
    x, y = np.atleast_1d(np.asarray(x, float)), np.atleast_1d(np.asarray(y, float))
    zmax = np.sqrt(max(r_out**2 - 1.0, 1.0))
    z = np.arange(-zmax, zmax + step / 2, step)
    pb, tb = np.zeros(len(x)), np.zeros(len(x))
    norm = K_THOMSON / (1.0 - u / 3.0) * step
    for i in range(0, len(x), chunk):
        xs, ys = x[i:i + chunk], y[i:i + chunk]
        rho2 = (xs * xs + ys * ys)[:, None]
        r = np.sqrt(rho2 + z[None, :] ** 2)
        ok = (r <= r_out) & (r > 1.0) & ~((rho2 < 1.0) & (z[None, :] < 0))
        pts = xs[:, None, None] * ex + ys[:, None, None] * ey + z[None, :, None] * ez
        ne = np.zeros(r.shape)
        ne[ok] = ne_at(pts[ok])
        A, B, C, D = vdh_coefficients(np.where(ok, r, 2.0))
        ab = (1 - u) * A + u * B
        cd = (1 - u) * C + u * D
        sin2 = rho2 / r**2
        pb[i:i + chunk] = norm * np.sum(np.where(ok, ne * ab * sin2, 0.0), axis=1)
        tb[i:i + chunk] = norm * np.sum(np.where(ok, ne * (2 * cd - ab * sin2), 0.0), axis=1)
    return pb, tb


def stokes_pb(images, angles_deg, k=0.5):
    """(tB, pB) from three images taken through ideal linear polarisers at the given angles [deg]: each image
    is k (I + Q cos 2t + U sin 2t), with k = 1/2 for SECCHI (cor_mueller of secchi_prep) and k = 1 for LASCO,
    whose calibration factor of the polariser already includes its transmission (make_xinv of the LASCO
    library)."""
    a = np.radians(np.asarray(angles_deg, float))
    M = np.linalg.inv(k * np.c_[np.ones(3), np.cos(2 * a), np.sin(2 * a)])
    im = [np.asarray(x, float) for x in images]
    I = sum(M[0, k] * im[k] for k in range(3))
    Q = sum(M[1, k] * im[k] for k in range(3))
    U = sum(M[2, k] * im[k] for k in range(3))
    return I, np.sqrt(Q * Q + U * U)


# ------------------------------------------------------------------ downloads of calibration files
def _get(url, dest):
    dest = Path(dest)
    if dest.exists():
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as o:
        while chunk := r.read(1 << 20):
            o.write(chunk)
    tmp.replace(dest)
    return dest


def calib_dir(home):
    return Path(home) / "data" / "calib"


def cor1_flat(home, sc):
    return _get(f"{SSW}/stereo/secchi/calibration/20090723_flatfd_fCc1{sc}.fts", calib_dir(home) / f"20090723_flatfd_fCc1{sc}.fts")


def c2_vignetting(home):
    return _get(f"{SSW}/soho/lasco/idl/data/calib/c2vig_final.fts", calib_dir(home) / "c2vig_final.fts")


def cor1_background(home, sc, when, angle, kind="monthly_roll"):
    """Path of the COR1 background of the SECCHI archive nearest in time to `when` for one polariser angle
    (monthly_roll: the monthly minimum combined with the roll minimum; monthly_min), downloaded once."""
    import re
    from datetime import timedelta
    prefix = {"monthly_roll": "mrc1", "monthly_min": "mc1"}[kind]
    sub = {"monthly_roll": "monthly_roll", "monthly_min": "monthly_min"}[kind]
    best = None
    for month in sorted({f"{when:%Y%m}", f"{(when.replace(day=1) - timedelta(days=1)):%Y%m}"}):
        url = f"{SSC_BKG}/{sc.lower()}/{sub}/{month}/"
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                names = set(re.findall(rf'href="({prefix}{sc}_p{angle:03d}_(\d{{6}})\.fts)"', r.read().decode()))
        except OSError:
            continue
        for name, d in names:
            dt = abs((datetime.strptime(d, "%y%m%d") - when).total_seconds())
            if best is None or dt < best[0]:
                best = (dt, url + name, name)
    if best is None:
        raise FileNotFoundError(f"no COR1-{sc} {kind} background near {when:%Y-%m-%d} in {SSC_BKG}")
    return _get(best[1], calib_dir(home) / best[2])


# ------------------------------------------------------------------ observed pB: STEREO/SECCHI COR1
SEBIP_FACTOR = {1: 2.0, 16: 64.0, 17: 64.0, 50: 4.0, 53: 4.0, 82: 2.0, 83: 4.0, 84: 8.0, 85: 16.0, 86: 32.0, 87: 64.0,
                88: 128.0, 118: 3.0}


def sebip_factor(hdr):
    """Factor that undoes the onboard divisions of the SECCHI Electronics Box (scc_sebip): IP_00_19 holds
    twenty three-digit function codes; a division by 2 that the ground system already corrected (DIV2CORR)
    is not counted twice. A square root (code 2) is not supported."""
    ip = str(hdr.get("IP_00_19", "")).rjust(60)
    codes = [int(ip[3 * k:3 * k + 3]) for k in range(20) if ip[3 * k:3 * k + 3].strip()]
    if 117 in codes:
        codes = codes[codes.index(117):]
    if 2 in codes:
        raise ValueError("square-root onboard compression is not supported")
    n1 = codes.count(1) - (1 if str(hdr.get("DIV2CORR", "F")).upper() in ("T", "TRUE") or hdr.get("DIV2CORR") is True else 0)
    f = 2.0 ** max(n1, 0)
    for c in codes:
        if c != 1 and c in SEBIP_FACTOR:
            f *= SEBIP_FACTOR[c]
    return f


def secchi_bias(hdr):
    """Bias [DN] of a SECCHI image as get_biasmean: zero when it was subtracted onboard (code 103) or on the
    ground (OFFSETCR), otherwise BIASMEAN times the number of pixels summed in the SEB."""
    ip = str(hdr.get("IP_00_19", ""))
    if " 103" in " " + ip or float(hdr.get("OFFSETCR", 0) or 0) > 0:
        return 0.0
    ipsum = float(hdr.get("IPSUM", 1) or 1)
    return float(hdr["BIASMEAN"]) * ((2.0 ** (ipsum - 1)) ** 2 if ipsum > 1 else 1.0)


def cor1_calfac(hdr):
    """Calibration factor of COR1 [B_sun per DN/s] with its decline in time, divided by the number of pixels
    summed in the SEB (get_calfac of secchi_prep)."""
    sc = hdr["OBSRVTRY"][-1]
    calfac, t0, rate = {"A": (6.578e-11, "2007-12-01T03:41:48.174", 0.00648),
                        "B": (7.080e-11, "2008-01-17T02:20:15.717", 0.00258)}[sc]
    t = datetime.fromisoformat(str(hdr.get("DATE-AVG") or hdr["DATE-OBS"])[:19])
    years = (t - datetime.fromisoformat(t0[:19])).total_seconds() / (3600 * 24 * 365.25)
    ipsum = float(hdr.get("IPSUM", 1) or 1)
    return calfac / (1 - rate * years) / ((2.0 ** (ipsum - 1)) ** 2 if ipsum > 1 else 1.0)


def _rebin(a, shape):
    sy, sx = a.shape[0] // shape[0], a.shape[1] // shape[1]
    return a[:shape[0] * sy, :shape[1] * sx].reshape(shape[0], sy, shape[1], sx).mean(axis=(1, 3))


def _hpc_axes(m):
    """x, y [R_sun] of every pixel of a sunpy map (helioprojective, solar north up), from its WCS."""
    import astropy.units as u
    ny, nx = m.data.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    c = m.pixel_to_world(xx * u.pix, yy * u.pix)
    rs = m.rsun_obs.to_value(u.arcsec)
    return c.Tx.to_value(u.arcsec) / rs, c.Ty.to_value(u.arcsec) / rs


def cor1_triplet(home, files, background="monthly_roll"):
    """Calibrated tB and pB [B_sun] of one COR1 polarisation triplet with the sky-plane x, y [R_sun] of its
    pixels; returns dict(tb, pb, x, y, time, notes)."""
    import sunpy.map
    from astropy.io import fits
    imgs, angles, notes = [], [], []
    hdr0 = None
    for f in files:
        with fits.open(f) as h:
            hdr, raw = h[0].header, h[0].data.astype(float)
        hdr0 = hdr0 or hdr
        sc = hdr["OBSRVTRY"][-1]
        im = (raw * sebip_factor(hdr) - secchi_bias(hdr)) / float(hdr["EXPTIME"])
        when = datetime.fromisoformat(hdr["DATE-OBS"][:19])
        angle = int(round(float(hdr["POLAR"])))
        if background != "none":
            bf = cor1_background(home, sc, when, angle, background)
            with fits.open(bf) as hb:
                bkg = hb[0].data.astype(float) * 4.0 ** (float(hdr.get("IPSUM", 1)) - float(hb[0].header.get("IPSUM", 1)))
            im = im - (bkg if bkg.shape == im.shape else _rebin(bkg, im.shape))
            notes.append(bf.name)
        with fits.open(cor1_flat(home, sc)) as hf:
            flat = hf[0].data.astype(float)
        flat = _rebin(flat, im.shape)
        with np.errstate(divide="ignore", invalid="ignore"):
            im = np.where(flat > 0.05, im / flat, np.nan) * cor1_calfac(hdr)
        imgs.append(im)
        angles.append(angle)
    tb, pb = stokes_pb(imgs, angles)
    m = sunpy.map.Map(files[0])
    x, y = _hpc_axes(m)
    return dict(tb=tb, pb=pb, x=x, y=y, time=datetime.fromisoformat(hdr0["DATE-OBS"][:19]), notes=notes,
                instrument=f"COR1-{hdr0['OBSRVTRY'][-1]}")


# ------------------------------------------------------------------ observed pB: SOHO/LASCO C2
def c2_calfac(hdr, polarised=True):
    """Calibration factor of C2 [B_sun per DN/s] for the orange filter (c2_calfactor): a linear function of
    the modified Julian date, divided by 0.25256 for a polariser, and divided by the number of summed pixels."""
    from astropy.time import Time
    if str(hdr.get("FILTER", "")).strip().upper() != "ORANGE":
        raise ValueError(f"C2 filter {hdr.get('FILTER')}: only the orange filter is calibrated here")
    mjd = Time(hdr["DATE-OBS"].replace("/", "-") + "T" + hdr["TIME-OBS"]).mjd
    cf = 4.60403e-07 * int(mjd) + 0.0374116
    if polarised:
        cf /= 0.25256
    for k in ("SUMCOL", "SUMROW"):
        if int(hdr.get(k, 0)) > 0:
            cf /= int(hdr[k])
    for k in ("LEBXSUM", "LEBYSUM"):
        if int(hdr.get(k, 1)) > 1:
            cf /= int(hdr[k])
    return cf * 1e-10


def c2_image(home, f):
    """Calibrated C2 image [B_sun] (c2_calibrate: (DN - bias) calfac / exptime * vignetting) and its header;
    the bias of an image summed in the LEB is the header OFFSET times the number of summed pixels."""
    from astropy.io import fits
    with fits.open(f) as h:
        hdr, raw = h[0].header, h[0].data.astype(float)
    with fits.open(c2_vignetting(home)) as hv:
        vig = np.clip(hv[0].data.astype(float), 0.0, 100.0)
    x1, x2, y1, y2 = int(hdr["R1COL"]) - 20, int(hdr["R2COL"]) - 20, int(hdr["R1ROW"]) - 1, int(hdr["R2ROW"]) - 1
    vig = vig[y1:y2 + 1, x1:x2 + 1]
    if vig.shape != raw.shape:
        vig = _rebin(vig, raw.shape)
    pol = str(hdr.get("POLAR", "Clear")).strip().upper() != "CLEAR"
    bias = float(hdr.get("OFFSET", 0.0)) * max(int(hdr.get("LEBXSUM", 1)), 1) * max(int(hdr.get("LEBYSUM", 1)), 1)
    im = (raw - bias) * c2_calfac(hdr, pol) / float(hdr["EXPTIME"]) * vig
    return im, hdr


def c2_axes(hdr, shape):
    """x, y [R_sun] of the pixels of a C2 image, from the world coordinate system that sunpy builds from the header
    (LASCOMap: plate scale, Sun centre and orientation, as for COR1 in `_hpc_axes`), so that the position angles are
    those of the sky with solar north up; the date is the date and time of the exposure."""
    import sunpy.map
    h = {k: v for k, v in dict(hdr).items() if k not in ("COMMENT", "HISTORY", "")}
    t = hdr["DATE-OBS"].replace("/", "-") + "T" + hdr["TIME-OBS"]
    h["DATE-OBS"] = t
    h.pop("TIME-OBS", None)
    m = sunpy.map.Map(np.zeros(shape), h)
    return _hpc_axes(m)


C2_ANGLE = {"+60 DEG": 60.0, "0 DEG": 0.0, "-60 DEG": -60.0}


def c2_triplet(home, files):
    """Calibrated tB and pB [B_sun] of one C2 polarisation sequence (+60, 0, -60 deg), with x, y [R_sun]."""
    from datetime import timedelta
    imgs, angles, hdrs = [], [], []
    for f in files:
        im, hdr = c2_image(home, f)
        imgs.append(im)
        angles.append(C2_ANGLE[str(hdr["POLAR"]).strip().upper()])
        hdrs.append(hdr)
    tb, pb = stokes_pb(imgs, angles, k=1.0)
    x, y = c2_axes(hdrs[0], tb.shape)
    t = datetime.strptime(hdrs[0]["DATE-OBS"] + " " + hdrs[0]["TIME-OBS"][:8], "%Y/%m/%d %H:%M:%S")
    return dict(tb=tb, pb=pb, x=x, y=y, time=t, notes=[Path(f).name for f in files], instrument="C2")


# ------------------------------------------------------------------ the sequences in the folder of the white-light images
def sequences(wl_dir, inst, sc="B"):
    """[(time, [three files])] of the polarisation sequences of COR1 (sc) or C2 in a folder, earliest first."""
    from astropy.io import fits
    d = Path(wl_dir)
    out = []
    if inst == "cor1":
        from .whitelight import _frames
        return _frames(d, "cor1", sc)
    cur = []
    for f in sorted(d.glob("2*.fts")):
        try:
            h = fits.getheader(f)
        except Exception:
            continue
        if h.get("DETECTOR") != "C2" or str(h.get("POLAR", "")).strip().upper() not in C2_ANGLE:
            continue
        t = datetime.strptime(h["DATE-OBS"] + " " + h["TIME-OBS"][:8], "%Y/%m/%d %H:%M:%S")
        if cur and ((t - cur[0][0]).total_seconds() > 1800 or str(h["POLAR"]).strip().upper() in
                    [str(x[2]).strip().upper() for x in cur]):
            cur = []
        cur.append((t, f, h["POLAR"]))
        if len(cur) == 3:
            out.append((cur[0][0], [c[1] for c in cur]))
            cur = []
    return out


def fetch_c2_polarised(home, day, max_mib=1.5, log=print):
    """Download the LASCO C2 images of one day that are smaller than the full-frame total-brightness images
    (max_mib), which holds the summed 512 x 512 images of the polarisation sequences, from the Virtual Solar
    Observatory into data/wl, one at a time; returns the number of files that failed."""
    import astropy.units as u
    from sunpy.net import Fido, attrs as a
    from .paths import wl_dir
    r = Fido.search(a.Time(f"{day:%Y-%m-%d} 00:00", f"{day:%Y-%m-%d} 23:59:59"), a.Instrument("LASCO"), a.Detector("C2"))
    t = r[0]
    small = t[np.array([q.to_value(u.Mibyte) < max_mib for q in t["Size"]])]
    log(f"pb: {len(small)} small C2 files on {day:%Y-%m-%d} (candidate polarisation sequences)")
    if not len(small):
        return 0
    res = Fido.fetch(small, path=str(wl_dir(home) / "{file}"), max_conn=1, progress=False)
    return len(res.errors)


# ------------------------------------------------------------------ comparison and fit
def polar_bins(obs, r_range, r_step, pa_step):
    """Median of the observed pB in bins of distance and position angle (counter-clockwise from north);
    returns arrays r, pa [deg], pb, n of the bins with valid pixels."""
    r = np.hypot(obs["x"], obs["y"])
    pa = np.degrees(np.arctan2(-obs["x"], obs["y"])) % 360.0
    v = obs["pb"]
    ok = np.isfinite(v) & (v > 0) & (r >= r_range[0]) & (r < r_range[1])
    rb = np.arange(r_range[0], r_range[1] + 1e-9, r_step)
    pb_ = np.arange(0.0, 360.0 + 1e-9, pa_step)
    ir = np.digitize(r[ok], rb) - 1
    ip = np.digitize(pa[ok], pb_) - 1
    out = []
    vals = v[ok]
    key = ir * len(pb_) + ip
    order = np.argsort(key)
    key, vals = key[order], vals[order]
    edges = np.flatnonzero(np.diff(key)) + 1
    for k, grp in zip(np.split(key, edges), np.split(vals, edges)):
        if not len(grp) or len(grp) < 3:
            continue
        i, j = divmod(int(k[0]), len(pb_))
        if 0 <= i < len(rb) - 1 and 0 <= j < len(pb_) - 1:
            out.append((rb[i] + r_step / 2, pb_[j] + pa_step / 2, float(np.median(grp)), len(grp)))
    a = np.array(out, float).reshape(-1, 4)
    return dict(r=a[:, 0], pa=a[:, 1], pb=a[:, 2], n=a[:, 3].astype(int))


def fit_scale(ratio, pa, r=None, ranges=None):
    """Density scale from the ratios observed/model pB of the bins: the median of ln(ratio) over the bins
    (exp), with the standard error from the scatter of the per-position-angle medians (1.4826 MAD / sqrt n)
    and their 16 to 84 percent range. With `ranges` ([[r0, r1], ...]) the same for each radial range."""
    lr = np.log(np.asarray(ratio, float))
    ok = np.isfinite(lr)

    def one(sel):
        sel = sel & ok
        if sel.sum() < 3:
            return dict(n=int(sel.sum()))
        per = np.array([np.median(lr[sel & (pa == a)]) for a in np.unique(pa[sel])])
        mad = 1.4826 * np.median(np.abs(per - np.median(per)))
        s = float(np.exp(np.median(lr[sel])))
        return dict(scale=s, err=float(s * mad / np.sqrt(len(per))), p16=float(np.exp(np.percentile(per, 16))),
                    p84=float(np.exp(np.percentile(per, 84))), n=int(sel.sum()), n_pa=int(len(per)))
    out = one(np.ones(len(lr), bool))
    if ranges is not None and r is not None:
        out["ranges"] = [dict(r0=a, r1=b, **one((r >= a) & (r < b))) for a, b in ranges]
    return out


def model_sampler(home, p):
    """ne_at(points) of the model density without scaling (density_scale 1), with model_dlon."""
    from . import mast_model
    from .heights import ne_sampler
    from .paths import model_dir
    m = mast_model.load(model_dir(home, p), density_scale=1.0, va_density="model")
    return ne_sampler(m, dict(p, density_scale=1.0, va_density="model"))


def analyse(home, p=None, log=print):
    """Observed and model pB of every instrument of pb.instruments at the chosen time, the fitted scales,
    and the tables; returns (summary, rows)."""
    from .paths import config, merge, wl_dir
    from .whitelight import utc
    p = merge(config(home), p)
    s = settings(p)
    t_win = (utc(p["wl_window"][0]), utc(p["wl_window"][1]))
    when = utc(s["time"]) if s["time"] else t_win[0]
    ne_at = model_sampler(home, p)
    rows, per_inst = [], {}
    for inst in s["instruments"]:
        if inst == "c2" and not sequences(wl_dir(home), inst) and s["fetch"]:
            try:
                fetch_c2_polarised(home, when.date(), log=log)
            except Exception as e:
                log(f"pb: C2 polarisation sequences not fetched ({e})")
        seqs = [q for q in sequences(wl_dir(home), inst, p["stereo_sc"]) if t_win[0] <= q[0] <= t_win[1]] or \
            sequences(wl_dir(home), inst, p["stereo_sc"])
        if not seqs:
            log(f"pb: no {inst.upper()} polarisation sequence in {wl_dir(home)}")
            per_inst[inst] = dict(n=0, reason="no polarisation sequence in the white-light folder")
            continue
        before = [q for q in seqs if q[0] <= utc(p["time"])] if not s["time"] else []
        t, files = min(before or seqs, key=lambda q: abs((q[0] - when).total_seconds()))   # the corona before the CME
        try:
            obs = cor1_triplet(home, files, s["background"]) if inst == "cor1" else c2_triplet(home, files)
        except Exception as e:
            log(f"pb: {inst} failed ({e})")
            per_inst[inst] = dict(n=0, reason=str(e))
            continue
        b = polar_bins(obs, s["r_range"][inst], s["r_step"], s["pa_step"])
        observer = p["stereo"] if inst == "cor1" else p["earth"]
        th = np.radians(b["pa"])
        mp, mt = thomson(ne_at, -b["r"] * np.sin(th), b["r"] * np.cos(th), observer, s["u"], s["los_step"], s["los_r_out"])
        ratio = b["pb"] / mp
        fit = fit_scale(ratio, b["pa"], b["r"], s["radial_ranges"])
        fit.update(instrument=obs["instrument"], time=obs["time"].isoformat(), files=[Path(f).name for f in files],
                   background=obs["notes"] if inst == "cor1" else "none (pB)")
        per_inst[inst] = fit
        for r_, pa_, o_, m_, tm_ in zip(b["r"], b["pa"], b["pb"], mp, mt):
            rows.append(dict(instrument=obs["instrument"], r_rsun=r_, pa_deg=pa_, pb_obs=o_, pb_model=m_, tb_model=tm_,
                             ratio=o_ / m_))
        log(f"pb: {obs['instrument']} {obs['time']:%H:%M:%S}: density scale {fit.get('scale', np.nan):.2f} "
            f"+/- {fit.get('err', np.nan):.2f} ({len(b['r'])} bins)")
    if rows:
        allr = np.array([q["ratio"] for q in rows])
        pa = np.array([q["pa_deg"] for q in rows])
        rr = np.array([q["r_rsun"] for q in rows])
        comb = fit_scale(allr, pa, rr, s["radial_ranges"])
    else:
        comb = dict(n=0)
    summ = dict(time=when.isoformat(), u=s["u"], background=s["background"], combined=comb, instruments=per_inst,
                current_density_scale=p["density_scale"], model=p["run"])
    return summ, rows


def figure(rows, summ, out):
    """Radial profiles of the observed pB and of the scaled model pB, median over position angle, per instrument,
    and the ratio against position angle."""
    from matplotlib.figure import Figure
    fig = Figure(figsize=(11, 4.4))
    a1, a2 = fig.add_subplot(1, 2, 1), fig.add_subplot(1, 2, 2)
    s = summ["combined"].get("scale", np.nan)
    for k, inst in enumerate(sorted({q["instrument"] for q in rows})):
        R = [q for q in rows if q["instrument"] == inst]
        r = np.array([q["r_rsun"] for q in R])
        rs = np.unique(r)
        med = lambda key: np.array([np.median([q[key] for q in R if q["r_rsun"] == x]) for x in rs])
        col = f"C{k}"
        a1.semilogy(rs, med("pb_obs"), "o", color=col, label=f"{inst} observed")
        a1.semilogy(rs, med("pb_model"), ":", color=col, label=f"{inst} model x1")
        a1.semilogy(rs, s * med("pb_model"), "-", color=col, label=f"{inst} model x{s:.2f}")
        pa = np.array([q["pa_deg"] for q in R])
        a2.plot(pa, [q["ratio"] for q in R], ".", color=col, ms=3, label=inst)
    a1.set_xlabel("distance from Sun centre (R_sun)")
    a1.set_ylabel("pB (mean solar brightness)")
    a1.legend(fontsize=7)
    a2.axhline(s, color="k", lw=1)
    a2.set_yscale("log")
    a2.set_xlabel("position angle (deg, from north through east)")
    a2.set_ylabel("observed / model pB")
    a2.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    return out


def write(out, summ, rows):
    import csv
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    if rows:
        with open(out / "pb_profiles.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        figure(rows, summ, out / "pb_profiles.png")
    (out / "pb_fit.json").write_text(json.dumps(summ, indent=1, default=str) + "\n", encoding="utf-8")
    return out


def store(home, summ):
    """Write the fitted scale into config.json as `pb_fit`, which `density_scale: "pb"` uses."""
    c = summ["combined"]
    if "scale" not in c:
        return None
    f = Path(home) / "config.json"
    user = json.loads(f.read_text(encoding="utf-8"))
    user["pb_fit"] = dict(scale=round(c["scale"], 4), err=round(c["err"], 4), p16=round(c["p16"], 4), p84=round(c["p84"], 4),
                          time=summ["time"], instruments={k: round(v["scale"], 4) for k, v in summ["instruments"].items()
                                                          if "scale" in v})
    f.write_text(json.dumps(user, indent=2) + "\n", encoding="utf-8")
    return user["pb_fit"]


def resolve(p):
    """Copy of the configuration with `density_scale` "pb" replaced by the scale of `pb_fit`."""
    if p.get("density_scale_pending_pb"):                  # set by paths.config when pb_fit was missing
        raise ValueError('density_scale "pb" needs the fit of `corona3d pb` (key pb_fit of config.json); run `corona3d pb` first')
    if p.get("density_scale") != "pb":
        return p
    fit = p.get("pb_fit") or {}
    if "scale" not in fit:
        raise ValueError('density_scale "pb" needs the fit of `corona3d pb` (key pb_fit of config.json); run it first')
    return dict(p, density_scale=float(fit["scale"]), density_scale_source="pb")


def report_lines(summ):
    c = summ["combined"]
    L = []
    if "scale" not in c:
        return ["- No polarised-brightness sequence could be compared with the model; density_scale is unchanged."]
    L.append(f"- The polarised brightness of the corona, observed at {summ['time'][:16]} UT and forward-modelled from the "
             f"MAS density by Thomson scattering (limb darkening u = {summ['u']:g}), gives a density scale of "
             f"{c['scale']:.2f} +/- {c['err']:.2f} (16 to 84 percent of the position angles {c['p16']:.2f} to {c['p84']:.2f}), "
             f"whereas the configuration uses {summ['current_density_scale']}.")
    for k, v in summ["instruments"].items():
        if "scale" in v:
            rr = "; ".join(f"{q['r0']:g} to {q['r1']:g} R_sun: {q['scale']:.2f}" for q in v.get("ranges", []) if "scale" in q)
            L.append(f"- {v['instrument']} at {v['time'][11:19]} UT: {v['scale']:.2f} +/- {v['err']:.2f} ({v['n']} bins; {rr}).")
        else:
            L.append(f"- {k.upper()}: not used ({v.get('reason', 'no data')}).")
    L.append('- With `"density_scale": "pb"` in config.json the reconstruction uses the fitted scale.')
    return L


def run(home, cfg, a):
    """The `corona3d pb` command."""
    from .report_sections import update
    summ, rows = analyse(home, dict(pb=dict(cfg.get("pb") or {}, **({"time": a.time} if a.time else {}))))
    out = write(Path(home) / a.out, summ, rows)
    if not a.no_store:
        store(home, summ)
    lines = report_lines(summ)
    print("\n".join(lines))
    update(home, "Density from the polarised brightness", lines)
    print(f"wrote {out}")
