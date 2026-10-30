"""Dynamic spectra of the type II burst, its lanes, the shock speed from their drift, and the Alfven Mach number
from their band splitting (`corona3d spectra`).

Sources of the spectra (public, no login):
    e-CALLISTO     the FITS files of the stations of the network (Benz et al. 2005), one file of 15 min per station and
                   focus code, from the archive of FHNW; a station is used when its files cover the burst interval
    OVRO-LWA       the daily total-power spectrum (Stokes I and V, 15 to 85 MHz, 0.29 s) of the public bucket of the
                   OVRO-LWA solar data, read with HTTP range requests for the burst interval only

Lanes. Each spectrum is reduced to 1 s and to a logarithmic frequency grid, the median of every channel over the
interval is subtracted (in dB), channels dominated by interference are removed, and the result is normalised by its
robust scatter. A lane of a type II burst is close to a straight line in (t, ln f), since the plasma frequency falls
roughly exponentially with height for a shock of constant speed. All lines ln f = ln f_ref + D (t - t_ref) with
drift rates D from `drift_range` are scored by the mean normalised intensity along them, and the best line is refined
by the frequency of the brightest pixel within 8 percent of it at each time, to which a straight line is fitted again.
A second lane parallel to it is searched at frequency ratios of 1.08 to 1.45 (band splitting, the upper and lower
branch of one lane) and of 1.7 to 2.3 (the harmonic of a fundamental lane). The main lane is accepted when its score
exceeds the median of the scores of all candidate lines by `min_z` times their scatter, and a second lane when its mean
over the stretch of the main lane exceeds the median over the scanned ratios by `min_z` times their scatter and reaches
a quarter of the score of the main lane. Lanes may also be given by hand
(`spectra.lanes`) as two points each, which are then refined in the same way.

Shock speed. With n_e(r) the density model and f = f_p (fundamental) or 2 f_p (harmonic), the drift gives
    V = 2 (d ln f / dt) / (d ln n_e / dr),
evaluated at the height of the frequency of the lane at the middle of its interval, for the models newkirk1 to
newkirk6, saito, and mas (the scaled MAS density along the axis of the CME shell).

Band splitting. If the upper and lower branches are emitted ahead of and behind the shock, the density compression is
X = (f_U / f_L)^2, and for a perpendicular shock in a plasma of low beta (Vrsnak et al. 2002)
    M_A = sqrt( X (X + 5) / (2 (4 - X)) ).
"""
import json
import re
import urllib.request
import warnings
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

CALLISTO = "http://soleil80.cs.technik.fhnw.ch/solarradio/data/2002-20yy_Callisto/{d:%Y/%m/%d}/"
LWA_SPEC = "https://ovro-lwa-solar.s3.us-west-2.amazonaws.com/spec_fits/{d:%Y}/{d:%Y%m%d}.fits"
CNAME = re.compile(r'href="(([A-Za-z0-9\-]+)_(\d{8})_(\d{6})_(\d{2})\.fit\.gz)"')
PREFERRED = ["BIR", "BLEN7M", "BLENSW", "GLASGOW", "HUMAIN", "DARO", "TRIEST", "HURBANOVO", "MRT2", "MRT3",
             "ALASKA-COHOE", "ALASKA-HAARP", "ALASKA-ANCHORAGE", "USA-ARIZONA-ERAU", "ROSWELL-NM", "MEXART",
             "GERMANY-DLR", "AUSTRIA-UNIGRAZ", "SPAIN-PERALEJOS", "SWISS-Landschlacht", "DENMARK", "USA-BOSTON"]
DEFAULTS = dict(
    sources=["callisto", "ovro_lwa"],
    stations=None,                 # e-CALLISTO stations (names or name_focus); default: PREFERRED then the others
    max_files=8,                   # e-CALLISTO files downloaded at most
    window=None,                   # [start, end] of the type II burst; default: SWPC list, else the radio interval
    pad_min=2.0,                   # minutes added before and after the window
    f_range=[20.0, 450.0],         # MHz, frequencies searched for lanes
    drift_range=[-0.01, -0.0008],  # 1/s, d ln f / dt of the lane lines (type II drift rates of the metric range)
    background_s=300.0,            # s, running median of every channel that is subtracted (radio-frequency interference, continua)
    lane_min_s=120.0,              # s, length of the stretch over which a lane is scored
    min_z=4.0,                     # acceptance of a lane: score above the median of all lines by min_z times their scatter
    lanes=None,                    # [{"t0": iso, "f0": MHz, "t1": iso, "f1": MHz, "spectrum": name}] given by hand
    emission="auto",               # the strongest lane: "F", "H", or "auto" (H unless a lane at half its frequency is found)
    density_model="mas",           # model of the reported speed; every model is listed
    models=["newkirk1", "newkirk2", "newkirk4", "newkirk6", "saito", "mas"],
)


def settings(p):
    return dict(DEFAULTS, **(p.get("spectra") or {}))


def _get(url, dest, headers=None):
    dest = Path(dest)
    if dest.exists():
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers=headers or {})
    with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "wb") as o:
        while chunk := r.read(1 << 20):
            o.write(chunk)
    tmp.replace(dest)
    return dest


def spectra_dir(home):
    return Path(home) / "data" / "spectra"


# ------------------------------------------------------------------ e-CALLISTO
def callisto_list(day, opener=urllib.request.urlopen):
    """[(station, focus, start datetime, file name)] of the e-CALLISTO archive for one day."""
    with opener(CALLISTO.format(d=day), timeout=60) as r:
        text = r.read().decode("utf-8", errors="replace")
    out = []
    for name, st, d, t, fc in set(CNAME.findall(text)):
        out.append((st, fc, datetime.strptime(d + t, "%Y%m%d%H%M%S"), name))
    return sorted(out, key=lambda x: (x[0], x[1], x[2]))


def callisto_choose(listing, t0, t1, stations=None, max_files=8):
    """File names of the stations whose 15 min files cover [t0, t1], at most max_files, grouped by station and focus,
    in the order of `stations` (or PREFERRED) and then by name."""
    groups = {}
    for st, fc, ts, name in listing:
        if ts <= t1 and ts + timedelta(minutes=15) >= t0:
            groups.setdefault((st, fc), []).append((ts, name))
    order = stations or PREFERRED

    def rank(k):
        for i, s in enumerate(order):
            if s in (k[0], f"{k[0]}_{k[1]}"):
                return (i, k)
        return (len(order), k)
    out = []
    for k in sorted(groups, key=rank):
        if stations and rank(k)[0] == len(order):
            continue
        names = [n for _, n in sorted(groups[k])]
        if len(out) + len(names) > max_files:
            break
        out.append((f"{k[0]}_{k[1]}", names))
    return out


def read_callisto(paths):
    """Spectrum of one station and focus from its consecutive files: dict(name, time [datetime], freq [MHz], data
    (n_freq, n_time) [digits]), channels sorted by frequency and duplicates averaged."""
    from astropy.io import fits
    parts = []
    for p in paths:
        with fits.open(p) as h:
            hd, d = h[0].header, h[0].data.astype(float)
            tt = np.asarray(h[1].data["TIME"][0], float)
            ff = np.asarray(h[1].data["FREQUENCY"][0], float)
        t_start = datetime.strptime(hd["DATE-OBS"] + " " + hd["TIME-OBS"][:8], "%Y/%m/%d %H:%M:%S")
        parts.append(([t_start + timedelta(seconds=float(s)) for s in tt], ff, d))
    ff = parts[0][1]
    times = sum((p[0] for p in parts), [])
    data = np.concatenate([p[2] if np.array_equal(p[1], ff) else np.full((len(ff), len(p[0])), np.nan) for p in parts], axis=1)
    u, inv = np.unique(np.round(ff, 3), return_inverse=True)
    out = np.zeros((len(u), data.shape[1]))
    np.add.at(out, inv, data)
    out /= np.bincount(inv)[:, None]
    return dict(time=times, freq=u, data=out, units="digits")


def fetch_callisto(home, t0, t1, stations=None, max_files=8, log=print):
    """Download the chosen e-CALLISTO files (sequentially) into data/spectra/callisto; returns
    [(name, [paths])]."""
    d = spectra_dir(home) / "callisto"
    chosen = callisto_choose(callisto_list(t0.date()), t0, t1, stations, max_files)
    out = []
    for name, files in chosen:
        paths = []
        for f in files:
            try:
                paths.append(_get(CALLISTO.format(d=t0.date()) + f, d / f))
            except OSError as e:
                log(f"spectra: {f} failed ({e})")
        if paths:
            out.append((name, paths))
    log(f"spectra: e-CALLISTO {', '.join(n for n, _ in out) or 'none'}")
    return out


# ------------------------------------------------------------------ OVRO-LWA total-power spectrum
def _range(url, a, b):
    req = urllib.request.Request(url, headers={"Range": f"bytes={a}-{b}"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read()


def ovro_spectrum(home, t0, t1, ch_step=2, avg_s=1.0, log=print):
    """Stokes I and V [Jy] of the OVRO-LWA daily spectrum between t0 and t1, averaged to avg_s seconds and every
    ch_step-th channel, read by range requests and kept as data/spectra/ovro_lwa_<start>_<end>.npz."""
    from astropy.io import fits
    import io
    cache = spectra_dir(home) / f"ovro_lwa_{t0:%Y%m%dT%H%M%S}_{t1:%H%M%S}.npz"
    if cache.exists():
        z = np.load(cache, allow_pickle=True)
        return dict(time=list(z["time"]), freq=z["freq"], data=z["I"], V=z["V"], units="Jy", name="OVRO-LWA")
    url = LWA_SPEC.format(d=t0)
    raw = b""                                              # header blocks of 2880 bytes up to the one with the END card
    while True:
        raw += _range(url, len(raw), len(raw) + 2879)
        cards = [raw[i:i + 80] for i in range(len(raw) - 2880, len(raw), 80)]
        if any(c[:8].rstrip() == b"END" for c in cards) or len(raw) > 2880 * 50:
            break
    head = fits.Header.fromstring(raw.decode("ascii", errors="replace"))
    nt, nf, npol = int(head["NAXIS1"]), int(head["NAXIS2"]), int(head["NAXIS4"])
    n3 = int(head.get("NAXIS3", 1))
    bitpix = int(head["BITPIX"])
    dt = {-64: ">f8", -32: ">f4", 8: "u1", 16: ">i2", 32: ">i4", 64: ">i8"}[bitpix]
    isz = abs(bitpix) // 8
    bscale, bzero = float(head.get("BSCALE", 1.0)), float(head.get("BZERO", 0.0))
    hlen = len(raw)
    with urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=60) as r:
        size = int(r.headers["Content-Length"])
    data_end = hlen + int(np.ceil(npol * n3 * nf * nt * isz / 2880)) * 2880
    tail = _range(url, data_end, size - 1)
    ext = fits.open(io.BytesIO(fits.PrimaryHDU().header.tostring().encode() + tail))
    freq = np.asarray(ext["SFREQ"].data["sfreq"], float) * 1e3
    ut = ext["UT"].data
    times = np.array([datetime(1858, 11, 17) + timedelta(days=int(m), milliseconds=int(ms)) for m, ms in zip(ut["mjd"], ut["time"])])
    i0, i1 = int(np.searchsorted(times, t0)), int(np.searchsorted(times, t1))
    if i1 <= i0:
        raise ValueError(f"the OVRO-LWA spectrum of {t0:%Y-%m-%d} covers {times[0]:%H:%M}-{times[-1]:%H:%M} UT, not the burst")
    chans = np.arange(0, nf, ch_step)
    nb = max(1, int(round(avg_s / np.median(np.diff([t.timestamp() for t in times[i0:i0 + 50]])))))
    n = (i1 - i0) // nb * nb
    out = np.empty((2, len(chans), n // nb))
    log(f"spectra: OVRO-LWA {len(chans)} channels x {n} samples by range requests")
    for p in range(min(npol, 2)):
        for j, c in enumerate(chans):
            off = hlen + (((p * n3) * nf + c) * nt + i0) * isz
            row = np.frombuffer(_range(url, off, off + n * isz - 1), dt).astype(float) * bscale + bzero
            out[p, j] = row.reshape(-1, nb).mean(axis=1)
    tt = times[i0:i0 + n:nb]
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, time=np.array(tt, dtype=object), freq=freq[chans], I=out[0], V=out[1])
    return dict(time=list(tt), freq=freq[chans], data=out[0], V=out[1], units="Jy", name="OVRO-LWA")


# ------------------------------------------------------------------ cleaning and lanes
def regrid(spec, t0, t1, dt=1.0, n_f=200, f_range=(20.0, 450.0)):
    """Spectrum on a time grid of dt seconds from t0 to t1 and a logarithmic frequency grid inside its band and
    f_range; returns (seconds since t0, frequencies [MHz], data (n_f, n_t), V or None)."""
    ts = np.array([(t - t0).total_seconds() for t in spec["time"]])
    tg = np.arange(0.0, (t1 - t0).total_seconds() + dt / 2, dt)
    f = np.asarray(spec["freq"], float)
    lo, hi = max(f.min(), f_range[0]), min(f.max(), f_range[1])
    fg = np.geomspace(lo, hi, n_f)
    order = np.argsort(f)
    out = []
    for key in ("data", "V"):
        if spec.get(key) is None:
            out.append(None)
            continue
        d = np.asarray(spec[key], float)[order]
        bins = np.clip(np.searchsorted(tg, ts) , 0, len(tg) - 1)
        acc = np.zeros((d.shape[0], len(tg)))
        cnt = np.zeros(len(tg))
        ok = (ts >= 0) & (ts <= tg[-1] + dt)
        np.add.at(acc.T, bins[ok], d[:, ok].T)
        np.add.at(cnt, bins[ok], 1)
        with np.errstate(invalid="ignore"):
            acc /= np.where(cnt > 0, cnt, np.nan)
        g = np.array([np.interp(np.log(fg), np.log(f[order]), acc[:, k]) for k in range(len(tg))]).T
        out.append(g)
    return tg, fg, out[0], out[1]


def normalise(d, units):
    """Background-subtracted, normalised spectrum: dB (digits are already logarithmic) minus the median of each
    channel, channels whose scatter is far above the others set to NaN, divided by the robust scatter."""
    x = np.asarray(d, float)
    if units != "digits":
        x = 10 * np.log10(np.clip(x - np.nanpercentile(x, 1, axis=1, keepdims=True) + 1e-3 * np.nanmax(np.abs(x)), 1e-30, None))
    x = x - np.nanmedian(x, axis=1, keepdims=True)
    s = 1.4826 * np.nanmedian(np.abs(x), axis=1)
    bad = s > 3 * np.nanmedian(s)
    x[bad] = np.nan
    return x / (1.4826 * np.nanmedian(np.abs(x[np.isfinite(x)])) + 1e-30)


def _along(z, tg, fg, lnf_ref, drift, t_ref):
    """Values of z along the lines ln f = lnf_ref + drift (t - t_ref) (lnf_ref an array): (n_lines, n_t), NaN outside
    the band."""
    lf = np.log(fg)
    step = lf[1] - lf[0]
    line = np.asarray(lnf_ref, float)[:, None] + drift * (tg[None, :] - t_ref)
    idx = np.rint((line - lf[0]) / step).astype(int)
    inside = (idx >= 0) & (idx < len(fg))
    vals = np.where(inside, z[np.clip(idx, 0, len(fg) - 1), np.arange(len(tg))[None, :]], np.nan)
    return np.clip(vals, -2, 10)


def _running(vals, w):
    """Running mean over w samples along the last axis, NaN where fewer than w/2 samples are valid."""
    ok = np.isfinite(vals)
    v = np.where(ok, vals, 0.0)
    k = np.ones(w)
    num = np.apply_along_axis(lambda a: np.convolve(a, k, "same"), -1, v)
    den = np.apply_along_axis(lambda a: np.convolve(a, k, "same"), -1, ok.astype(float))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den >= w / 2, num / den, np.nan)


def _score(z, tg, fg, lnf_ref, drift, t_ref, w):
    """(best running mean of z along each line over w samples, the time index of that best)."""
    r = _running(_along(z, tg, fg, lnf_ref, drift, t_ref), w)
    good = np.isfinite(r).any(axis=1)
    best = np.full(r.shape[0], np.nan)
    at = np.zeros(r.shape[0], int)
    best[good] = np.nanmax(r[good], axis=1)
    at[good] = np.nanargmax(r[good], axis=1)
    return best, at


def search(z, tg, fg, drift_range, w, n_drift=50):
    """Best line over the drift rates and reference frequencies; returns (score, ln f_ref, drift, time index of the
    best stretch, all scores)."""
    t_ref = tg[len(tg) // 2]
    lnf = np.log(fg)
    allsc, best = [], (-np.inf, None, None, 0)
    for D in np.linspace(drift_range[0], drift_range[1], n_drift):
        s, at = _score(z, tg, fg, lnf, D, t_ref, w)
        allsc.append(s)
        if np.isfinite(s).any():
            k = int(np.nanargmax(s))
            if s[k] > best[0]:
                best = (float(s[k]), float(lnf[k]), float(D), int(at[k]))
    return best[0], best[1], best[2], best[3], np.concatenate(allsc)


def refine(z, tg, fg, lnf_ref, drift, t_span=None, width=0.08, z_min=2.0):
    """Frequencies of the brightest pixels within `width` (relative) of the line, at the times of t_span (indices)
    where they exceed z_min, and a straight line through them in (t, ln f) with 3-sigma clipping; returns dict(t, f,
    lnf_ref, drift, drift_err, t_ref) or None."""
    t_ref = tg[len(tg) // 2]
    lf = np.log(fg)
    pts = []
    ks = range(len(tg)) if t_span is None else range(max(t_span[0], 0), min(t_span[1], len(tg)))
    for k in ks:
        t = tg[k]
        c = lnf_ref + drift * (t - t_ref)
        sel = np.abs(lf - c) <= width
        if sel.sum() < 2 or not np.isfinite(z[sel, k]).any():
            continue
        col = np.where(sel & np.isfinite(z[:, k]), z[:, k], -np.inf)
        j = int(np.argmax(col))
        if col[j] >= z_min:
            pts.append((t, lf[j]))
    if len(pts) < 5:
        return None
    t, y = np.array(pts).T
    keep = np.ones(len(t), bool)
    for _ in range(5):
        coef = np.polyfit(t[keep] - t_ref, y[keep], 1)
        res = y - np.polyval(coef, t - t_ref)
        sd = np.std(res[keep])
        new = np.abs(res) <= 3 * sd + 1e-12
        if (new == keep).all() or new.sum() < 5:
            break
        keep = new
    if keep.sum() < 5:
        return None
    coef, cov = np.polyfit(t[keep] - t_ref, y[keep], 1, cov=True)
    return dict(t=t[keep], f=np.exp(y[keep]), lnf_ref=float(coef[1]), drift=float(coef[0]),
                drift_err=float(np.sqrt(cov[0, 0])), t_ref=float(t_ref))


def given_lanes(z, tg, fg, t0, given):
    """Lanes given by hand, each as two points (t0, f0) and (t1, f1) read off the spectrum, with an optional kind
    (main, split_upper, split_lower, harmonic, fundamental); every lane is refined by the brightest pixels within
    8 percent of its line between t0 and t1. Returns (lanes, statistics) as `lanes`."""
    from .whitelight import utc
    out = []
    t_ref = tg[len(tg) // 2]
    for g in given:
        a, b = (utc(g["t0"]) - t0).total_seconds(), (utc(g["t1"]) - t0).total_seconds()
        D = (np.log(g["f1"]) - np.log(g["f0"])) / (b - a)
        lnf = np.log(g["f0"]) + D * (t_ref - a)
        span = (int(np.searchsorted(tg, a)), int(np.searchsorted(tg, b)) + 1)
        fit = refine(z, tg, fg, lnf, D, span)
        vals = _along(z, tg, fg, [lnf], D, t_ref)[0, span[0]:span[1]]
        out.append(dict(kind=g.get("kind", "main"), ratio=float(g.get("ratio", 1.0)), score=float(np.nanmean(vals)),
                        z=np.nan, fit=fit, given=True, emission=g.get("emission")))
    mains = [l for l in out if l["kind"] == "main"]
    for l in out:                                     # ratios of the partners from the fitted lines
        if mains and l["kind"] != "main" and l["fit"] and mains[0]["fit"]:
            m = mains[0]["fit"]
            l["ratio"] = float(np.exp(l["fit"]["lnf_ref"] - m["lnf_ref"]))
    return out, dict(given=len(out))


def clean(z, tg, s):
    """Subtract the running median of every channel over background_s and then the median over frequency at
    every time."""
    from scipy.ndimage import median_filter
    dt = float(np.median(np.diff(tg)))
    nan = ~np.isfinite(z)
    z = np.where(nan, 0.0, z)
    z = z - median_filter(z, size=(1, max(3, int(round(s["background_s"] / dt)))), mode="nearest")
    with np.errstate(invalid="ignore"):
        z = z - np.median(z, axis=0, keepdims=True)
    z[nan] = np.nan
    return z


def lanes(z, tg, fg, s):
    """The main lane and its partners: band split (ratio 1.08 to 1.45 above or below) and harmonic (1.7 to 2.3).
    The running median of every channel over background_s (interference and continua that last longer than a lane
    stays in one channel) and the intensity common to all frequencies at a time (broadband bursts) are removed first. Returns (list of lanes
    dict(kind, ratio, score, z, fit), statistics of the scores)."""
    dt = float(np.median(np.diff(tg)))
    w = max(5, int(round(s["lane_min_s"] / dt)))
    sc, lnf, D, at, allsc = search(z, tg, fg, s["drift_range"], w)
    med = float(np.nanmedian(allsc))
    mad = float(1.4826 * np.nanmedian(np.abs(allsc - med)))
    zsc = (sc - med) / (mad + 1e-30)
    out = []
    if not np.isfinite(zsc) or zsc < s["min_z"]:
        return out, dict(median=med, scatter=mad, best=sc, z_best=float(zsc))
    t_ref = tg[len(tg) // 2]
    # the stretch of the lane: where the running mean along the best line stays above half its peak
    run = _running(_along(z, tg, fg, [lnf], D, t_ref), w)[0]
    above = np.flatnonzero(np.nan_to_num(run, nan=-np.inf) >= 0.5 * sc)
    seg = np.split(above, np.flatnonzero(np.diff(above) > 1) + 1)
    span = next((g for g in seg if g[0] <= at <= g[-1]), np.array([at]))
    span = (int(span[0] - w // 2), int(span[-1] + w // 2 + 1))
    main = refine(z, tg, fg, lnf, D, span)
    if main is None or not 2 * s["drift_range"][0] <= main["drift"] <= 0.5 * s["drift_range"][1]:
        return out, dict(median=med, scatter=mad, best=sc, z_best=float(zsc), rejected="refined drift outside the range")
    out.append(dict(kind="main", ratio=1.0, score=sc, z=zsc, fit=main))
    for kind, rng in (("split_upper", (1.08, 1.45)), ("split_lower", (1 / 1.45, 1 / 1.08)),
                      ("harmonic", (1.7, 2.3)), ("fundamental", (1 / 2.3, 1 / 1.7))):
        ratios = np.geomspace(*rng, 40)
        lnr = main["lnf_ref"] + np.log(ratios)
        vals = _along(z, tg, fg, lnr, main["drift"], t_ref)[:, span[0]:span[1]] if span[1] > span[0] else None
        if vals is None or not np.isfinite(vals).any():
            continue
        with np.errstate(invalid="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            s_ = np.where(np.isfinite(vals).sum(axis=1) >= max(3, 0.5 * vals.shape[1]), np.nanmean(vals, axis=1), np.nan)
        if not np.isfinite(s_).any():
            continue
        k = int(np.nanargmax(s_))
        bg = float(np.nanmedian(s_))
        sd = float(1.4826 * np.nanmedian(np.abs(s_ - bg)))
        zz = (s_[k] - bg) / (sd + 1e-30)                   # against the other ratios of the scan
        if zz >= s["min_z"] and s_[k] >= 0.25 * sc and 0 < k < len(ratios) - 1:   # a peak inside the scan, not its edge
            fit = refine(z, tg, fg, float(lnr[k]), main["drift"], span, width=0.04)
            out.append(dict(kind=kind, ratio=float(ratios[k]), score=float(s_[k]), z=float(zz), fit=fit))
    return out, dict(median=med, scatter=mad, best=sc, z_best=float(zsc))


# ------------------------------------------------------------------ speeds and Mach numbers
def mach_bandsplit(ratio):
    """Alfven Mach number of a perpendicular shock at low beta from the band-split ratio f_U / f_L (Vrsnak et al.
    2002): X = ratio^2, M_A = sqrt(X (X + 5) / (2 (4 - X))); NaN for X >= 4."""
    X = np.asarray(ratio, float) ** 2
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(X < 4, np.sqrt(X * (X + 5) / (2 * (4 - X))), np.nan)


def model_profiles(home, p):
    """{model name: ne(r) [cm^-3]} for the Newkirk models, Saito, and the scaled MAS density along the CME axis."""
    from . import heights
    from .scattering import newkirk_profile
    prof = {f"newkirk{k}": newkirk_profile(k) for k in (1, 2, 4, 6)}
    prof["saito"] = heights.saito_ne
    try:
        from . import mast_model
        from .geometry import unit_vec
        from .paths import model_dir
        from .scattering import radial_profile
        m = mast_model.load(model_dir(home, p), density_scale=p["density_scale"], va_density=p["va_density"])
        g = p["gcs"]
        prof["mas"] = radial_profile(heights.ne_sampler(m, p), unit_vec(g["lon_stonyhurst"] + p["earth"][0], g["lat_deg"]),
                                     float(m["r"][-1]))
    except Exception:
        pass
    return prof


def height_of(ne_of_r, n, r_lo=1.0, r_hi=30.0):
    """Outermost distance [R_sun] at which ne_of_r reaches n (scanning inward from r_hi), or NaN."""
    rr = np.geomspace(r_lo + 1e-3, r_hi, 3000)
    ne = ne_of_r(rr)
    above = np.flatnonzero(ne >= n)
    if not len(above):
        return np.nan
    i = above[-1]
    if i + 1 >= len(rr):
        return float(rr[i])
    a, b = np.log(ne[i]), np.log(ne[i + 1])
    return float(rr[i] + (rr[i + 1] - rr[i]) * (a - np.log(n)) / (a - b)) if a != b else float(rr[i])


def shock_speed(f_mhz, dlnf_dt, ne_of_r, harmonic=True):
    """(speed [km/s], height [R_sun]) of a lane at f with drift d ln f / dt [1/s] in the density ne_of_r:
    V = 2 (d ln f/dt) / (d ln n/dr)."""
    n = (f_mhz * 1e6 / (2 if harmonic else 1) / 8980.0) ** 2
    r = height_of(ne_of_r, n)
    if not np.isfinite(r):
        return np.nan, np.nan
    h = 1e-3 * r
    dln = (np.log(ne_of_r(np.array([r + h]))[0]) - np.log(ne_of_r(np.array([r - h]))[0])) / (2 * h)
    return float(2 * dlnf_dt / dln * 695700.0), r


# ------------------------------------------------------------------ the analysis of a working directory
def window(home, p, s, log=print):
    """(t0, t1, source) of the type II burst: spectra.window, else the type II of the SWPC list nearest to the event
    time, else the interval of the radio images."""
    from .whitelight import utc
    if s["window"]:
        return utc(s["window"][0]), utc(s["window"][1]), "spectra.window"
    t = utc(p["time"])
    try:
        from . import bursts
        rows = [r for r in bursts.parse(bursts.fetch(t.date(), home), t.date()) if r["type"] == "II"]
        if rows:
            r = min(rows, key=lambda r: abs((r["begin"] - t).total_seconds()))
            if abs((r["begin"] - t).total_seconds()) < 3600:
                end = r["end"] or r["begin"] + timedelta(minutes=10)
                return r["begin"], end, f"SWPC type II {r['begin']:%H:%M}-{end:%H:%M} UT ({r['obs']}" + (
                    f", {r['speed_km_s']:g} km/s" if r["speed_km_s"] else "") + ")"
    except Exception as e:
        log(f"spectra: SWPC list not used ({e})")
    ri = p["radio_images"]
    if ri:
        return utc(min(e["start"] for e in ri)), utc(max(e["end"] for e in ri)), "interval of the radio images"
    return t - timedelta(minutes=5), t + timedelta(minutes=10), "event time -5 to +10 min"


def analyse_one(name, spec, t0, t1, s, prof, emission_default):
    """Lanes, speeds, and Mach numbers of one spectrum."""
    tg, fg, d, v = regrid(spec, t0, t1, f_range=s["f_range"])
    z = clean(normalise(d, spec["units"]), tg, s)
    given = [g for g in (s["lanes"] or []) if g.get("spectrum", "") in name]
    found, stats = given_lanes(z, tg, fg, t0, given) if given else lanes(z, tg, fg, s)
    res = dict(name=name, f_min=float(fg[0]), f_max=float(fg[-1]), n_t=len(tg), lanes=[], stats=stats)
    kinds = [l["kind"] for l in found]
    main_emission = s["emission"] if s["emission"] in ("F", "H") else ("F" if "harmonic" in kinds else
                                                                       "H" if "fundamental" in kinds else emission_default)
    for l in found:
        fit = l["fit"]
        if fit is None:
            continue
        tm = float(np.mean(fit["t"])) if len(fit["t"]) else fit["t_ref"]
        f_mid = float(np.exp(fit["lnf_ref"] + fit["drift"] * (tm - fit["t_ref"])))
        em = l.get("emission") or (main_emission if l["kind"] in ("main", "split_upper", "split_lower") else (
            "H" if l["kind"] == "harmonic" else "F"))
        sp = {k: shock_speed(f_mid, fit["drift"], fn, em == "H") for k, fn in prof.items()}
        dcp = None
        if v is not None and len(fit["t"]):
            ii = np.searchsorted(tg, fit["t"]).clip(0, len(tg) - 1)
            jj = np.abs(np.log(np.asarray(fit["f"], float))[:, None] - np.log(fg)[None, :]).argmin(axis=1)   # nearest channel in log f
            with np.errstate(invalid="ignore", divide="ignore"):
                dcp = float(np.nanmedian(v[jj, ii] / d[jj, ii]))
        res["lanes"].append(dict(kind=l["kind"], ratio=l["ratio"], z=l["z"], emission=em, given=bool(l.get("given")), f_mid_mhz=f_mid,
                                 t_mid=(t0 + timedelta(seconds=float(np.mean(fit["t"])) if len(fit["t"]) else fit["t_ref"])).isoformat(timespec="seconds"),
                                 drift_per_s=fit["drift"], drift_err=fit["drift_err"],
                                 df_dt_mhz_s=f_mid * fit["drift"], n_points=int(len(fit["t"])),
                                 speed_km_s={k: v_[0] for k, v_ in sp.items()}, height_rsun={k: v_[1] for k, v_ in sp.items()},
                                 dcp=dcp, t=[float(x) for x in fit["t"]], f=[float(x) for x in fit["f"]]))
    split = [l for l in res["lanes"] if l["kind"].startswith("split")]
    if split:
        l = max(split, key=lambda q: q["z"])
        ratio = l["ratio"] if l["kind"] == "split_upper" else 1.0 / l["ratio"]
        res["band_split"] = dict(ratio=ratio, X=ratio**2, mach_a=float(mach_bandsplit(ratio)))
    res["_plot"] = (tg, fg, z)
    return res


def figure(res, t0, out):
    from matplotlib.figure import Figure
    tg, fg, z = res["_plot"]
    fig = Figure(figsize=(9, 4.5))
    ax = fig.add_subplot(1, 1, 1)
    ax.pcolormesh(tg / 60.0, fg, np.nan_to_num(z, nan=0.0), cmap="viridis", vmin=-1, vmax=6, shading="auto")
    ax.set_yscale("log")
    ax.invert_yaxis()
    for l in res["lanes"]:
        if l["t"]:
            ax.plot(np.array(l["t"]) / 60.0, l["f"], ".", ms=3, color="w")
            tt = np.array([min(l["t"]), max(l["t"])])
            tm = (datetime.fromisoformat(l["t_mid"]) - t0).total_seconds()
            ax.plot(tt / 60.0, l["f_mid_mhz"] * np.exp(l["drift_per_s"] * (tt - tm)), "-", color="r", lw=1.5)
            ax.text(tt[0] / 60.0, l["f_mid_mhz"] * np.exp(l["drift_per_s"] * (tt[0] - tm)), f" {l['kind']} ({l['emission']})",
                    color="r", fontsize=8)
    ax.set_xlabel(f"minutes after {t0:%Y-%m-%d %H:%M:%S} UT")
    ax.set_ylabel("frequency (MHz)")
    ax.set_title(res["name"])
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    return out


def analyse(home, p=None, log=print):
    """Spectra, lanes, speeds, and band-split Mach numbers of a working directory; returns the summary."""
    from .paths import config, merge
    p = merge(config(home), p)
    from .pb import resolve
    p = resolve(p)
    s = settings(p)
    t0, t1, src = window(home, p, s, log)
    pad = timedelta(minutes=s["pad_min"])
    t0, t1 = t0 - pad, t1 + pad
    log(f"spectra: burst window {t0:%H:%M:%S} to {t1:%H:%M:%S} UT ({src})")
    spectra = []
    if "callisto" in s["sources"]:
        try:
            for name, paths in fetch_callisto(home, t0, t1, s["stations"], s["max_files"], log):
                try:
                    sp = read_callisto(paths)
                    sp["name"] = f"e-CALLISTO {name}"
                    spectra.append(sp)
                except Exception as e:
                    log(f"spectra: {name} unreadable ({e})")
        except OSError as e:
            log(f"spectra: e-CALLISTO archive not reachable ({e})")
    if "ovro_lwa" in s["sources"] and t0 >= datetime(2023, 1, 1):
        try:
            spectra.append(ovro_spectrum(home, t0, t1, log=log))
        except Exception as e:
            log(f"spectra: OVRO-LWA spectrum not used ({e})")
    prof = model_profiles(home, p)
    default_em = str(p.get("radio_emission", "H")).upper()
    default_em = default_em if default_em in ("F", "H") else "H"
    results = []
    for sp in spectra:
        try:
            r = analyse_one(sp["name"], sp, t0, t1, s, prof, default_em)
        except Exception as e:
            log(f"spectra: {sp['name']} failed ({e})")
            continue
        results.append(r)
        n = len(r["lanes"])
        log(f"spectra: {sp['name']}: {n} lane(s)" + (f", main lane drift {r['lanes'][0]['df_dt_mhz_s']:.3f} MHz/s at "
                                                      f"{r['lanes'][0]['f_mid_mhz']:.1f} MHz" if n else ""))
    pairs = pair_across(results, prof, s)
    sp = [l["speed_km_s"].get(s["density_model"], np.nan) for r in results for l in r["lanes"] if l["kind"] in ("main", "harmonic", "fundamental")]
    sp = np.array([v for v in sp if v is not None and np.isfinite(v)])
    mach = np.array([r["band_split"]["mach_a"] for r in results if r.get("band_split") and np.isfinite(r["band_split"]["mach_a"])])
    combined = dict(speed_median_km_s=float(np.median(sp)) if len(sp) else None,
                    speed_p16=float(np.percentile(sp, 16)) if len(sp) else None, speed_p84=float(np.percentile(sp, 84)) if len(sp) else None,
                    n_lanes=int(len(sp)), mach_a_bandsplit_median=float(np.median(mach)) if len(mach) else None, n_split=int(len(mach)))
    summ = dict(window=[t0.isoformat(), t1.isoformat()], window_source=src, pairs=pairs, combined=combined, density_scale=p["density_scale"],
                model=s["density_model"], spectra=results, reference=reference_values(home, p))
    return summ


def pair_across(results, prof, s, tol=0.3):
    """Fundamental-harmonic pairs among the main lanes of different spectra: two lanes at the same time (middle times
    within 60 s) with drift rates d ln f/dt within `tol` of each other and frequencies in the ratio 1.7 to 2.3. Unless
    spectra.emission fixes it, the lower lane is then F and the upper H, and their speeds are recomputed."""
    from datetime import datetime as dt_
    mains = [(r, l) for r in results for l in r["lanes"] if l["kind"] == "main" and not l.get("given")]
    pairs = []
    for i, (ra, a) in enumerate(mains):
        for rb, b in mains[i + 1:]:
            if ra is rb:
                continue
            lo, hi = (a, b) if a["f_mid_mhz"] < b["f_mid_mhz"] else (b, a)
            ratio = hi["f_mid_mhz"] / lo["f_mid_mhz"]
            same_t = abs((dt_.fromisoformat(a["t_mid"]) - dt_.fromisoformat(b["t_mid"])).total_seconds()) <= 60
            same_d = abs(a["drift_per_s"] - b["drift_per_s"]) <= tol * abs(max(a["drift_per_s"], b["drift_per_s"], key=abs))
            if 1.7 <= ratio <= 2.3 and same_t and same_d:
                pairs.append(dict(fundamental=f"{(ra if lo is a else rb)['name']} {lo['f_mid_mhz']:.1f} MHz",
                                  harmonic=f"{(rb if lo is a else ra)['name']} {hi['f_mid_mhz']:.1f} MHz", ratio=ratio))
                if s["emission"] == "auto":
                    for l, em in ((lo, "F"), (hi, "H")):
                        l["emission"], l["emission_reason"] = em, f"F-H pair across spectra (ratio {ratio:.2f})"
                        sp = {k: shock_speed(l["f_mid_mhz"], l["drift_per_s"], fn, em == "H") for k, fn in prof.items()}
                        l["speed_km_s"] = {k: v[0] for k, v in sp.items()}
                        l["height_rsun"] = {k: v[1] for k, v in sp.items()}
    return pairs


def reference_values(home, p):
    """Values of the reconstruction to compare with: the apex speed of the GCS fit, the speed of the SWPC list, and
    the Alfven Mach numbers of corona3d shock (analysis/shock_summary.json) when present."""
    out = {}
    g = p.get("gcs_fit") or {}
    if g.get("speed_km_s"):
        out["gcs_apex_speed_km_s"] = float(g["speed_km_s"])
    else:
        f = Path(home) / "gcs_fit.json"
        try:
            v = (json.loads(f.read_text(encoding="utf-8")).get("fit") or {}).get("speed_km_s") if f.exists() else None
        except ValueError:
            v = None
        if v:
            out["gcs_apex_speed_km_s"] = float(v)
    f = Path(home) / "analysis" / "shock_summary.json"
    if f.exists():
        try:
            sm = json.loads(f.read_text(encoding="utf-8"))
            out["shock_apex_mach_a"] = (sm.get("apex") or {}).get("mach_a")
            out["shock_radio_mach_a_median"] = ((sm.get("radio") or {}).get("mach_a") or {}).get("median")
            out["shock_flank_mach_a_median"] = ((sm.get("flank") or {}).get("mach_a") or {}).get("median")
        except ValueError:
            pass
    return out


def write(out, summ, t0):
    import csv
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for r in summ["spectra"]:
        safe = re.sub(r"[^A-Za-z0-9_-]+", "_", r["name"])
        try:
            figure(r, t0, out / f"spectrum_{safe}.png")
        except Exception:
            pass
        for l in r["lanes"]:
            row = dict(spectrum=r["name"], kind=l["kind"], ratio=l["ratio"], z=l["z"], emission=l["emission"],
                       t_mid=l["t_mid"], f_mid_mhz=l["f_mid_mhz"], df_dt_mhz_s=l["df_dt_mhz_s"],
                       drift_per_s=l["drift_per_s"], drift_err=l["drift_err"], n_points=l["n_points"], dcp=l["dcp"])
            row.update({f"speed_{k}": v for k, v in l["speed_km_s"].items()})
            row.update({f"height_{k}": v for k, v in l["height_rsun"].items()})
            rows.append(row)
    if rows:
        keys = list(dict.fromkeys(k for r in rows for k in r))
        with open(out / "spectra_lanes.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, keys)
            w.writeheader()
            w.writerows(rows)
    clean = dict(summ, spectra=[{k: v for k, v in r.items() if k != "_plot"} for r in summ["spectra"]])
    (out / "spectra_summary.json").write_text(json.dumps(clean, indent=1, default=str) + "\n", encoding="utf-8")
    return out


def _f(x, fmt=".0f"):
    return "n/a" if x is None or not np.isfinite(x) else format(x, fmt)


def report_lines(summ):
    m = summ["model"]
    L = [f"- Burst window {summ['window'][0][11:19]} to {summ['window'][1][11:19]} UT ({summ['window_source']}); spectra: "
         + (", ".join(r["name"] for r in summ["spectra"]) or "none could be read") + "."]
    ref = summ["reference"]
    for r in summ["spectra"]:
        if not r["lanes"]:
            L.append(f"- {r['name']} ({r['f_min']:.0f} to {r['f_max']:.0f} MHz): no lane above the acceptance level.")
            continue
        for l in r["lanes"]:
            if l["kind"] not in ("main", "harmonic", "fundamental"):
                continue
            sp = l["speed_km_s"]
            others = ", ".join(f"{k} {_f(v)}" for k, v in sp.items() if k != m)
            pol = f", V/I {l['dcp']:+.2f}" if l.get("dcp") is not None and np.isfinite(l["dcp"]) else ""
            L.append(f"- {r['name']}, {l['kind']} lane ({l['emission']}) at {l['f_mid_mhz']:.1f} MHz around {l['t_mid'][11:19]} UT: "
                     f"drift {l['df_dt_mhz_s']:.3f} MHz/s (d ln f/dt {l['drift_per_s']:.2e} +/- {l['drift_err']:.1e} per s){pol}, "
                     f"shock speed {_f(sp.get(m))} km/s in the {m} density (other models: {others}).")
        if r.get("band_split"):
            b = r["band_split"]
            L.append(f"- {r['name']}: band splitting f_U/f_L = {b['ratio']:.3f}, compression X = {b['X']:.2f}, Alfven Mach number "
                     f"{_f(b['mach_a'], '.2f')} (perpendicular shock, low beta).")
    c = summ.get("combined") or {}
    if c.get("n_lanes"):
        L.append(f"- Over the {c['n_lanes']} lanes the shock speed in the {m} density is {_f(c['speed_median_km_s'])} km/s (median; 16 to 84 "
                 f"percent {_f(c['speed_p16'])} to {_f(c['speed_p84'])} km/s)" + (
                     f", and the band splitting of {c['n_split']} spectra gives a median M_A of {_f(c['mach_a_bandsplit_median'], '.2f')}."
                     if c.get("n_split") else "; no band splitting was found."))
    for q in summ.get("pairs", []):
        L.append(f"- Fundamental and harmonic: {q['fundamental']} and {q['harmonic']} drift together at the frequency ratio "
                 f"{q['ratio']:.2f}, so that the lower lane is taken as F and the upper as H.")
    if ref:
        names = dict(gcs_apex_speed_km_s="apex speed of the GCS fit {} km/s", shock_apex_mach_a="M_A at the apex of the shell {}",
                     shock_radio_mach_a_median="median M_A at the radio sources {}", shock_flank_mach_a_median="median M_A on the flanks {}")
        L.append("- For comparison (corona3d gcs and shock): " + ", ".join(
            names.get(k, k + " {}").format(_f(v, '.2f') if 'mach' in k else _f(v)) for k, v in ref.items() if v is not None) + ".")
    return L


def run(home, cfg, a):
    """The `corona3d spectra` command."""
    from .report_sections import update
    summ = analyse(home)
    out = write(Path(home) / a.out, summ, datetime.fromisoformat(summ["window"][0]))
    lines = report_lines(summ)
    print("\n".join(lines))
    update(home, "Dynamic spectra", lines)
    print(f"wrote {out}")
