"""Solar Orbiter as an observer: Metis, EUI/FSI and SoloHI images from the Solar Orbiter Archive (SOAR).

`corona3d fetch --solo` downloads, through the SOAR client of SunPy and one file at a time, the level-2
total-brightness (tB) and polarised-brightness (pB) images of the visible-light channel of Metis
(Antonucci et al. 2020) into data/wl, the 174 A images of the Full Sun Imager (FSI) of EUI (Rochus et
al. 2020) into data/euv, and with --solohi the SoloHI images (Howard et al. 2020) into data/solohi.
The position of the spacecraft is then read from the FITS header of a downloaded image (HGLN_OBS,
HGLT_OBS, DSUN_OBS) and written to config.json as the Carrington longitude and latitude `solo` and the
distance `solo_dsun_au`, so that a panel can be drawn from Solar Orbiter ("from": "solo").

Times: Solar Orbiter is closer to the Sun than Earth, so that it records an event earlier than an
observer at 1 AU. Every Solar Orbiter frame is therefore given the time at which its light would reach
1 AU, its time plus (1 AU - DSUN_OBS) / c, which is the time convention of LASCO, AIA and the radio data
(156 s on 2024 May 9, at 0.688 AU), and a Metis image, which is summed over a polarisation sequence of
about 30 min, is given the middle of that sequence (DATE-AVG). DATE_EAR of the headers is not used,
because the SoloHI files of 2024 May 9 give it 76 s before DATE-OBS.

Field of view of Metis: the inner and outer edge of the annulus in degrees (INN_FOV, OUT_FOV of the
header, or the nominal 1.6 and 2.9 deg) divided by the angular radius of the Sun seen from the
spacecraft (RSUN_ARC), so that it scales with the distance of Solar Orbiter from the Sun.
"""
import json
from datetime import datetime, timedelta
from pathlib import Path

AU_M = 149597870700.0
C_M_S = 299792458.0
METIS_FOV_DEG = (1.6, 2.9)              # nominal inner and outer edge of the Metis field [deg]
METIS_SEED = 0.25                       # seed of the envelope at this fraction of the annulus above its inner edge
PAD_MIN = 60.0                          # Metis is searched from this many minutes before wl_window (30 min sequences)
METIS_GLOB = "solo_L2_metis-vl-tb_*.fits"


def _iso(text):
    return datetime.fromisoformat(str(text).replace("Z", "")[:23])


def frame_time(header):
    """Time of a Solar Orbiter frame as seen from Earth: the middle of the exposure or sequence (DATE-AVG,
    else DATE-OBS) plus the light travel time (D_Earth - DSUN_OBS) / c, where D_Earth is the actual Sun-Earth
    distance at that time (sunpy get_earth; 1 AU if the ephemeris is not available), which differs from 1 AU by
    up to 1.7 percent (8 s of light time)."""
    t = _iso(header["DATE-AVG"]) if header.get("DATE-AVG") else _iso(header["DATE-OBS"])
    try:
        from sunpy.coordinates import get_earth
        d_earth = float(get_earth(t).radius.to_value("m"))
    except Exception:
        d_earth = AU_M
    t += timedelta(seconds=(d_earth - float(header["DSUN_OBS"])) / C_M_S)
    return (t + timedelta(microseconds=500000)).replace(microsecond=0)        # to the nearest second


def header_of(f):
    """Header of the image of a Solar Orbiter file (the first extension of a compressed EUI file)."""
    from astropy.io import fits
    with fits.open(f) as hd:
        for h in hd:
            if h.header.get("NAXIS", 0) >= 2 or "ZNAXIS" in h.header:
                return h.header.copy()
        return hd[0].header.copy()


def metis_fov(header):
    """(r_in, r_out, r_seed) [R_sun] of a Metis image: the edges of the field in degrees (INN_FOV,
    OUT_FOV, else the nominal 1.6 and 2.9 deg) divided by the angular solar radius (RSUN_ARC, else
    from DSUN_OBS)."""
    import numpy as np
    rs = (float(header["RSUN_ARC"]) / 3600.0 if header.get("RSUN_ARC")
          else float(np.degrees(np.arcsin(695700e3 / float(header["DSUN_OBS"])))))
    lo = float(header.get("INN_FOV") or METIS_FOV_DEG[0])
    hi = float(header.get("OUT_FOV") or METIS_FOV_DEG[1])
    r_in, r_out = lo / rs, hi / rs
    return round(r_in, 3), round(r_out, 3), round(r_in + METIS_SEED * (r_out - r_in), 3)


def load_map(f):
    """SunPy map of the image of a Solar Orbiter file. The level-2 images are calibrated (MSB for Metis,
    DN/s for EUI), so that they are not divided by the exposure time; the pixels of a Metis image that are
    exactly zero lie behind the occulter or outside the field stop and become NaN."""
    import numpy as np
    import sunpy.map
    from astropy.io import fits
    with fits.open(f) as hd:
        x = next((h for h in hd if h.data is not None and np.ndim(h.data) == 2), hd[0])
        d, hdr = x.data.astype(float), x.header.copy()
    if str(hdr.get("INSTRUME", "")).upper() == "METIS":
        d[d == 0] = np.nan
    return sunpy.map.Map(d, hdr)


def metis_frames(wl):
    """[(Earth time, [file])] of the Metis tB images in a folder, earliest first."""
    return sorted(((frame_time(header_of(f)), [f]) for f in Path(wl).glob(METIS_GLOB)), key=lambda x: x[0])


def observer(header, time=None):
    """(Carrington longitude, latitude [deg], distance [AU]) of the spacecraft from HGLN_OBS, HGLT_OBS
    and DSUN_OBS of a header, in the Carrington frame of corona3d (observer "self", as corona3d new), at
    `time` (default DATE-OBS; the Stonyhurst position changes by less than 0.1 deg in a day)."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from sunpy.coordinates import HeliographicCarrington, HeliographicStonyhurst
    t = time or header["DATE-OBS"]
    c = SkyCoord(float(header["HGLN_OBS"]) * u.deg, float(header["HGLT_OBS"]) * u.deg, float(header["DSUN_OBS"]) * u.m,
                 frame=HeliographicStonyhurst, obstime=t)
    h = c.transform_to(HeliographicCarrington(observer="self", obstime=t))
    return round(float(h.lon.deg) % 360.0, 3), round(float(h.lat.deg), 3), round(float(header["DSUN_OBS"]) / AU_M, 4)


def files(h):
    """Solar Orbiter image files of a working directory (Metis, then EUI, then SoloHI)."""
    h = Path(h)
    return (sorted((h / "data" / "wl").glob("solo_L2_metis-vl-*.fits")) + sorted((h / "data" / "euv").glob("solo_L2_eui-*.fits"))
            + sorted((h / "data" / "solohi").glob("solo_L2_solohi-*.fits")))


def set_observer(h, f=None):
    """Write `solo` (Carrington lon, lat at the time of the reconstruction) and `solo_dsun_au` into
    config.json from the header of the first Solar Orbiter file (or of f); returns the values, or None
    without a file."""
    fs = [f] if f else files(h)
    if not fs:
        return None
    hd = header_of(fs[0])
    cf = Path(h) / "config.json"
    user = json.loads(cf.read_text(encoding="utf-8"))
    lon, lat, au = observer(hd, user.get("time"))
    user.update(solo=[lon, lat], solo_dsun_au=au)
    cf.write_text(json.dumps(user, indent=2) + "\n", encoding="utf-8")
    print(f"Solar Orbiter at Carrington lon {lon:.3f}, lat {lat:.3f} deg, {au:.3f} AU (Stonyhurst lon "
          f"{float(hd['HGLN_OBS']):.1f} deg, from {Path(fs[0]).name}); written to {cf}")
    return lon, lat, au


def _readable(f):
    try:
        import numpy as np
        from astropy.io import fits
        with fits.open(f) as hd:
            for x in hd:
                if x.data is not None:
                    np.asarray(x.data).sum()       # reads the whole array: a truncated file fails
        return True
    except Exception:
        return False


def fetch(h, cfg, solohi=False):
    """Download the Metis tB and pB images (from PAD_MIN minutes before wl_window to its end), the EUI/FSI
    174 A images of wl_window, and with `solohi` the SoloHI images of wl_window, one file at a time; then
    write the observer into config.json. Returns the number of files that failed or were damaged."""
    from sunpy.net import Fido, attrs as a
    from .whitelight import utc
    h = Path(h)
    t0, t1 = utc(cfg["wl_window"][0]), utc(cfg["wl_window"][1])
    wl, euv, hi = h / "data" / "wl", h / "data" / "euv", h / "data" / "solohi"
    jobs = [(a.Time(t0 - timedelta(minutes=PAD_MIN), t1), (a.Instrument("METIS"), a.soar.Product("metis-vl-tb")), wl),
            (a.Time(t0 - timedelta(minutes=PAD_MIN), t1), (a.Instrument("METIS"), a.soar.Product("metis-vl-pb")), wl),
            (a.Time(t0, t1), (a.Instrument("EUI"), a.soar.Product("eui-fsi174-image")), euv)]
    if solohi:
        jobs.append((a.Time(t0, t1), (a.Instrument("SOLOHI"),), hi))
    failed = 0
    for tr, q, d in jobs:
        d.mkdir(parents=True, exist_ok=True)
        try:
            r = Fido.search(tr, *q, a.Level(2))
        except Exception as e:                       # an archive that does not answer: the others still run
            print(f"SOAR search failed ({e}); run corona3d fetch --solo again")
            failed += 1
            continue
        r = r["soar"] if "soar" in r.keys() else r   # SoloHI is also listed by the VSO; the SOAR copy is used
        print(r)
        res = Fido.fetch(r, path=str(d / "{file}"), max_conn=1, progress=False)
        for e in res.errors:
            print(f"download failed: {e}")
        failed += len(res.errors)
    bad = [f for f in files(h) if not _readable(f)]
    for f in bad:
        print(f"{f.name}: damaged, deleted")
        f.unlink()
    failed += len(bad)
    set_observer(h)
    return failed
