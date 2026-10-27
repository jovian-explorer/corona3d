"""The configuration of a new event from its time and flare position (`corona3d new`).

The observers are placed with SunPy (Earth) and JPL Horizons (STEREO-A, or STEREO-B before
its loss of contact in October 2014), the PSI MAS/MAST run is the one of the Carrington rotation
of the event, and the time windows of the coronagraph and radio images are set around the event.
The GCS shell starts above the flare site and is meant to be fitted with `corona3d gcs`.
"""
import re
import urllib.request
from datetime import datetime, timedelta

PSI_RUNS = "https://www.predsci.com/data/runs/"
OVRO_LWA_START = datetime(2024, 4, 1)
STEREO_B_END = datetime(2014, 10, 1)


def parse_flare(text):
    """'S14W28' -> (Stonyhurst longitude, latitude) = (28.0, -14.0) [deg]."""
    m = re.fullmatch(r"([NS])(\d+(?:\.\d+)?)([EW])(\d+(?:\.\d+)?)", text.strip().upper())
    if not m:
        raise ValueError(f"flare position {text!r}: give it as e.g. S14W28")
    lat = float(m.group(2)) * (1 if m.group(1) == "N" else -1)
    lon = float(m.group(4)) * (1 if m.group(3) == "W" else -1)
    return lon, lat


def observers(t, sc):
    """Carrington longitude and latitude [deg] of Earth and of STEREO-`sc` at time t."""
    from sunpy.coordinates import HeliographicCarrington, get_earth, get_horizons_coord
    out = []
    for c in (get_earth(t), get_horizons_coord(f"STEREO-{sc}", t)):
        h = c.transform_to(HeliographicCarrington(observer="self", obstime=t))
        out.append([round(float(h.lon.deg), 3), round(float(h.lat.deg), 3)])
    return out


def psi_run(t, listing=None):
    """Path of the PSI MAS/MAST run of the Carrington rotation of t, preferring the medium
    resolution; `listing` is the text of the PSI run index (downloaded when not given)."""
    from sunpy.coordinates import sun
    cr = int(sun.carrington_rotation_number(t))
    if listing is None:
        with urllib.request.urlopen(PSI_RUNS, timeout=60) as r:
            listing = r.read().decode()
    for res in ("medium", "high"):
        if f"cr{cr}-{res}" in listing:
            return f"cr{cr}-{res}/hmi_mast_mas_std_0201"
    raise ValueError(f"no PSI run listed for Carrington rotation {cr} at {PSI_RUNS}")


def new_config(time, flare, earth_stereo=None, run=None):
    """Configuration keys of a new event; `earth_stereo` and `run` can be given to work offline."""
    t = datetime.fromisoformat(time)
    lon, lat = parse_flare(flare)
    sc = "B" if t < STEREO_B_END else "A"
    earth, stereo = earth_stereo or observers(time, sc)
    iso = lambda dt: dt.strftime("%Y-%m-%dT%H:%M:%S")
    shell = dict(kappa=0.4, alpha_deg=30.0, lon_stonyhurst=lon, lat_deg=lat, tilt_deg=-90.0)
    radio = [dict(instrument="nrh", folder="data/nrh", freqs=[150.9, 228.0], start=iso(t),
                  end=iso(t + timedelta(minutes=5)), step_s=30, level=0.5, min_tb=1e8, contours=False)]
    if t >= OVRO_LWA_START:
        radio.append(dict(instrument="ovro_lwa", folder="data/ovro_lwa", freqs=[34.1, 52.4, 70.8, 84.6],
                          start=iso(t), end=iso(t + timedelta(minutes=5)), step_s=30, level=0.5,
                          min_tb=1e8, contours=False))
    return {
        "event": f"event of {t:%Y-%m-%d %H:%M} UT, flare {flare.upper()} (set up with corona3d new)",
        "time": iso(t),
        "run": run or psi_run(time),
        "model_dlon": 0.0,
        "earth": earth,
        "stereo": stereo,
        "stereo_sc": sc,
        "gcs": dict(shell, height=1.5),
        "gcs_wl": dict(shell, height=5.0, time=iso(t + timedelta(minutes=60))),
        "radio_images": radio,
        "wl_window": [iso(t - timedelta(minutes=60)), iso(t + timedelta(minutes=150))],
        "wl_time": iso(t + timedelta(minutes=60)),
        "views": [{"from": "earth", "frame": [1.6, 0.0, 0.0]}, {"from": "north", "frame": [1.6, 0.0, 0.0]},
                  {"from": "stereo", "frame": [1.6, 0.0, 0.0]}],
    }
