"""`corona3d doctor`: what is installed, where the data go, and which archives answer."""
import platform
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from importlib import metadata
from importlib.resources import files

from . import __version__, _extras, paths

CORE = ["numpy", "scipy", "matplotlib", "astropy", "sunpy", "pyhdf", "platformdirs"]

ARCHIVES = {                                   # name: (url that answers a GET, what it serves)
    "PSI MAS/MAST runs": ("https://www.predsci.com/data/runs/", "coronal model (fetch --model)"),
    "VSO (LASCO, SECCHI)": ("https://sdac.virtualsolar.org/cgi/search", "coronagraph images (fetch --wl)"),
    "OVRO-LWA (AWS)": ("https://ovro-lwa-solar.s3.us-west-2.amazonaws.com/?max-keys=1", "radio images (fetch --ovro-lwa)"),
    "NOAA NCEI (SWPC events)": ("https://www.ngdc.noaa.gov/stp/space-weather/swpc-products/daily_reports/"
                                 "solar_event_reports/", "burst list for a date (event --time DATE)"),
    "JPL Horizons": ("https://ssd.jpl.nasa.gov/api/horizons.api?format=json&COMMAND=%27399%27", "STEREO and Earth positions"),
}


def _ver(pkg):
    try:
        return metadata.version(pkg)
    except metadata.PackageNotFoundError:
        return None


def reach(url, timeout=8):
    """'ok (HTTP 200)' or the reason the archive did not answer."""
    try:
        req = urllib.request.Request(url, headers={"User-Agent": f"corona3d/{__version__}"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return f"ok (HTTP {r.status})"
    except Exception as e:
        return f"not reached ({getattr(e, 'reason', e) or type(e).__name__})"


def run(home_arg=None, network=True, out=print):
    """Print the report; returns 1 when a core package is missing, else 0."""
    out(f"corona3d {__version__}  Python {platform.python_version()}  {platform.system()} {platform.machine()}")
    out("\nCore packages")
    bad = 0
    for p in CORE:
        v = _ver(p)
        bad += v is None
        out(f"  {p:<12} {v or 'MISSING (reinstall corona3d)'}")
    out("\nOptional extras (pip install \"corona3d[NAME]\"; [all] installs every one)")
    for extra, pk in _extras.status().items():
        ok = all(pk.values())
        out(f"  [{extra}] {'installed' if ok else 'not installed'}: "
            + ", ".join(f"{n} {v or '-'}" for n, v in pk.items()))
    h, how = paths.home_source(home_arg)
    out(f"\nWorking directory: {h}  (from {how}; change with --home or CORONA3D_HOME)")
    out(f"  config.json: {'present' if (h / 'config.json').exists() else 'absent (corona3d init, new or event creates it)'}")
    if (h / "config.json").exists():
        cfg = paths.config(h)
        md = paths.model_dir(h, cfg)
        have = [v for v in ("rho002", "br002", "bt002", "bp002", "va002") if (md / f"{v}.hdf").exists()]
        out(f"  model run {cfg['run']}: {len(have)} of 5 files in {md}")
        wl = paths.wl_dir(h)
        out(f"  white-light images: {len(list(wl.glob('*.fts'))) if wl.exists() else 0} files in {wl}")
    out(f"  default home (no choice made): {paths.default_home()}")
    out(f"\nPackage data: {files('corona3d')}  presets: {', '.join(paths.events())}")
    if network:
        out("\nArchives (HTTP check, 8 s limit each)")
        with ThreadPoolExecutor(len(ARCHIVES)) as ex:
            res = dict(zip(ARCHIVES, ex.map(lambda k: reach(ARCHIVES[k][0]), ARCHIVES)))
        for k, r in res.items():
            out(f"  {k:<26} {r}  - {ARCHIVES[k][1]}")
        out("  NRH 2-D images: no automatic download; request them at https://rsdb.obs-nancay.fr/")
    else:
        out("\nArchives: not checked (--offline)")
    return 1 if bad else 0
