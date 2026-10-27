"""Downloads of the public data: the PSI MAS/MAST run, and the LASCO C2 and STEREO-B COR1 and
COR2 images (with --cor1a also the STEREO-A COR1 images used by scripts/gcs_fit.py) through
the Virtual Solar Observatory."""
import urllib.request

from .paths import model_dir, wl_dir

PSI = "https://www.predsci.com/data/runs/{run}/corona/{var}.hdf"
VARS = ("rho002", "br002", "bt002", "bp002", "va002")


def model(h, cfg):
    """Download the five files of the run; returns the number of files that failed."""
    d = model_dir(h, cfg)
    d.mkdir(parents=True, exist_ok=True)
    failed = 0
    for v in VARS:
        f = d / f"{v}.hdf"
        if f.exists():
            print(f"{f.name}: present")
            continue
        url = PSI.format(run=cfg["run"], var=v)
        print(f"{v}: {url}")
        tmp = f.with_suffix(".part")
        try:
            with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as o:
                while chunk := r.read(1 << 20):
                    o.write(chunk)
            tmp.replace(f)                     # the final name only for a complete file
        except OSError as e:
            tmp.unlink(missing_ok=True)
            print(f"{v}: failed ({e})")
            failed += 1
    return failed


def _readable(f):
    from astropy.io import fits
    try:
        with fits.open(f) as hd:
            hd[0].data                         # reads the whole array: a truncated file fails
        return True
    except Exception:
        return False


def whitelight(h, cfg, cor1a=False):
    """Download the images of wl_window; damaged files are deleted so that the next run fetches
    them again. Returns the number of files that failed or were damaged."""
    from sunpy.net import Fido, attrs as a

    d = wl_dir(h)
    d.mkdir(parents=True, exist_ok=True)
    t = a.Time(*cfg["wl_window"])
    queries = [(a.Instrument("LASCO"), a.Detector("C2")),
               (a.Instrument("SECCHI"), a.Detector("COR1"), a.Source("STEREO_B")),
               (a.Instrument("SECCHI"), a.Detector("COR2"), a.Source("STEREO_B"))]
    if cor1a:
        queries.append((a.Instrument("SECCHI"), a.Detector("COR1"), a.Source("STEREO_A")))
    failed = 0
    for q in queries:
        r = Fido.search(t, *q)
        print(r)
        res = Fido.fetch(r, path=str(d / "{file}"), max_conn=4, progress=False)
        failed += len(res.errors)
    bad = [f for f in d.glob("*.fts") if not _readable(f)]
    for f in bad:
        print(f"{f.name}: damaged, deleted")
        f.unlink()
    failed += len(bad)
    if failed:
        print(f"{failed} files failed; run `corona3d fetch --wl` again to fetch them")
    return failed


def ovro_lwa(h, cfg):
    """Download the OVRO-LWA band-averaged (mfs) image files of every radio_images entry with
    instrument "ovro_lwa", between its start and end, from the public bucket into its folder;
    returns the number of files that failed."""
    from datetime import timedelta

    from . import ovro_lwa as lwa
    from .whitelight import utc
    failed = 0
    for rs in cfg["radio_images"]:
        if rs["instrument"] != "ovro_lwa":
            continue
        t0, t1 = utc(rs["start"]) - timedelta(minutes=1), utc(rs["end"]) + timedelta(minutes=1)
        d = h / rs["folder"]
        d.mkdir(parents=True, exist_ok=True)
        keys = [k for day in {t0.date(), t1.date()} for k in lwa.list_remote(day)]
        keys = [k for k in keys if t0 <= lwa.file_time(k) <= t1]
        print(f"OVRO-LWA: {len(keys)} files between {t0:%H:%M} and {t1:%H:%M} UT")
        for k in keys:
            f = d / k.split("/")[-1]
            if f.exists():
                continue
            tmp = f.with_suffix(".part")
            try:
                with urllib.request.urlopen(f"{lwa.S3}/{k}", timeout=120) as r, open(tmp, "wb") as o:
                    while chunk := r.read(1 << 20):
                        o.write(chunk)
                tmp.replace(f)
            except OSError as e:
                tmp.unlink(missing_ok=True)
                print(f"{f.name}: failed ({e})")
                failed += 1
    return failed
