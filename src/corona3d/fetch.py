"""Downloads of the public data: the PSI MAS/MAST run and the LASCO C2 and STEREO-B COR1 and
COR2 images through the Virtual Solar Observatory."""
import urllib.request

from .paths import model_dir, wl_dir

PSI = "https://www.predsci.com/data/runs/{run}/corona/{var}.hdf"
VARS = ("rho002", "br002", "bt002", "bp002", "va002")


def model(h, cfg):
    d = model_dir(h, cfg)
    d.mkdir(parents=True, exist_ok=True)
    for v in VARS:
        f = d / f"{v}.hdf"
        if f.exists():
            print(f"{f.name}: present")
            continue
        url = PSI.format(run=cfg["run"], var=v)
        print(f"{v}: {url}")
        tmp = f.with_suffix(".part")
        urllib.request.urlretrieve(url, tmp)
        tmp.replace(f)
    return d


def whitelight(h, cfg):
    from sunpy.net import Fido, attrs as a

    d = wl_dir(h)
    d.mkdir(parents=True, exist_ok=True)
    t = a.Time(*cfg["wl_window"])
    queries = [(a.Instrument("LASCO"), a.Detector("C2")),
               (a.Instrument("SECCHI"), a.Detector("COR1"), a.Source("STEREO_B")),
               (a.Instrument("SECCHI"), a.Detector("COR2"), a.Source("STEREO_B"))]
    for q in queries:
        r = Fido.search(t, *q)
        print(r)
        res = Fido.fetch(r, path=str(d / "{file}"), max_conn=4, progress=False)
        if res.errors:
            print(f"{len(res.errors)} files failed; run `corona3d fetch --wl` again to retry them")
    return d
