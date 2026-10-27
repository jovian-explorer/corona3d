"""corona3d command line: init, fetch, render."""
import argparse
import json
import shutil
import sys
from importlib.resources import files

from . import __version__, fetch, paths
from .layers import LAYERS, ON


def main(argv=None):
    ap = argparse.ArgumentParser(prog="corona3d", description="Three-dimensional reconstruction of a CME "
                                 "and the corona from public MHD, coronagraph, and radio data")
    ap.add_argument("--version", action="version", version=f"corona3d {__version__}")
    ap.add_argument("--home", help="working directory (default: CORONA3D_HOME or the current folder)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("init", help="create a working directory with the event of 2013-10-25 in config.json")
    f = sub.add_parser("fetch", help="download the MAS/MAST run and the white-light images")
    f.add_argument("--model", action="store_true", help="only the MAS/MAST run")
    f.add_argument("--wl", action="store_true", help="only the LASCO C2 and STEREO-B COR1, COR2 images")
    f.add_argument("--cor1a", action="store_true", help="also the STEREO-A COR1 images (for scripts/gcs_fit.py)")
    f.add_argument("--ovro-lwa", action="store_true",
                   help="only the OVRO-LWA images of the radio_images entries (public AWS bucket)")
    r = sub.add_parser("render", help="draw the three-panel figure",
                       epilog="layers: " + ", ".join(LAYERS))
    for k in LAYERS:
        r.add_argument(f"--{k}", action=argparse.BooleanOptionalAction, default=ON[k])
    r.add_argument("--out", default="figures/corona3d.png", help="figure path, relative to the working directory")
    r.add_argument("--half", type=float, help="panel half-height [R_sun] for all three panels")
    r.add_argument("--wl-time", help="white-light frame, e.g. 2013-10-25T15:24")
    r.add_argument("--view", action="append", default=[], metavar="PANEL=VIEW",
                   help="viewpoint of a panel: a=earth, b=north, c=stereo_b, or a Carrington "
                        "longitude and latitude of the viewing direction with an optional roll, "
                        "e.g. b=300,20 or c=120,-10,45 (repeatable)")
    r.add_argument("--radio-sources", help='JSON list, e.g. [{"freq": 173, "x": -1.1, "y": -0.1}]')
    a = ap.parse_args(argv)

    if a.cmd == "init":
        h = paths.home(a.home, must_exist=False)
        h.mkdir(parents=True, exist_ok=True)
        if (h / "config.json").exists():
            print(f"{h / 'config.json'} exists; left unchanged")
        else:
            shutil.copy(files("corona3d").joinpath("default_config.json"), h / "config.json")
            print(f"created {h / 'config.json'}; next: corona3d fetch")
        return
    h = paths.home(a.home)
    cfg = paths.config(h)
    if a.cmd == "fetch":
        failed = 0
        everything = not (a.model or a.wl or a.ovro_lwa)
        if a.model or everything:
            failed += fetch.model(h, cfg)
        if a.wl or everything:
            failed += fetch.whitelight(h, cfg, cor1a=a.cor1a)
        if a.ovro_lwa or everything:
            failed += fetch.ovro_lwa(h, cfg)
        if failed:
            sys.exit(1)
        return
    params = {}
    if a.wl_time:
        params["wl_time"] = a.wl_time
    if a.radio_sources:
        params["radio_sources"] = json.loads(a.radio_sources)
    from .reconstruction import make_figure, parse_view
    if a.view:
        views = [dict(v) for v in cfg["views"]]
        for text in a.view:
            i, v = parse_view(text)
            views[i] = dict(v, frame=views[i].get("frame", [1.6, 0.0, 0.0]) if "from" in v and
                            v["from"] == views[i].get("from") else [1.6, 0.0, 0.0])
        params["views"] = views
    make_figure(a.out, h, a.half, params, **{k: getattr(a, k) for k in LAYERS})


if __name__ == "__main__":
    main()
