"""corona3d command line: init, fetch, render."""
import argparse
import json
import shutil
from importlib.resources import files

from . import __version__, fetch, paths
from .reconstruction import LAYERS, ON, make_figure


def main(argv=None):
    ap = argparse.ArgumentParser(prog="corona3d", description="Three-dimensional reconstruction of a CME and the corona "
                                 "from public MHD, coronagraph and radio data")
    ap.add_argument("--version", action="version", version=f"corona3d {__version__}")
    ap.add_argument("--home", help="working directory (default: CORONA3D_HOME or the current folder)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("init", help="create a working directory with the event of 2013-10-25 in config.json")
    f = sub.add_parser("fetch", help="download the MAS/MAST run and the white-light images")
    f.add_argument("--model", action="store_true", help="only the MAS/MAST run")
    f.add_argument("--wl", action="store_true", help="only the LASCO C2 and STEREO-B COR1, COR2 images")
    r = sub.add_parser("render", help="draw the three-panel figure")
    for k in LAYERS:
        r.add_argument(f"--{k}", action=argparse.BooleanOptionalAction, default=ON[k])
    r.add_argument("--out", default="figures/corona3d.png")
    r.add_argument("--half", type=float, help="panel half-height [R_sun] for all three panels")
    r.add_argument("--wl-time", help="white-light frame, e.g. 2013-10-25T15:24")
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
        if a.model or not a.wl:
            fetch.model(h, cfg)
        if a.wl or not a.model:
            fetch.whitelight(h, cfg)
        return
    if a.wl_time:
        cfg["wl_time"] = a.wl_time
    if a.radio_sources:
        cfg["radio_sources"] = json.loads(a.radio_sources)
    make_figure(h / a.out, h, a.half, cfg, **{k: getattr(a, k) for k in LAYERS})


if __name__ == "__main__":
    main()
