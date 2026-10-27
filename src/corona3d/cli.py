"""corona3d command line: init, new, fetch, gcs, render."""
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
    i = sub.add_parser("init", help="create a working directory with an event in config.json")
    i.add_argument("--event", default="2013-10-25", help="event preset: " + ", ".join(paths.events()))
    n = sub.add_parser("new", help="write config.json for a new event from its time and flare position")
    n.add_argument("--time", required=True, help="time of the reconstruction, e.g. 2024-05-09T17:32")
    n.add_argument("--flare", required=True, help="flare position, e.g. S14W28")
    n.add_argument("--force", action="store_true", help="replace an existing config.json")
    f = sub.add_parser("fetch", help="download the MAS/MAST run and the white-light images")
    f.add_argument("--model", action="store_true", help="only the MAS/MAST run")
    f.add_argument("--wl", action="store_true", help="only the LASCO C2 and STEREO COR1, COR2 images")
    f.add_argument("--cor1a", action="store_true", help="also the STEREO-A COR1 images (for corona3d gcs)")
    f.add_argument("--ovro-lwa", action="store_true",
                   help="only the OVRO-LWA images of the radio_images entries (public AWS bucket)")
    g = sub.add_parser("gcs", help="fit the CME shell to the coronagraph images and write it into config.json")
    g.add_argument("--wl", help="folder of the images (default: data/wl of the working directory)")
    g.add_argument("--t-ref", help="time of the six-parameter fit (default: gcs_wl time of config.json)")
    g.add_argument("--t-out", help="time of the reconstruction (default: time of config.json)")
    g.add_argument("--flare", nargs=2, type=float, metavar=("LON", "LAT"),
                   help="Stonyhurst longitude and latitude of the flare site [deg] that seeds the envelopes")
    g.add_argument("--views", default="c2,cor1a,cor1b", help="instruments: c2, cor1a, cor1b, cor2a, cor2b")
    g.add_argument("--fix-direction", action="store_true",
                   help="keep the shell axis at the flare site; fit height, aspect ratio, half-angle, and tilt")
    g.add_argument("--front", help="instrument of the height-time profile (default: the limb view)")
    g.add_argument("--running-diff", type=float, metavar="MINUTES",
                   help="envelopes from running differences (frame minus the frame MINUTES earlier) "
                        "instead of differences from the first frame")
    g.add_argument("--height-out", type=float, metavar="R_SUN",
                   help="apex height at --t-out given by hand; step 2 is skipped and the height is recorded as assumed")
    g.add_argument("--height-note", help="reason for --height-out, written into the provenance record")
    g.add_argument("--skip-step1", action="store_true", help="take step 1 from the gcs_fit.json of the working directory")
    r = sub.add_parser("render", help="draw the three-panel figure",
                       epilog="layers: " + ", ".join(LAYERS))
    for k in LAYERS:
        r.add_argument(f"--{k}", action=argparse.BooleanOptionalAction, default=ON[k])
    r.add_argument("--out", default="figures/corona3d.png", help="figure path, relative to the working directory")
    r.add_argument("--half", type=float, help="panel half-height [R_sun] for all three panels")
    r.add_argument("--wl-time", help="white-light frame, e.g. 2013-10-25T15:24")
    r.add_argument("--view", action="append", default=[], metavar="PANEL=VIEW",
                   help="viewpoint of a panel: a=earth, b=north, c=stereo, or a Carrington "
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
            if a.event == "2013-10-25":
                shutil.copy(files("corona3d").joinpath("default_config.json"), h / "config.json")
            elif a.event in paths.events():
                cfg = paths.merge(paths.default_config(), paths.event(a.event))
                (h / "config.json").write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")
            else:
                sys.exit(f"unknown event {a.event}; presets: " + ", ".join(paths.events()))
            print(f"created {h / 'config.json'} ({a.event}); next: corona3d fetch")
        return
    if a.cmd == "new":
        from .events import new_config
        h = paths.home(a.home, must_exist=False)
        h.mkdir(parents=True, exist_ok=True)
        if (h / "config.json").exists() and not a.force:
            sys.exit(f"{h / 'config.json'} exists; give --force to replace it")
        cfg = paths.merge(paths.default_config(), new_config(a.time, a.flare))
        (h / "config.json").write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")
        print(f"created {h / 'config.json'} for {a.time} ({cfg['run']}, STEREO-{cfg['stereo_sc']}); "
              "next: corona3d fetch, then fit the shell with corona3d gcs")
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
    if a.cmd == "gcs":
        from . import gcs
        gcs.run(h, cfg, a)
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
