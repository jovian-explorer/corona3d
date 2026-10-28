"""corona3d command line: doctor, init, new, fetch, radio-interval, gcs, event, render, export."""
import argparse
import json
import shutil
import sys
from importlib.resources import files

from . import __version__, fetch, paths
from .layers import LAYERS, ON


def gcs_views(text):
    """Instrument list of --views with the legacy names stereo_a and stereo_b expanded to the COR1 and
    COR2 of that spacecraft; returns the comma-separated text."""
    out = []
    for k in text.split(","):
        k = k.strip()
        out += [f"cor1{k[-1]}", f"cor2{k[-1]}"] if k in ("stereo_a", "stereo_b") else [k]
    return ",".join(dict.fromkeys(out))


def main(argv=None):
    from ._extras import MissingExtra
    try:
        return _main(argv)
    except MissingExtra as e:
        sys.exit(f"corona3d: {e}")


def _main(argv=None):
    ap = argparse.ArgumentParser(prog="corona3d", description="Three-dimensional reconstruction of a CME "
                                 "and the corona from public MHD, coronagraph, and radio data")
    ap.add_argument("--version", action="version", version=f"corona3d {__version__}")
    ap.add_argument("--home", help="working directory (default: CORONA3D_HOME, the current folder if it holds a "
                    "config.json, else the platform cache directory; `corona3d doctor` shows it)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("doctor", help="report the installed extras, the data directories, and the reachable archives")
    d.add_argument("--offline", action="store_true", help="do not contact the archives")
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
    f.add_argument("--nrh", action="store_true",
                   help="check the NRH images of the radio_images entries and print how to obtain them "
                        "(they cannot be downloaded automatically)")
    ri = sub.add_parser("radio-interval", help="set the burst interval of the radio entries from the images")
    ri.add_argument("--search", type=float, default=10.0, metavar="MINUTES",
                    help="search the images within MINUTES of the event time (default 10)")
    ri.add_argument("--force", action="store_true", help="also replace an interval given by hand")
    g = sub.add_parser("gcs", help="fit the CME shell to the coronagraph images and write it into config.json")
    g.add_argument("--wl", help="folder of the images (default: data/wl of the working directory)")
    g.add_argument("--t-ref", help="time of the six-parameter fit (default: gcs_wl time of config.json)")
    g.add_argument("--t-out", help="time of the reconstruction (default: time of config.json)")
    g.add_argument("--flare", nargs=2, type=float, metavar=("LON", "LAT"),
                   help="Stonyhurst longitude and latitude of the flare site [deg] that seeds the envelopes")
    g.add_argument("--views", default="c2,cor1a,cor1b", type=gcs_views,
                   help="instruments: c2, cor1a, cor1b, cor2a, cor2b (legacy stereo_a and stereo_b stand for "
                        "cor1a,cor2a and cor1b,cor2b)")
    g.add_argument("--fix-direction", action="store_true",
                   help="keep the shell axis at the flare site; fit height, aspect ratio, half-angle, and tilt")
    g.add_argument("--front", help="instrument of the height-time profile (default: the limb view)")
    g.add_argument("--running-diff", type=float, metavar="MINUTES",
                   help="envelopes from running differences (frame minus the frame MINUTES earlier) "
                        "instead of differences from the first frame")
    g.add_argument("--height-out", type=float, metavar="R_SUN",
                   help="apex height at --t-out given by hand; step 2 is skipped and the height is recorded as assumed")
    g.add_argument("--height-note", help="reason for --height-out, written into the provenance record")
    g.add_argument("--drop-poor-views", dest="drop_poor_views", action="store_true", default=True,
                    help="refit without a view whose IoU is far below the others and keep the refit when it is better "
                         "(default; rule in docs/METHODS.md)")
    g.add_argument("--keep-all-views", dest="drop_poor_views", action="store_false", help="fit all views as given")
    g.add_argument("--onset", help="flare onset (ISO time) that fixes the base window of the envelopes (default: the "
                    "reconstruction time)")
    g.add_argument("--base-minutes", type=float, default=10.0, metavar="N",
                    help="the base frame is the last frame at or before onset minus N minutes (default 10), "
                         "so that the fit does not depend on when the download starts")
    g.add_argument("--skip-step1", action="store_true", help="take step 1 from the gcs_fit.json of the working directory")
    ev = sub.add_parser("event", help="new, fetch, radio-interval, gcs, and render for an event, with a report",
                        description="Runs every step for one event in the working directory and writes event_report.md; "
                                    "steps already done are skipped.")
    ev.add_argument("--home", default=argparse.SUPPRESS, help="working directory (same as the option before the command)")
    ev.add_argument("--time", required=True, help="time of the reconstruction, e.g. 2024-05-09T17:32; with only a date "
                    "(2024-05-09) the type II and III bursts of that day are listed from the SWPC event list")
    ev.add_argument("--burst", type=int, metavar="N",
                    help="with a date as --time: use the begin time of burst number N of the list")
    ev.add_argument("--flare", required=True, help="flare position, e.g. S14W28")
    ev.add_argument("--fit", choices=("auto", "fixed", "skip"), default="auto",
                    help="auto: free GCS fit, then the fixed-direction fit when it fails or is poor; "
                         "fixed: only the fixed-direction fit; skip: keep the shell of config.json")
    ev.add_argument("--iou-min", type=float, default=0.5,
                    help="mean intersection over union below which a fit is poor (default 0.5)")
    ev.add_argument("--height-out", type=float, metavar="R_SUN", help="apex height at --time given by hand")
    ev.add_argument("--height-note", help="reason for --height-out (required with it)")
    ev.add_argument("--views", type=gcs_views, help="instruments of the fit (default: c2, cor1 and cor2 of the STEREO spacecraft)")
    ev.add_argument("--t-ref", help="time of the six-parameter fit (default: the gcs_wl time of config.json)")
    ev.add_argument("--drop-poor-views", dest="drop_poor_views", action="store_true", default=True,
                    help="refit without a view whose IoU is far below the others and keep the refit when it is better "
                         "(default; rule in docs/METHODS.md)")
    ev.add_argument("--keep-all-views", dest="drop_poor_views", action="store_false", help="fit all views as given")
    ev.add_argument("--onset", help="flare onset (ISO time) that fixes the base window of the envelopes (default: the "
                    "reconstruction time)")
    ev.add_argument("--base-minutes", type=float, default=10.0, metavar="N",
                    help="the base frame is the last frame at or before onset minus N minutes (default 10), "
                         "so that the fit does not depend on when the download starts")
    ev.add_argument("--radio-search", type=float, default=10.0, metavar="MINUTES",
                    help="radio images are searched within MINUTES of --time (default 10)")
    ev.add_argument("--force", action="store_true", help="repeat the downloads, the fit, and the figures")
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
    e = sub.add_parser("export", help="write every enabled layer to data files with a JSON metadata sidecar",
                       epilog="layers: " + ", ".join(LAYERS))
    for k in LAYERS:
        if k != "sphere":
            e.add_argument(f"--{k}", action=argparse.BooleanOptionalAction, default=ON[k])
    e.add_argument("--format", choices=("npz", "csv", "vtk"), default="npz")
    e.add_argument("--out", default="export", help="output folder, relative to the working directory")
    e.add_argument("--radio-sources", help='JSON list, as for render')
    a = ap.parse_args(argv)

    if a.cmd == "doctor":
        from . import doctor
        sys.exit(doctor.run(a.home, network=not a.offline))
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
    if a.cmd == "event":
        from . import bursts, pipeline
        if bursts.is_date(a.time):
            try:
                t = bursts.choose(a.time.strip(), a.burst, paths.home(a.home, must_exist=False))
            except RuntimeError as e:
                sys.exit(str(e))
            if t is None:
                return
            a.time = t
            print(f"--burst {a.burst}: time {t}")
        elif a.burst is not None:
            sys.exit("--burst needs a date without a time as --time, e.g. --time 2024-05-09")
        if a.height_out is not None and not a.height_note:
            sys.exit("--height-out needs --height-note, the reason for the assumed height")
        sys.exit(pipeline.run(paths.home(a.home, must_exist=False), a))
    h = paths.home(a.home)
    cfg = paths.config(h)
    if a.cmd == "radio-interval":
        from . import radio
        radio.run(h, cfg, a)
        return
    if a.cmd == "fetch":
        failed, parts = 0, {}
        everything = not (a.model or a.wl or a.ovro_lwa or a.nrh)
        if a.model or everything:
            parts["model"] = fetch.model(h, cfg)
        if a.wl or everything:
            parts["white light"] = fetch.whitelight(h, cfg, cor1a=a.cor1a)
        if a.ovro_lwa or everything:
            parts["OVRO-LWA"] = fetch.ovro_lwa(h, cfg)
        if a.nrh or everything:
            fetch.nrh(h, cfg)                  # a report only: no automatic download exists
        failed = sum(parts.values())
        if failed:
            sys.exit(f"corona3d fetch: {failed} files failed or were damaged ("
                     + ", ".join(f"{k} {n}" for k, n in parts.items() if n)
                     + "), usually a transient server or network error; complete files are kept, so give the same "
                       "command again")
        return
    if a.cmd == "gcs":
        from . import gcs
        gcs.run(h, cfg, a)
        return
    if a.cmd == "export":
        from .exporting import export
        params = {"radio_sources": json.loads(a.radio_sources)} if a.radio_sources else None
        d = export(h, a.out, a.format, params=params,
                   **{k: getattr(a, k) for k in LAYERS if k != "sphere"})
        print("exported to", d)
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
            sc = v.pop("sc", None)                     # legacy stereo_a / stereo_b
            if sc and sc != cfg["stereo_sc"]:
                sys.exit(f"--view {text}: this event uses STEREO-{cfg['stereo_sc']} (config.json stereo_sc)")
            views[i] = dict(v, frame=views[i].get("frame", [1.6, 0.0, 0.0]) if "from" in v and
                            v["from"] == views[i].get("from") else [1.6, 0.0, 0.0])
        params["views"] = views
    make_figure(a.out, h, a.half, params, **{k: getattr(a, k) for k in LAYERS})


if __name__ == "__main__":
    main()
