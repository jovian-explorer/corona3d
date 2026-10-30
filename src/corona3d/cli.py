"""corona3d command line: doctor, init, new, fetch, radio-interval, gcs, event, render, export, shock, pb, sources,
spectra."""
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


OPACITY_KEYS = ("sphere", "surface", "open", "closed", "streamer", "streamer_with_radio", "gcs", "radio",
                "radio_contours", "wl")


def floats(text):
    """'0.3,0.5,0.9' -> [0.3, 0.5, 0.9]."""
    try:
        return [float(x) for x in text.split(",")]
    except ValueError:
        raise argparse.ArgumentTypeError(f"{text!r} is not a comma list of numbers")


def parse_opacity(text):
    """'radio_contours=0.2,1' -> ('radio_contours', [0.2, 1.0]); 'wl=0.5' -> ('wl', 0.5);
    'surface=linear' -> ('surface', 'linear')."""
    k, _, v = text.partition("=")
    k = k.strip()
    if k not in OPACITY_KEYS or not v:
        sys.exit(f"--opacity {text}: give LAYER=VALUE with LAYER one of {', '.join(OPACITY_KEYS)}")
    try:
        vals = floats(v)
    except argparse.ArgumentTypeError:
        return k, v.strip()
    return k, vals[0] if len(vals) == 1 else vals


def figure_options(r):
    """The layers, the panels, and the other options that corona3d render and corona3d animate share."""
    for k in LAYERS:
        r.add_argument(f"--{k}", action=argparse.BooleanOptionalAction, default=ON[k])
    r.add_argument("--half", type=float, help="panel half-height [R_sun] for all panels")
    r.add_argument("--shell-on-top", action="store_true",
                   help="draw the GCS shell over the model surface, so that a shell close to the Sun stays visible")
    r.add_argument("--view", action="append", default=[], metavar="PANEL=VIEW",
                   help="viewpoint of a panel: a=earth, b=north, c=stereo, c=solo (Solar Orbiter), or a Carrington "
                        "longitude and latitude of the viewing direction with an optional roll, "
                        "e.g. b=300,20 or c=120,-10,45 (repeatable)")
    r.add_argument("--radio-sources", help='JSON list, e.g. [{"freq": 173, "x": -1.1, "y": -0.1}]')
    r.add_argument("--opacity", action="append", default=[], metavar="LAYER=VALUE",
                   help="opacity of a layer from 0 (transparent) to 1: one number, or for radio, radio_contours, "
                        "surface and wl a comma list (per blob from early to late, per contour level from the lowest "
                        "to the highest, or over the colour range), or for surface and wl a transfer function such as "
                        "linear or sigmoid; layers: " + ", ".join(OPACITY_KEYS) + " (repeatable)")
    r.add_argument("--contours", action="store_true", help="draw the contours of the last radio image of every frequency")
    r.add_argument("--shell-color", choices=("mach_a", "mach_f", "theta_bn", "vn", "va"),
                   help="colour the GCS shell by the Alfven or fast-magnetosonic Mach number, theta_Bn, the normal speed, "
                        "or the upstream Alfven speed (corona3d shock; needs the apex speed of the fit or --shell-speed)")
    r.add_argument("--shell-clim", type=floats, metavar="LO,HI", help="colour limits of --shell-color")
    r.add_argument("--shell-speed", type=float, metavar="KM_S", help="apex speed for --shell-color when the configuration has no fit")
    r.add_argument("--radio-height-range", action="store_true",
                   help="draw the range of heights of every radio position over the density models of shock.height_models "
                        "as a segment along the line of sight")
    r.add_argument("--radio-surface", choices=("sphere", "model"),
                   help="place the radio sources on the sphere of the four-fold Newkirk height (default) or where the line of "
                        "sight meets the plasma-frequency surface of the scaled MAS density (config radio_deproject)")
    r.add_argument("--emission", choices=("H", "F", "auto"),
                   help="harmonic or fundamental emission of every radio source, or auto from the NRH polarisation "
                        "(config radio_emission and the emission key of each radio_images entry)")
    r.add_argument("--scattering", action="store_true",
                   help="correct the radio positions for the radial shift by scattering (config scattering.apply)")
    r.add_argument("--contour-levels", type=floats, metavar="L1,L2,...",
                   help="contour levels as fractions of the peak of the image, e.g. 0.3,0.5,0.7,0.9 (implies --contours)")
    r.add_argument("--panel", action="append", default=[], metavar="VIEW",
                   help="a view to draw, in this order (repeatable, or a comma list): earth, north, stereo, solo, or "
                        "Carrington degrees lon,lat[,roll] of the viewing direction; default: the three views of "
                        "config.json. A view that is not in the configuration takes the default framing")
    r.add_argument("--layout", choices=("landscape", "portrait", "separate"), default="landscape",
                   help="the panels in one row (landscape, default), one column (portrait), or one file per view "
                        "(separate: the name of the view is appended to --out)")


def render_params(a, cfg):
    """The configuration keys that the options of corona3d render and corona3d animate change."""
    params = {}
    if getattr(a, "wl_time", None):
        params["wl_time"] = a.wl_time
        params["wl_time_given"] = a.wl_time
    if getattr(a, "wl_at", None) == "burst":
        params.update(wl_at="burst", wl_height=a.wl_height, wl_height_note=a.wl_height_note)
    if getattr(a, "shell_on_top", False):
        params["shell_on_top"] = True
    if a.radio_sources:
        params["radio_sources"] = json.loads(a.radio_sources)
    if a.opacity:
        params["opacity"] = dict(parse_opacity(t) for t in a.opacity)
        if "wl" in params["opacity"]:
            params["wl_opacity"] = params["opacity"]["wl"]   # replaces a wl_opacity of config.json
    if a.shell_color:
        params["shell_color"] = a.shell_color
        if a.shell_clim:
            params["shell_clim"] = a.shell_clim
        if a.shell_speed is not None:
            params["shock_speed_km_s"] = a.shell_speed
    if a.radio_height_range:
        params["radio_height_range"] = True
    if a.radio_surface:
        params["radio_deproject"] = a.radio_surface
    if a.emission:
        params["radio_emission"] = a.emission
        params["radio_images"] = [{k: v for k, v in e.items() if k != "emission"} for e in cfg["radio_images"]]
    if a.scattering:
        params["scattering"] = {"apply": True}
    if a.contours or a.contour_levels:
        params["radio_images"] = [dict(e, contours=True, **({"contour_levels": a.contour_levels}
                                                            if a.contour_levels else {}))
                                  for e in params.get("radio_images", cfg["radio_images"])]
    from .reconstruction import parse_view
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
    return params


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
    f.add_argument("--model-extra", action="store_true",
                   help="only the temperature and flow files of the run (t002, vr002, vt002, vp002), used by corona3d shock")
    f.add_argument("--wl", action="store_true", help="only the LASCO C2 and STEREO COR1, COR2 images")
    f.add_argument("--cor1a", action="store_true", help="also the STEREO-A COR1 images (for corona3d gcs)")
    f.add_argument("--euv", action="store_true",
                   help="only the SDO/AIA 193 and 211 A and STEREO/EUVI 195 A images of euv_window (for corona3d gcs --euv)")
    f.add_argument("--solo", action="store_true",
                   help="only the Solar Orbiter images from the SOAR: Metis tB and pB (data/wl) and EUI/FSI 174 A (data/euv) "
                        "of wl_window, and the position of the spacecraft from their headers (config.json solo)")
    f.add_argument("--solohi", action="store_true", help="with --solo, also the SoloHI images (data/solohi)")
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
    g.add_argument("--t-ref", help="time of the six-parameter fit (default: chosen from the images; 2013-10-25T15:24:30 for the default event)")
    g.add_argument("--t-out", help="time of the reconstruction (default: time of config.json)")
    g.add_argument("--flare", nargs=2, type=float, metavar=("LON", "LAT"),
                   help="Stonyhurst longitude and latitude of the flare site [deg] that seeds the envelopes")
    g.add_argument("--views", default="c2,cor1a,cor1b", type=gcs_views,
                   help="instruments: c2, cor1a, cor1b, cor2a, cor2b, metis (Solar Orbiter; corona3d fetch --solo) (legacy stereo_a and stereo_b stand for "
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
    g.add_argument("--direction-cone", type=float, default=15.0, metavar="DEG",
                    help="with one view, or views whose observers are closer than 30 deg, the direction is not free: "
                         "bounded within DEG of the flare site (default 15; 0 fixes it at the flare site)")
    g.add_argument("--skip-step1", action="store_true", help="take step 1 from the gcs_fit.json of the working directory")
    g.add_argument("--euv", action="store_true",
                    help="measure the apex height at the reconstruction time in the EUV images of data/euv (corona3d fetch "
                         "--euv) with the shape of the white-light fit, in place of the extrapolated height (both are recorded)")
    g.add_argument("--euv-views", type=lambda t: t.split(","), metavar="LIST",
                    help="EUV views for --euv: aia193, aia211, euvi195a, euvi195b, fsi174 (Solar Orbiter/EUI) (default: aia193, "
                         "the EUVI and the FSI images present)")
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
    ev.add_argument("--t-ref", help="time of the six-parameter fit (default: chosen from the images, the earliest coronagraph "
                    "frame after --time at which the CME is well developed; the report states it). The radio and shell "
                    "time is always --time")
    ev.add_argument("--wl-at", choices=("fit", "burst"), default="fit",
                    help="time of the coronagraph figure: the fit time (default) or the burst time --time, with the shell "
                         "at t_out and each coronagraph's nearest frame after its base frame (labelled when more than 30 min away)")
    ev.add_argument("--strict", action="store_true",
                    help="exit with status 1 when any step failed (default: 0 for a partial result, with the failures "
                         "under Warnings in event_report.md; 1 only when no figure could be drawn)")
    ev.add_argument("--drop-poor-views", dest="drop_poor_views", action="store_true", default=True,
                    help="refit without a view whose IoU is far below the others and keep the refit when it is better "
                         "(default; rule in docs/METHODS.md)")
    ev.add_argument("--keep-all-views", dest="drop_poor_views", action="store_false", help="fit all views as given")
    ev.add_argument("--onset", help="flare onset (ISO time) that fixes the base window of the envelopes (default: the "
                    "reconstruction time)")
    ev.add_argument("--base-minutes", type=float, default=10.0, metavar="N",
                    help="the base frame is the last frame at or before onset minus N minutes (default 10), "
                         "so that the fit does not depend on when the download starts")
    ev.add_argument("--direction-cone", type=float, default=15.0, metavar="DEG",
                    help="with one view, or views whose observers are closer than 30 deg, the direction is not free: "
                         "bounded within DEG of the flare site (default 15; 0 fixes it at the flare site)")
    ev.add_argument("--radio-search", type=float, default=10.0, metavar="MINUTES",
                    help="radio images are searched within MINUTES of --time (default 10)")
    ev.add_argument("--euv", action="store_true",
                    help="measure the apex height at the reconstruction time in the EUV images of data/euv (corona3d fetch "
                         "--euv) with the shape of the white-light fit, in place of the extrapolated height (both are recorded)")
    ev.add_argument("--euv-views", type=lambda t: t.split(","), metavar="LIST",
                    help="EUV views for --euv: aia193, aia211, euvi195a, euvi195b, fsi174 (Solar Orbiter/EUI) (default: aia193, "
                         "the EUVI and the FSI images present)")
    ev.add_argument("--pb", action="store_true",
                    help="also fit the density scale to the polarised brightness (corona3d pb; section of event_report.md)")
    ev.add_argument("--sources", action="store_true",
                    help="also write the radio source positions under every placement "
                         "(corona3d sources)")
    ev.add_argument("--ovro-offset", action="store_true",
                    help="with --sources: also measure the OVRO-LWA quiet-Sun position offset (downloads about 18 image files)")
    ev.add_argument("--spectra", action="store_true",
                    help="also analyse the dynamic spectra of the type II burst (corona3d spectra)")
    ev.add_argument("--shock", action="store_true",
                    help="also run corona3d shock (Mach numbers, theta_Bn, speeds, radio height ranges, polarisation), "
                         "draw the shell coloured by the fast-magnetosonic Mach number, and add the section Shock "
                         "diagnostics to event_report.md")
    ev.add_argument("--force", action="store_true", help="repeat the downloads, the fit, and the figures")
    r = sub.add_parser("render", help="draw the three-panel figure",
                       epilog="layers: " + ", ".join(LAYERS))
    r.add_argument("--out", default="figures/corona3d.png", help="figure path, relative to the working directory")
    r.add_argument("--wl-time", help="white-light frame, e.g. 2013-10-25T15:24")
    r.add_argument("--wl-at", choices=("fit", "burst"), default="fit",
                   help="time of the coronagraph figure: the fit time wl_time (default) or the burst time, with the shell at "
                        "t_out and each coronagraph's nearest frame after its base frame, labelled when more than 30 min away")
    r.add_argument("--wl-height", type=float, help="apex height [R_sun] of the shell at the burst time (--wl-at burst), "
                   "for a configuration without gcs_fit; drawn as assumed")
    r.add_argument("--wl-height-note", help="reason for --wl-height (written on the figure)")
    figure_options(r)
    an = sub.add_parser("animate", help="GIF and MP4 of the reconstruction through a time interval",
                        description="One frame per time: the shell has the apex height of a height law and the radio "
                                    "images accumulate up to the frame time (docs/USER_GUIDE.md, Section 4.1). Without "
                                    "--start and --end the interval of the radio images of config.json is used.",
                        epilog="layers: " + ", ".join(LAYERS))
    an.add_argument("--out", default="figures/corona3d_animation", help="file name without extension, relative to the working directory")
    an.add_argument("--start", help="first frame, e.g. 2013-10-25T14:59:00")
    an.add_argument("--end", help="last frame")
    an.add_argument("--step", type=float, metavar="S", help="seconds between frames (default: step_s of the first radio entry, else 60)")
    an.add_argument("--fps", type=float, default=2.0, help="frames per second (default 2)")
    an.add_argument("--gif", action="store_true", help="write a GIF (default when --mp4 is not given)")
    an.add_argument("--mp4", action="store_true", help="write an MP4 (H.264, yuv420p; needs imageio-ffmpeg or ffmpeg)")
    an.add_argument("--gif-scale", type=float, default=0.5, metavar="F", help="size of the GIF frames relative to the figure (default 0.5)")
    an.add_argument("--height-fits", metavar="TIME=H,TIME=H",
                    help="apex height [R_sun] at two times, for example the two GCS fits; the shell rises on the straight line "
                         "through them (default: the shells gcs and gcs_wl of config.json)")
    an.add_argument("--speed", type=float, metavar="KM_S", help="constant apex speed from --height-at (default: the shell gcs at its time)")
    an.add_argument("--height-at", metavar="TIME=H", help="the height at one time, with --speed")
    an.add_argument("--heights", metavar="FILE", help="measured apex heights (EUV, Solar Orbiter) as a CSV file of time and height "
                    "[R_sun]; the frames must lie between the first and the last time")
    figure_options(an)
    sk = sub.add_parser("shock", help="Mach numbers, theta_Bn, and speeds of the CME shell, and heights, speeds, and "
                        "polarisation of the radio sources",
                        description="Samples the MAS/MAST corona on the GCS shell and at the radio sources and writes the "
                                    "tables and summary to the folder --out (method in docs/METHODS.md, Section 12).")
    sk.add_argument("--speed", type=float, metavar="KM_S", help="apex speed given by hand (default: the height-time fit of corona3d gcs)")
    sk.add_argument("--speed-note", help="reason for --speed, written into the summary")
    sk.add_argument("--out", default="analysis", help="output folder, relative to the working directory")
    sk.add_argument("--no-trace", action="store_true", help="do not trace the field line through every shell point (open or closed)")
    sk.add_argument("--no-radio", action="store_true", help="leave out the radio sources")
    sk.add_argument("--standoff", action="store_true",
                    help="also place a shock surface ahead of the shell at the stand-off distance of the local fast-mode Mach "
                         "number and compute the diagnostics on it (METHODS.md, Section 16)")
    sk.add_argument("--standoff-radius", choices=("min", "mean"),
                    help="radius of curvature of the stand-off: 1 / largest principal curvature (min, default) or 1 / mean curvature")
    sk.add_argument("--mc", type=int, metavar="N", help="ensemble of N members over apex height, speed, density scale, and model runs")
    sk.add_argument("--mc-seed", type=int, metavar="SEED", help="seed of the ensemble (default 0)")
    sk.add_argument("--no-connectivity", action="store_true", help="do not trace the field lines of the radio sources to the shock")
    sk.add_argument("--no-traps", action="store_true", help="do not search for closed loops crossed twice by the shock")
    sk.add_argument("--no-evolution", action="store_true", help="leave out the diagnostics at a series of times through the burst")
    sk.add_argument("--evolution-step", type=float, metavar="S", help="time step of the evolution in seconds (default 60)")
    pbp = sub.add_parser("pb", help="fit the density scale of the model to the polarised brightness of COR1 and C2",
                         description="Forward-models the Thomson-scattered polarised brightness of the MAS density, compares it "
                                     "with the calibrated pB of the coronagraphs, and writes the fitted scale to config.json as "
                                     "pb_fit, which density_scale \"pb\" uses (docs/METHODS.md, Section 17).")
    pbp.add_argument("--time", help="time of the polarisation sequences (default: the first of wl_window, before the CME)")
    pbp.add_argument("--out", default="analysis/pb", help="output folder, relative to the working directory")
    pbp.add_argument("--no-store", action="store_true", help="do not write pb_fit into config.json")
    so = sub.add_parser("sources", help="heights and positions of the radio sources on the Newkirk sphere and on the MAS "
                        "plasma surface, for F and H emission, with the scattering correction",
                        description="Writes analysis/radio_sources.csv (docs/METHODS.md, Sections 18 to 21).")
    so.add_argument("--out", default="analysis", help="output folder, relative to the working directory")
    so.add_argument("--ovro-offset", action="store_true",
                    help="measure the position offset of OVRO-LWA on the quiet Sun (downloads a few image files around the event)")
    spp = sub.add_parser("spectra", help="dynamic spectra of the type II burst (e-CALLISTO, OVRO-LWA), its lanes, shock speed, "
                         "and band-split Mach number", description="docs/METHODS.md, Section 22.")
    spp.add_argument("--out", default="analysis/spectra", help="output folder, relative to the working directory")
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
        everything = not (a.model or a.wl or a.ovro_lwa or a.nrh or a.model_extra or a.euv or a.solo)
        if a.model or everything:
            parts["model"] = fetch.model(h, cfg)
        if a.model_extra:
            parts["model extra"] = fetch.model(h, cfg, fetch.EXTRA)
        if a.euv:
            from . import euv
            parts["EUV"] = euv.fetch(h, cfg)
        if a.solo:
            from . import solo
            parts["Solar Orbiter"] = solo.fetch(h, cfg, solohi=a.solohi)
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
    if a.cmd == "shock":
        from . import shock
        if a.speed is not None and not a.speed_note:
            sys.exit("--speed needs --speed-note, the reason for the assumed speed")
        if fetch.model(h, cfg, fetch.EXTRA):
            print("some temperature or flow files failed; the flow is then not subtracted, and without t002 the "
                  "sound and fast-magnetosonic speeds need shock.temperature_k in config.json")
        shock.run(h, cfg, a)
        return
    if a.cmd in ("pb", "sources", "spectra"):
        from importlib import import_module
        import_module(f"corona3d.{a.cmd}").run(h, cfg, a)
        return
    if a.cmd == "export":
        from .exporting import export
        params = {"radio_sources": json.loads(a.radio_sources)} if a.radio_sources else None
        d = export(h, a.out, a.format, params=params,
                   **{k: getattr(a, k) for k in LAYERS if k != "sphere"})
        print("exported to", d)
        return
    params = render_params(a, cfg)
    panels = a.panel or None
    if a.cmd == "animate":
        from . import animate
        from .animate import height_law, parse_height_points, read_heights
        try:
            pm = paths.merge(cfg, params)
            law = height_law(pm, fits=parse_height_points(a.height_fits) if a.height_fits else None,
                             table=read_heights(a.heights) if a.heights else None, speed_km_s=a.speed,
                             anchor=parse_height_points(a.height_at)[0] if a.height_at else None)
            animate.animate(a.out, h, a.start, a.end, a.step, a.fps, a.gif, a.mp4, a.half, params, panels, a.layout,
                            law, a.gif_scale, **{k: getattr(a, k) for k in LAYERS})
        except (ValueError, RuntimeError) as e:
            sys.exit(f"corona3d animate: {e}")
        return
    from .reconstruction import make_figure
    make_figure(a.out, h, a.half, params, panels=panels, layout=a.layout, **{k: getattr(a, k) for k in LAYERS})


if __name__ == "__main__":
    main()
