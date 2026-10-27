"""`corona3d event`: from a time and a flare position to the figures and a report, in one command.

The steps are new (config.json), fetch (model, white light, OVRO-LWA, NRH check), radio-interval
(burst interval from the images), gcs (free fit, falling back to the fixed-direction fit), render
(base, coronagraph layers, radio layer, each when its data exist), and the report
`event_report.md`. Each step skips work that is already done: the file `event_state.json` of the
working directory records the downloads and the fit, and a figure is drawn again only when it is
older than `config.json`. `--force` repeats everything except writing config.json.
"""
import json
from datetime import datetime, timedelta, timezone
from pathlib import Path
from types import SimpleNamespace

from . import __version__, fetch, gcs, paths, radio
from .layers import LAYERS, ON

STATE = "event_state.json"
REPORT = "event_report.md"
FIGURES = {"base": "figures/corona3d.png", "wl": "figures/wl/corona3d.png", "radio": "figures/radio/corona3d.png"}


def wl_available(wl, sc):
    """{'c2', 'cor1', 'cor2': bool} for the white-light images in a folder (STEREO spacecraft sc)."""
    d, sc = Path(wl), sc.upper()
    return dict(c2=any(f.stem.isdigit() for f in d.glob("2*.fts")), cor1=any(d.glob(f"*_[sn][45]c1{sc}.fts")),
                cor2=any(d.glob(f"*_d4c2{sc}.fts"))) if d.exists() else dict(c2=False, cor1=False, cor2=False)


def wl_gaps(wl, cfg, wl_on, sc):
    """{instrument: text} for the instruments with images that have no frame within
    whitelight.WL_MAX_GAP_MIN of the white-light time (a coronagraph gap, e.g. LASCO C2 on 2024-05-09
    from 14:48 to 18:52 UT); the figure is drawn with the other instruments."""
    from . import whitelight as wlt
    t = wlt.utc(cfg.get("wl_time") or cfg["time"])
    out = {}
    for k, ok in wl_on.items():
        if not ok:
            continue
        fr = wlt._frames(wl, k, sc.upper())
        if not fr or min(abs((f[0] - t).total_seconds()) for f in fr) > 60 * wlt.WL_MAX_GAP_MIN:
            out[k] = f"no frame within {wlt.WL_MAX_GAP_MIN:g} min of wl_time {t:%Y-%m-%dT%H:%M} ({wlt.nearest_times(fr, t)})"
    return out


def radio_available(h, cfg):
    """[(entry index, instrument, folder has files)] of the radio_images entries."""
    return [(i, rs["instrument"], (h / rs["folder"]).exists() and any((h / rs["folder"]).iterdir()))
            for i, rs in enumerate(cfg["radio_images"])]


def _iso(t):
    return datetime.fromisoformat(t) if isinstance(t, str) else t


def _load_state(h):
    f = h / STATE
    return json.loads(f.read_text(encoding="utf-8")) if f.exists() else {}


def _save_state(h, st):
    (h / STATE).write_text(json.dumps(st, indent=1), encoding="utf-8")


def fit_step(h, cfg, a, flare, views, log):
    """Fit the shell: with fit "auto" the free fit first, and the fixed-direction fit when it fails or
    is poor (see gcs.quality); with "fixed" only the latter. Writes gcs_fit.json and config.json for
    the accepted result. Returns the record kept in the state file."""
    wl = paths.wl_dir(h)
    t_ref, t_out = a.t_ref or cfg["gcs_wl"]["time"], cfg["time"]
    modes = [False, True] if a.fit == "auto" else [True]
    attempts, accepted = [], None
    for fixed in modes:
        name = "fixed direction" if fixed else "free"
        try:
            out = gcs.fit(wl, t_ref, t_out, flare, views, fixed, None, None, None, a.height_out, a.height_note)
        except (SystemExit, Exception) as e:                     # SystemExit: the fit found no envelope or front
            attempts.append(dict(mode=name, status="failed", why=[str(e) or type(e).__name__]))
            log(f"gcs: {name} fit failed ({attempts[-1]['why'][0]})")
            continue
        why = gcs.quality(out, flare, fixed, a.iou_min)
        attempts.append(dict(mode=name, status="poor" if why else "good", why=why))
        log(f"gcs: {name} fit {'poor: ' + '; '.join(why) if why else 'accepted'}")
        if accepted is None or not why or fixed:                  # the fixed fit replaces a poor free fit
            accepted = (out, fixed)
        if not why:
            break
    rec = dict(attempts=attempts, flare=list(flare), views=list(views), t_ref=t_ref, t_out=t_out,
               height_out=a.height_out, iou_min=a.iou_min, accepted=None)
    if accepted:
        out, fixed = accepted
        gcs.write_result(h, out, flare, views, fixed)
        rec["accepted"] = "fixed direction" if fixed else "free"
    return rec


def _same(a, b):
    return a is not None and b is not None and abs((_iso(a) - _iso(b)).total_seconds()) < 1


def run(h, a, log=print):
    """The `corona3d event` command; returns the exit code (1 when a download, the fit, or a figure failed)."""
    from .events import new_config, parse_flare
    h = Path(h)
    h.mkdir(parents=True, exist_ok=True)
    steps, failed = [], []

    def step(name, status):
        steps.append((name, status))
        log(f"[{name}] {status}")

    # new
    f = h / "config.json"
    if f.exists() and not _same(paths.config(h)["time"], a.time):
        raise SystemExit(f"{f} is for {paths.config(h)['time']}, not {a.time}; use another --home, or "
                         "corona3d new --force there")
    if f.exists():
        step("new", "config.json present; kept")
    else:
        cfg = paths.merge(paths.default_config(), new_config(a.time, a.flare))
        f.write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")
        step("new", f"wrote config.json ({cfg['run']}, STEREO-{cfg['stereo_sc']})")
    cfg = paths.config(h)
    flare = tuple(parse_flare(a.flare))
    st = {} if a.force else _load_state(h)
    done = st.setdefault("done", {})
    sc = cfg["stereo_sc"]
    views = tuple(a.views.split(",")) if a.views else ("c2", f"cor1{sc.lower()}", f"cor2{sc.lower()}")

    # fetch
    n = fetch.model(h, cfg)
    step("fetch_model", f"{n} files failed" if n else "model files present")
    if n:
        failed.append("model download")
    if done.get("fetch_wl"):
        step("fetch_wl", "done earlier; skipped")
    else:
        try:
            n = fetch.whitelight(h, cfg, cor1a=("cor1a" in views and sc != "A"))
        except Exception as e:
            n = 1
            log(f"white-light download: {e}")
        step("fetch_wl", f"{n} files failed" if n else "downloaded")
        if n:
            failed.append("white-light download")
        else:
            done["fetch_wl"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    # the OVRO-LWA files of the search window, for the entries whose interval is chosen from the images
    t = _iso(cfg["time"])
    wide = dict(cfg, radio_images=[
        dict(e, start=(t - timedelta(minutes=a.radio_search)).isoformat(),
             end=(t + timedelta(minutes=a.radio_search)).isoformat()) if e.get("interval_source") in ("default", "data") else e
        for e in cfg["radio_images"]])
    if any(e["instrument"] == "ovro_lwa" for e in cfg["radio_images"]):
        key = [[e["start"], e["end"]] for e in wide["radio_images"] if e["instrument"] == "ovro_lwa"]
        if done.get("fetch_ovro_lwa") == key:
            step("fetch_ovro_lwa", "done earlier; skipped")
        else:
            try:
                n = fetch.ovro_lwa(h, wide)
            except Exception as e:
                n = 1
                log(f"OVRO-LWA download: {e}")
            step("fetch_ovro_lwa", f"{n} files failed" if n else "downloaded")
            if n:
                failed.append("OVRO-LWA download")
            else:
                done["fetch_ovro_lwa"] = key
    else:
        step("fetch_ovro_lwa", "no OVRO-LWA entry (before 2024 April); nothing to fetch")
    if any(e["instrument"] == "nrh" for e in cfg["radio_images"]):
        n = fetch.nrh(h, cfg)
        step("fetch_nrh", f"{n} frequencies without files; no automatic download exists" if n
             else "NRH files present")
    _save_state(h, st)

    # radio interval
    notes = radio.select(h, cfg, a.radio_search)
    st["radio"] = [m for _, m in notes]
    for _, m in notes:
        log(m)
    step("radio_interval", f"{sum('one image every' in m for m in st['radio'])} entries set from the images")
    cfg = paths.config(h)

    # gcs
    avail = wl_available(paths.wl_dir(h), sc)
    views = tuple(k for k in views if avail.get(k[:4]))
    old = cfg.get("gcs_fit") or {}
    if a.fit == "skip":
        step("gcs", "--fit skip; the shell of config.json is used")
    elif not views:
        step("gcs", "no coronagraph images; the shell of config.json is used")
    elif (old and not a.force and list(old.get("flare", [])) == list(flare) and _same(old.get("t_out"), cfg["time"])
          and (a.t_ref is None or _same(old.get("t_ref"), a.t_ref))
          and (a.height_out is None or old.get("height_t_out") == a.height_out)):
        step("gcs", "fit of this event present; skipped")
    else:
        st["gcs"] = fit_step(h, cfg, a, flare, views, log)
        if st["gcs"]["accepted"]:
            step("gcs", f"{st['gcs']['accepted']} fit written")
        else:
            step("gcs", "no fit succeeded; the shell of config.json is used")
            failed.append("GCS fit")
    if a.t_ref:                                   # the white-light figure uses frames of the time of the fit
        gcs.set_wl_time(h, a.t_ref)
        step("wl_time", f"set to {a.t_ref} (--t-ref)")
    _save_state(h, st)
    cfg = paths.config(h)

    # render
    from . import reconstruction
    figures, wl_on = {}, {k: avail[k] for k in ("c2", "cor1", "cor2")}
    st["wl_gap"] = gaps = wl_gaps(paths.wl_dir(h), cfg, wl_on, sc)
    wl_on = {k: ok and k not in gaps for k, ok in wl_on.items()}
    model_ok = all((paths.model_dir(h, cfg) / f"{v}.hdf").exists() for v in fetch.VARS)
    radio_on = any(ok for _, _, ok in radio_available(h, cfg))
    cases = {"base": {}, "wl": wl_on if any(wl_on.values()) else None, "radio": {"radio": True} if radio_on else None}
    cfg_time = (h / "config.json").stat().st_mtime
    for name, layers in cases.items():
        out = h / FIGURES[name]
        if layers is None:
            step(f"render_{name}", "no data for this figure; skipped")
        elif not model_ok:
            step(f"render_{name}", "model files missing; skipped")
        elif out.exists() and out.stat().st_mtime >= cfg_time and not a.force:
            figures[name] = out
            step(f"render_{name}", "figure newer than config.json; skipped")
        else:
            try:
                reconstruction.make_figure(FIGURES[name], h, None, None, **dict(ON, **layers))
                figures[name] = out
                step(f"render_{name}", "drawn")
            except Exception as e:
                log(f"render {name}: {e}")
                step(f"render_{name}", f"failed ({e})")
                failed.append(f"figure {name} ({e})")
    st["figures"] = {k: str(v.relative_to(h)).replace("\\", "/") for k, v in figures.items()}
    _save_state(h, st)

    (h / REPORT).write_text(report(h, cfg, a, flare, views, avail, st, steps, failed), encoding="utf-8")
    step("report", f"wrote {REPORT}")
    return 1 if failed else 0


def report(h, cfg, a, flare, views, avail, st, steps, failed):
    """Text of event_report.md: what was measured, what was assumed and why, what data were missing,
    the fit quality per view, and the figures."""
    from .fetch import VARS
    g = cfg.get("gcs_fit") or {}
    fit = st.get("gcs") or {}
    sc = cfg["stereo_sc"]
    L = [f"# Event report: {cfg['time']} UT, flare {a.flare.upper()}", "",
         f"Written by corona3d {__version__} on {datetime.now(timezone.utc):%Y-%m-%d %H:%M} UTC. Model run {cfg['run']}, "
         f"STEREO-{sc}, Earth at Carrington longitude {cfg['earth'][0]:g} deg and latitude {cfg['earth'][1]:g} deg.", "",
         "## Measured", ""]
    measured, assumed, missing = [], [], []
    if g:
        mode = "fixed direction" if g["fixed"] else "free"
        s, w = cfg["gcs"], cfg["gcs_wl"]
        measured.append(f"The CME shell was fitted to the coronagraph images of {', '.join(g['views'])} at {g['t_ref']} "
                        f"({mode} fit): height {w['height']:g} R_sun, aspect ratio {w['kappa']:g}, half-angle "
                        f"{w['alpha_deg']:g} deg, direction Stonyhurst longitude {w['lon_stonyhurst']:g} deg and latitude "
                        f"{w['lat_deg']:g} deg.")
        if g.get("height_source") == "measured":
            measured.append(f"The apex height at {g['t_out']} is {g['height_t_out']:g} +/- {g['height_t_out_err']:g} R_sun, "
                            f"extrapolated from the front in {g['front']} ({g['n_frames']} frames, {g['speed_km_s']:g} km/s).")
        else:
            assumed.append(f"The apex height at {g['t_out']} is {g['height_t_out']:g} R_sun, given by hand: {g.get('height_note', '')}")
        if g["fixed"]:
            assumed.append(f"The shell axis is held at the flare site (Stonyhurst longitude {flare[0]:g} deg, latitude "
                           f"{flare[1]:g} deg), because the free fit "
                           + ("was not requested." if a.fit == "fixed" else "failed or was poor (see the fit quality below)."))
    else:
        assumed.append("The CME shell is the default one placed above the flare site (aspect ratio 0.4, half-angle 30 deg, "
                       f"height {cfg['gcs']['height']:g} R_sun at the time of the event and {cfg['gcs_wl']['height']:g} R_sun "
                       "at the image time); no fit was made or none succeeded, so its size and height are not measured.")
    radio_notes = [m for m in st.get("radio", []) if "MHz:" in m and "no image" not in m]
    measured += [f"Radio burst interval: {m}." for m in radio_notes]
    if not measured:
        measured.append("Nothing was measured from data.")
    L += [f"- {m}" for m in measured] + ["", "## Assumed", ""]
    assumed += [f"The corona is the steady-state MAS/MAST model run {cfg['run']} of the Carrington rotation of the event, "
                "not a model of the corona of that moment.",
                "Radio sources are placed at the height of the second harmonic in a four-fold Newkirk density model "
                "on the side of the plane of the sky that faces Earth."]
    for i, rs in enumerate(cfg["radio_images"]):
        if rs.get("interval_source") == "default":
            assumed.append(f"The {rs['instrument']} interval {rs['start']} to {rs['end']} is the default window after the "
                           "event time, because no image reached min_tb in the search window or no image was available.")
        elif "interval_source" not in rs:
            assumed.append(f"The {rs['instrument']} interval {rs['start']} to {rs['end']} was given in config.json.")
    L += [f"- {m}" for m in assumed] + ["", "## Missing data", ""]
    d = paths.model_dir(h, cfg)
    miss = [v for v in VARS if not (d / f"{v}.hdf").exists()]
    if miss:
        missing.append(f"Model files not downloaded: {', '.join(miss)}.")
    inst = {"c2": "LASCO C2", "cor1": f"STEREO-{sc} COR1", "cor2": f"STEREO-{sc} COR2"}
    missing += [f"No {inst[k]} images in data/wl." for k, ok in avail.items() if not ok]
    missing += [f"{inst[k]}: {m}; left out of the white-light figure." for k, m in (st.get("wl_gap") or {}).items()]
    if a.fit != "skip" and not views:
        missing.append("No coronagraph images for the GCS fit.")
    for i, instr, ok in radio_available(h, cfg):
        if not ok:
            how = " (the NRH 2-D images must be requested from https://rsdb.obs-nancay.fr/ and copied there)" if instr == "nrh" else ""
            missing.append(f"No {instr} files in {cfg['radio_images'][i]['folder']}{how}.")
    missing += [m for m in st.get("radio", []) if "no image at or above" in m or "no frequency reaches" in m]
    missing += [f"Failed: {x}." for x in failed]
    L += [f"- {m}" for m in missing] if missing else ["- None."]
    L += ["", "## Fit quality", ""]
    if fit.get("attempts"):
        for t in fit["attempts"]:
            L.append(f"- {t['mode']} fit: {t['status']}" + (f" ({'; '.join(t['why'])})" if t["why"] else "") + ".")
        if g.get("iou"):
            L.append(f"- Intersection over union of the accepted shell with the envelope, per view at {g['t_ref']}: "
                     + ", ".join(f"{k} {v:g}" for k, v in g["iou"].items())
                     + f" (a fit is poor below a mean of {fit.get('iou_min', gcs.IOU_MIN):g}).")
    else:
        L.append("- No fit was made in this run.")
    L += ["", "## Figures", ""]
    figs = st.get("figures", {})
    L += [f"- {FIGURES[k]}: " + {"base": "model, shell, and field lines", "wl": "with the coronagraph images",
                                  "radio": "with the radio sources"}[k] for k in figs] or ["- None."]
    L += ["", "## Steps of the last run", ""] + [f"- {n}: {s}" for n, s in steps]
    return "\n".join(L) + "\n"
