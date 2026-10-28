"""`corona3d event`: from a time and a flare position to the figures and a report, in one command.

The steps are new (config.json), fetch (model, white light, OVRO-LWA, NRH check), radio-interval
(burst interval from the images), gcs (free fit, falling back to the fixed-direction fit), render
(base, coronagraph layers, radio layer, each when its data exist), and the report
`event_report.md`. Each step skips work that is already done: the file `event_state.json` of the
working directory records the downloads and the fit, and a figure is drawn again only when it is
older than `config.json`. `--force` repeats everything except writing config.json.

Times: the reconstruction time t_out (the shell, the radio interval search, the shell at the radio
sources) is always `--time`, which is the begin time of the burst when `--burst` is used. The time of
the six-parameter fit t_ref is `--t-ref` when given; otherwise it is chosen from the images as the
earliest coronagraph frame after t_out at which the CME is well developed (gcs.pick_t_ref), and the
report states it. The white-light figure uses wl_time, which follows t_ref.

Exit code: 0 when at least one figure exists, with every failed step listed under "Warnings" in
event_report.md; 1 only when no figure could be drawn, or with `--strict` when any step failed.
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
    """{instrument: text} for the instruments that cannot give a white-light panel at wl_time: no frame
    within whitelight.WL_MAX_GAP_MIN of it (a coronagraph gap, e.g. LASCO C2 on 2024-05-09 from 14:48 to
    18:52 UT), or the nearest frame inside wl_window is the first frame of the window, which is the
    base frame itself and so has nothing to be compared with. The figure is drawn with the others."""
    from . import whitelight as wlt
    t = wlt.utc(cfg.get("wl_time") or cfg["time"])
    win = [wlt.utc(x) for x in cfg["wl_window"]] if cfg.get("wl_window") else None
    out = {}
    for k, ok in wl_on.items():
        if not ok:
            continue
        fr = wlt._frames(wl, k, sc.upper())
        use = [f for f in fr if win[0] <= f[0] <= win[1]] if win else fr
        near = min(use, key=lambda f: abs((f[0] - t).total_seconds())) if use else None
        if near is None or abs((near[0] - t).total_seconds()) > 60 * wlt.WL_MAX_GAP_MIN:
            out[k] = f"no frame within {wlt.WL_MAX_GAP_MIN:g} min of wl_time {t:%Y-%m-%dT%H:%M} ({wlt.nearest_times(fr, t)})"
        elif win and near[0] == use[0][0]:
            out[k] = (f"the frame nearest to wl_time {t:%Y-%m-%dT%H:%M} ({near[0]:%H:%M:%S}) is the first frame of "
                      "wl_window, the base frame, so there is no later frame to show")
    return out


def wl_nearest(wl, cfg, wl_on, sc):
    """(time, instrument) of the frame closest to wl_time that can be drawn (inside wl_window and not the
    first frame of it), over the instruments with images; None when there is none."""
    from . import whitelight as wlt
    t = wlt.utc(cfg.get("wl_time") or cfg["time"])
    win = [wlt.utc(x) for x in cfg["wl_window"]] if cfg.get("wl_window") else None
    best = None
    for k, ok in wl_on.items():
        if not ok:
            continue
        fr = wlt._frames(wl, k, sc.upper())
        use = [f for f in fr if win[0] <= f[0] <= win[1]] if win else fr
        for f in use[1:]:
            d = abs((f[0] - t).total_seconds())
            if best is None or d < best[0]:
                best = (d, f[0], k)
    return (best[1], best[2]) if best else None


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
    """Fit the shell. With fit "auto": (1) the free fit on all views, (2) the view rule inside that fit
    (a poor view is left out and the free fit repeated), (3) the quality of the result (gcs.quality),
    (4) only when it is poor or failed, the fixed-direction fit, also with the view rule, and (5) the
    better of the two (gcs.choose_fit). With "fixed" only the fixed-direction fit. Writes gcs_fit.json
    (with the key `decision`) and config.json for the accepted result. Returns the record kept in the
    state file."""
    wl = paths.wl_dir(h)
    t_ref, t_out = a.t_ref, cfg["time"]              # t_ref None: chosen from the images by gcs.fit
    attempts, results = [], {}

    def attempt(fixed):
        name = "fixed direction" if fixed else "free"
        try:
            out = gcs.fit(wl, t_ref, t_out, flare, views, fixed, None, None, None, a.height_out, a.height_note,
                          drop_poor_views=getattr(a, "drop_poor_views", True), iou_min=a.iou_min,
                          onset=getattr(a, "onset", None), base_minutes=getattr(a, "base_minutes", gcs.BASE_MIN),
                          direction_cone=getattr(a, "direction_cone", gcs.CONE))
        except (SystemExit, Exception) as e:                     # SystemExit: the fit found no envelope or front
            attempts.append(dict(mode=name, status="failed", why=[str(e) or type(e).__name__]))
            log(f"gcs: {name} fit failed ({attempts[-1]['why'][0]})")
            return None
        why = gcs.quality(out, flare, fixed, a.iou_min)
        m = gcs.mean_iou(out)
        attempts.append(dict(mode=name, status="poor" if why else "good", why=why, view_check=out.get("view_check"),
                             views=out.get("views_used"), direction=out.get("direction"), dropped=(out.get("view_check") or {}).get("dropped", []),
                             bounds=[h[0] for h in gcs.bound_hits(out["t_ref"], flare, fixed or bool((out.get("direction") or {}).get("constrained")),
                                                         (out.get("height_bound") or {}).get("max", gcs.HEIGHT_MAX))],
                             height_bound=out.get("height_bound"), mean_iou=None if m is None else round(m, 3)))
        log(f"gcs: {name} fit (views {', '.join(out.get('views_used') or [])}) "
            f"{'poor: ' + '; '.join(why) if why else 'accepted'}")
        results[name] = out
        return attempts[-1]

    free = attempt(False) if a.fit == "auto" else None
    fixed = attempt(True) if a.fit != "auto" or free is None or free["status"] != "good" else None
    pick = gcs.choose_fit(free, fixed)
    rec = dict(attempts=attempts, flare=list(flare), views=list(views), t_ref=t_ref, t_out=t_out,
               height_out=a.height_out, iou_min=a.iou_min, accepted=None)
    if pick:
        out = results[pick]
        out["decision"] = dict(
            order=["free fit on all views", "view rule and refit of the free fit", "quality of the free fit",
                   "fixed-direction fit with the view rule, only if the free fit is poor or failed", "better of the two"],
            criterion=f"a good fit beats a poor one; otherwise the higher mean IoU of the kept views wins, and a difference "
                      f"below {gcs.TIE_IOU:g} goes to the free fit",
            fixed_run=fixed is not None, attempts=[dict(x, view_check=None) for x in attempts], accepted=pick)
        gcs.write_result(h, out, flare, views, pick == "fixed direction")
        rec["accepted"] = pick
        rec["t_ref"] = out["t_ref"]["time"]
        rec["t_ref_choice"] = out.get("t_ref_choice") or dict(method="--t-ref")
        if fixed is None and a.fit == "auto":
            log("gcs: the free fit is good; the fixed-direction fit was not needed")
    return rec


def _same(a, b):
    return a is not None and b is not None and abs((_iso(a) - _iso(b)).total_seconds()) < 1


def run(h, a, log=print):
    """The `corona3d event` command; returns the exit code: 0 for a complete or partial result (failed steps
    are listed under Warnings in event_report.md), 1 when no figure exists or, with a.strict, when any step failed."""
    from ._extras import need
    from .events import new_config, parse_flare
    need("render", "pyvista", "the figures of `corona3d event`")      # before any download
    need("render", "PyThea.geometrical_models", "the CME shell of the figures")
    if a.fit != "skip":
        need("gcs", "skimage.morphology", "the GCS fit (--fit skip leaves it out)")
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
    elif (st.get("gcs") or {}).get("accepted") and cfg.get("wl_time") != (st.get("wl_moved") or {}).get("new"):
        tr = st["gcs"]["t_ref"]
        gcs.set_wl_time(h, tr, "t_ref chosen automatically")
        step("wl_time", f"set to {tr} (the t_ref of the fit)")
    _save_state(h, st)
    cfg = paths.config(h)

    # render
    from . import reconstruction
    figures, wl_on = {}, {k: avail[k] for k in ("c2", "cor1", "cor2")}
    st["wl_gap"] = gaps = wl_gaps(paths.wl_dir(h), cfg, wl_on, sc)
    if wl_on and any(wl_on.values()) and all(k in gaps for k, ok in wl_on.items() if ok):
        near = wl_nearest(paths.wl_dir(h), cfg, wl_on, sc)           # no instrument has a usable frame at wl_time
        if near:
            old = cfg.get("wl_time") or cfg["time"]
            gcs.set_wl_time(h, f"{near[0]:%Y-%m-%dT%H:%M:%S}", "nearest frame (no instrument had a frame at wl_time)")
            cfg = paths.config(h)
            st["wl_moved"] = dict(old=old, new=cfg["wl_time"], instrument=near[1])
            step("wl_time", f"moved from {old} to {cfg['wl_time']}, the nearest {near[1]} frame")
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
    st["failed"] = failed
    _save_state(h, st)
    if not figures:
        return 1                                  # a total failure: nothing could be drawn
    return 1 if failed and getattr(a, "strict", False) else 0


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
        how = g.get("t_ref_choice") or {}
        if how.get("method") == "auto":
            measured.append(f"The fit time {g['t_ref']} was chosen from the images: {how['rule']}; it is developed in "
                            f"{', '.join(how['developed'])}. The radio and shell time is {g['t_out']}, the time given with --time.")
        elif how.get("method") == "fallback":
            assumed.append(f"No coronagraph frame showed a developed CME ({how['why']}), so the fit time {g['t_ref']} is "
                           f"the reconstruction time plus 60 minutes. The radio and shell time is {g['t_out']}.")
        elif a.t_ref:
            measured.append(f"The fit time {g['t_ref']} was given with --t-ref. The radio and shell time is {g['t_out']}.")
        if g.get("height_source") == "measured":
            measured.append(f"The apex height at {g['t_out']} is {g['height_t_out']:g} +/- {g['height_t_out_err']:g} R_sun, "
                            f"extrapolated from the front in {g['front']} ({g['n_frames']} frames, {g['speed_km_s']:g} km/s).")
        else:
            assumed.append(f"The apex height at {g['t_out']} is {g['height_t_out']:g} R_sun, given by hand: {g.get('height_note', '')}")
        dc = g.get("direction_constraint")
        if dc:
            where = (f"held at the flare site (Stonyhurst longitude {flare[0]:g} deg, latitude {flare[1]:g} deg)" if dc["cone_deg"] <= 0
                     else f"limited to {dc['cone_deg']:g} deg around the flare site (Stonyhurst longitude {flare[0]:g} deg, "
                          f"latitude {flare[1]:g} deg)")
            assumed.append(f"{dc['text'][0].upper() + dc['text'][1:]}: the shell axis is {where}, because a fit to "
                           f"{', '.join(dc['views'])} cannot determine the direction along the line of sight "
                           f"(Thernisien et al. 2011; option --direction-cone).")
        elif g["fixed"]:
            assumed.append(f"The shell axis is held at the flare site (Stonyhurst longitude {flare[0]:g} deg, latitude "
                           f"{flare[1]:g} deg), "
                           + ("because the fixed-direction fit was requested." if a.fit == "fixed" else
                              "because the free fit was poor or failed and the fixed-direction fit was the better one (see the fit quality below)."))
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
    L += [f"- {m}" for m in missing] if missing else ["- None."]
    L += ["", "## Warnings", ""]
    L += [f"- Failed: {x}." for x in failed] or ["- None."]
    if not (st.get("figures") or {}):
        L.append("- No figure could be drawn (total failure).")
    if st.get("wl_moved"):
        w = st["wl_moved"]
        L.append(f"- No coronagraph had a frame usable at wl_time {w['old']}; wl_time was moved to {w['new']}, the nearest "
                 f"{w['instrument']} frame, so the coronagraph panel shows that frame and not the fit time.")
    L += ["", "## Fit quality", ""]
    if fit.get("attempts"):
        if a.fit == "auto":
            L.append("- Order: (1) free fit on all views, (2) view rule and refit of the free fit, (3) quality of the free fit, "
                     "(4) fixed-direction fit with the view rule only if the free fit is poor or failed, (5) the better of the two.")
        for t in fit["attempts"]:
            used = f"; views used: {', '.join(t['views'])}" if t.get("views") else ""
            L.append(f"- {t['mode']} fit: {t['status']}" + (f" ({'; '.join(t['why'])})" if t["why"] else "") + used + ".")
            if t.get("direction"):
                L.append(f"  - {t['direction']['text']}; cone {t['direction']['cone_deg']:g} deg around the flare site.")
            L += [f"  - {x}" for x in gcs.view_lines(t.get("view_check"))]
        if a.fit == "auto":
            names = [t["mode"] for t in fit["attempts"]]
            if "fixed direction" not in names:
                L.append("- The fixed-direction fit was not run, because the free fit was good.")
            L.append(f"- Kept: the {fit['accepted']} fit" if fit.get("accepted") else "- Kept: none, no fit succeeded")
            if fit.get("accepted") and "fixed direction" in names:
                L[-1] += f", the better of the two (a good fit beats a poor one; otherwise the higher mean IoU of the kept views, a difference below {gcs.TIE_IOU:g} going to the free fit)"
            L[-1] += "."
        acc = next((t for t in fit["attempts"] if t["mode"] == fit.get("accepted")), None)
        hb = (acc or {}).get("height_bound")
        if hb and hb.get("outer_fov"):
            L.append(f"- The apex height of the fit was bounded above by {hb['max']:g} R_sun, {hb['factor']:g} times the outer field of view "
                     f"({hb['outer_fov']:g} R_sun) of the fitted views (never above {hb['cap']:g} R_sun), because the images cannot constrain an apex far outside the field.")
        hit = [b for b in (acc or {}).get("bounds", []) if b in ("kappa", "alpha_deg", "height")]
        if hit:
            L.append(f"- The kept fit has {' and '.join(hit)} at its bound, so the images do not fix the shape of the shell "
                     "in that parameter and the fit stays poor.")
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
