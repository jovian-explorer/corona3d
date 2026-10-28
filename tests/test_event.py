"""`corona3d event` and `corona3d radio-interval`: the order of the steps, work that is skipped when
done, the fallback from the free to the fixed-direction fit, the content of the report, and the
choice of the burst interval from synthetic radio images. The network, the fit, and the rendering
are replaced by mocks; the radio images are real (synthetic) NRH files."""
import json
from datetime import datetime, timedelta
from pathlib import Path
from types import SimpleNamespace

import pytest
from test_nrh import T0, write_nrh

from corona3d import events, fetch, gcs, pipeline, radio, reconstruction
from corona3d.events import new_config
from corona3d.paths import config, default_config, merge, model_dir

FLARE = (-69.0, -6.0)
ORIGINAL_WL_GAPS = pipeline.wl_gaps
NRH_FILE = "nrh2_1509_h60_20131025_145800c00_i.fts"


def args(**kw):
    d = dict(time="2013-10-25T15:00", flare="S06E69", fit="auto", iou_min=0.5, height_out=None, height_note=None,
             views=None, t_ref=None, radio_search=10.0, force=False)
    return SimpleNamespace(**dict(d, **kw))


def result(iou=0.8, kappa=0.4, fixed=False):
    """A result of gcs.fit with the keys that gcs.config_update and the report read."""
    shell = dict(height=5.5, kappa=kappa, alpha_deg=25.0, lon_stonyhurst=-69.0 if fixed else -60.0,
                 lat_deg=-6.0 if fixed else -2.0, tilt_pythea_deg=-30.0)
    return {"t_ref": dict(shell, time="2013-10-25T16:00:00", iou={"c2": iou, "cor1b": iou}),
            "t_out": dict(shell, height=1.9, time="2013-10-25T15:00:00"),
            "height_time": dict(instrument="cor1b", projection=0.9, frames=[]),
            "fit": dict(speed_rsun_per_min=0.09, speed_km_s=1000.0, height_t_out=1.9, height_t_out_err=0.1, n_frames=6)}


@pytest.fixture
def world(tmp_path, monkeypatch):
    """Mocks of every network and heavy step; `calls` lists them in the order they run."""
    calls, fits, fit_results = [], [], []
    monkeypatch.setattr(events, "observers", lambda t, sc: ([357.4, 5.0], [10.0, 1.0]))
    monkeypatch.setattr(events, "psi_run", lambda t: "cr2143-medium/hmi_mast_mas_std_0201")

    def model(h, cfg):
        calls.append("model")
        d = model_dir(h, cfg)
        d.mkdir(parents=True, exist_ok=True)
        for v in fetch.VARS:
            (d / f"{v}.hdf").write_bytes(b"x")
        return 0

    def wl(h, cfg, cor1a=False):
        calls.append("wl")
        d = h / "data" / "wl"
        d.mkdir(parents=True, exist_ok=True)
        for n in ("20131025120000.fts", "20131025_s4c1B.fts", "20131025_d4c2B.fts"):
            (d / n).write_bytes(b"x")
        return 0

    def figure(out, h, half, params, **layers):
        calls.append(f"figure {out} " + ",".join(k for k in ("c2", "cor1", "cor2", "radio") if layers.get(k)))
        p = Path(h) / out
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(b"png")

    def fit(*a, **kw):
        calls.append("fit")
        fits.append(a)
        return fit_results.pop(0) if fit_results else result()

    monkeypatch.setattr(fetch, "model", model)
    monkeypatch.setattr(fetch, "whitelight", wl)
    monkeypatch.setattr(reconstruction, "make_figure", figure)
    monkeypatch.setattr(gcs, "fit", fit)
    monkeypatch.setattr(pipeline, "wl_gaps", lambda *a: {})      # the fake frames above are not FITS files
    # NRH 150.9 MHz from 14:58:00 UT, a burst in images 18 to 27 (15:00:58 to 15:02:27 UT); 228 MHz is missing
    nrh = tmp_path / "data" / "nrh"
    nrh.mkdir(parents=True)
    write_nrh(nrh / NRH_FILE, 150.9, T0, 60, lambda k: (-1.1, -0.2) if 18 <= k <= 27 else None)
    return SimpleNamespace(home=tmp_path, calls=calls, fits=fits, fit_results=fit_results)


def make_config(w):
    cfg = merge(default_config(), new_config("2013-10-25T15:00", "S06E69"))
    (w.home / "config.json").write_text(json.dumps(cfg, indent=2), encoding="utf-8")
    return config(w.home)


def run(w, **kw):
    log = []
    code = pipeline.run(w.home, args(**kw), log.append)
    return code, [x for x in log if x.startswith("[")]


def report(w):
    return (w.home / "event_report.md").read_text(encoding="utf-8")


def fixed_flags(w):
    return [bool(a[5]) for a in w.fits]


def test_steps_run_in_order_and_write_the_report(world):
    code, log = run(world)
    assert code == 0
    assert [x.split("]")[0][1:] for x in log] == ["new", "fetch_model", "fetch_wl", "fetch_ovro_lwa", "fetch_nrh",
                                                  "radio_interval", "gcs", "wl_time", "render_base", "render_wl",
                                                  "render_radio", "report"]
    assert world.calls[:3] == ["model", "wl", "fit"]
    assert world.calls[3].startswith("figure figures/corona3d.png") and world.calls[3].endswith(" ")
    assert world.calls[4].endswith("c2,cor1,cor2") and world.calls[5].endswith("radio")
    assert world.fits[0][4] == ("c2", "cor1b", "cor2b")
    cfg = config(world.home)
    assert cfg["gcs_fit"]["fixed"] == [] and cfg["gcs_wl"]["time"] == "2013-10-25T16:00:00"
    rs = cfg["radio_images"][0]
    assert rs["interval_source"] == "data" and (rs["start"], rs["end"]) == ("2013-10-25T15:00:58", "2013-10-25T15:02:27")


def test_second_run_skips_the_work_that_is_done(world):
    run(world)
    world.calls.clear()
    code, log = run(world)
    assert code == 0 and world.calls == ["model"]                    # fetch.model itself skips present files
    status = dict(x[1:].split("] ", 1) for x in log)
    for step in ("fetch_wl", "gcs", "render_base", "render_wl", "render_radio"):
        assert "skipped" in status[step]
    run(world, force=True)
    assert "fit" in world.calls and any(c.startswith("figure") for c in world.calls)


def test_poor_free_fit_falls_back_to_the_fixed_direction(world):
    world.fit_results += [result(iou=0.3), result(iou=0.8, fixed=True)]
    run(world)
    assert fixed_flags(world) == [False, True]
    cfg = config(world.home)
    assert cfg["gcs_fit"]["fixed"] == ["lon_stonyhurst", "lat_deg"] and cfg["gcs"]["lon_stonyhurst"] == -69.0
    text = report(world)
    assert "free fit: poor (mean IoU 0.30 is below 0.5)." in text and "fixed direction fit: good." in text
    assert "held at the flare site" in text and "cor1b 0.8" in text


def may9(iou, views, fixed=False, dropped=()):
    """A result in the form of the May 9 2024 fits: views kept after the view rule, with their IoU."""
    r = result(fixed=fixed)
    r["t_ref"]["iou"] = dict(iou)
    r["views_used"] = list(views)
    r["view_check"] = dict(enabled=True, dropped=list(dropped), kept="refit" if dropped else "full", refits=[],
                           full=dict(views=["cor1a", "cor2a"], iou=dict(cor1a=0.901, cor2a=0.0), mean_iou=0.45))
    return r


def test_free_fit_with_a_dropped_view_is_good_and_the_fixed_fit_is_not_run(world):
    """May 9 2024: free fit cor1a 0.90 and cor2a 0.0; cor2a is dropped, the free fit is judged on cor1a only."""
    world.fit_results.append(may9(dict(cor1a=0.902), ["cor1a"], dropped=["cor2a"]))
    run(world, views="c2,cor1b,cor2b")
    assert fixed_flags(world) == [False]
    assert config(world.home)["gcs_fit"]["views"] == ["cor1a"]            # no c2: it had no envelope
    text = report(world)
    assert "images of cor1a at" in text and "c2, cor1a" not in text and "views used: cor1a" in text
    assert "The fixed-direction fit was not run, because the free fit was good." in text
    assert "Kept: the free fit." in text and "(1) free fit on all views" in text
    d = json.loads((world.home / "gcs_fit.json").read_text())["decision"]
    assert d["accepted"] == "free" and d["fixed_run"] is False and len(d["order"]) == 5


def test_both_fits_poor_keeps_the_better_one_and_reports_every_step(world):
    world.fit_results += [may9(dict(cor1a=0.45, cor2a=0.20), ["cor1a", "cor2a"]),
                          may9(dict(cor1a=0.487, cor2a=0.268), ["cor1a", "cor2a"], fixed=True)]
    run(world)
    assert fixed_flags(world) == [False, True]
    text = report(world)
    assert "Kept: the fixed direction fit, the better of the two" in text
    assert "free fit: poor" in text and "fixed direction fit: poor" in text
    d = json.loads((world.home / "gcs_fit.json").read_text())["decision"]
    assert d["accepted"] == "fixed direction" and d["fixed_run"] and [x["mode"] for x in d["attempts"]] == ["free", "fixed direction"]
    world.fit_results += [may9(dict(cor1a=0.30), ["cor1a"]), may9(dict(cor1a=0.31), ["cor1a"], fixed=True)]
    run(world, force=True)
    assert "Kept: the free fit, the better of the two" in report(world)       # a tie goes to the free fit


def test_fit_at_a_bound_falls_back_and_a_good_fit_does_not(world):
    world.fit_results += [result(kappa=0.8), result(fixed=True)]       # kappa at its upper bound
    run(world)
    assert fixed_flags(world) == [False, True]
    world.fits.clear()
    world.fit_results.append(result())
    run(world, force=True)
    assert fixed_flags(world) == [False]


def test_failed_free_fit_falls_back(world, monkeypatch):
    def fit(*a, **kw):
        world.fits.append(a)
        if not a[5]:
            raise SystemExit("no CME envelope")
        return result(fixed=True)
    monkeypatch.setattr(gcs, "fit", fit)
    run(world)
    assert fixed_flags(world) == [False, True]
    assert "free fit: failed (no CME envelope)." in report(world)


def test_fit_skip_fixed_and_a_height_given_by_hand(world):
    run(world, fit="skip")
    assert world.fits == [] and config(world.home).get("gcs_fit") is None
    run(world, fit="fixed", height_out=1.77, height_note="kinematic estimate", force=True)
    assert fixed_flags(world) == [True] and world.fits[0][-2:] == (1.77, "kinematic estimate")


def test_when_every_fit_fails_the_shell_is_assumed_and_a_partial_result_exits_zero(world, monkeypatch):
    def fit(*a, **kw):
        raise SystemExit("fewer than two usable fronts")
    monkeypatch.setattr(gcs, "fit", fit)
    code, _ = run(world)
    assert code == 0                                                 # figures exist: a partial result
    text = report(world)
    assert "default one placed above the flare site" in text and "Failed: GCS fit." in text
    assert text.index("## Warnings") < text.index("Failed: GCS fit.") < text.index("## Fit quality")
    world.home.joinpath("event_state.json").unlink()
    assert run(world, strict=True)[0] == 1
    assert (world.home / "figures" / "corona3d.png").exists()


def test_report_lists_measured_assumed_missing_and_figures(world):
    (world.home / "data" / "nrh" / NRH_FILE).unlink()
    run(world)
    text = report(world)
    for h in ("## Measured", "## Assumed", "## Missing data", "## Fit quality", "## Figures"):
        assert h in text
    assert "1000 km/s" in text and "Newkirk" in text and "steady-state" in text
    assert "No nrh files in data/nrh" in text and "rsdb.obs-nancay.fr" in text
    assert "figures/corona3d.png" in text and "figures/wl/corona3d.png" in text and "figures/radio/corona3d.png" not in text
    assert "claude" not in text.lower() and "scene" not in text.lower()


def test_another_event_in_the_same_home_is_refused(world):
    run(world)
    with pytest.raises(SystemExit, match="not 2024-05-09T17:32"):
        pipeline.run(world.home, args(time="2024-05-09T17:32"), lambda s: None)


def test_new_marks_the_radio_intervals_as_default():
    c = new_config("2024-05-09T17:32", "S14W28", earth_stereo=([317.1, -3.2], [329.6, -1.7]), run="r")
    assert [r["interval_source"] for r in c["radio_images"]] == ["default", "default"]


def test_quality_flags_low_iou_and_parameters_at_a_bound():
    assert gcs.quality(result(), FLARE) == []
    assert "mean IoU 0.30" in gcs.quality(result(iou=0.3), FLARE)[0]
    assert any("kappa" in w for w in gcs.quality(result(kappa=0.8), FLARE))
    r = result()
    r["t_ref"]["lon_stonyhurst"] = FLARE[0] + 40                       # the free fit at its longitude bound
    assert any("lon_stonyhurst" in w for w in gcs.quality(r, FLARE))
    assert not any("lon_stonyhurst" in w for w in gcs.quality(r, FLARE, fix_direction=True))
    r["t_ref"]["iou"] = {}
    assert gcs.quality(r, FLARE, True)[0] == "no instrument had an envelope"


def test_burst_run_takes_the_run_with_the_brightest_frame():
    t = lambda s: datetime(2013, 10, 25, 15, 0, 0) + timedelta(seconds=s)
    series = [(t(10 * k), 1e5) for k in range(40)]
    series[3] = (t(30), 2e8)                                           # a single spike
    for k in range(20, 27):
        series[k] = (t(10 * k), 1e9 if k == 23 else 5e8)
    r = radio.burst_run(series, 1e8)
    assert (r["start"], r["end"], r["n"], r["t_peak"]) == (t(200), t(260), 7, t(230))
    assert radio.burst_run(series, 1e10) is None


def test_radio_interval_keeps_a_manual_interval_unless_forced(world):
    cfg = make_config(world)
    notes = radio.select(world.home, cfg, 10)
    rs = config(world.home)["radio_images"][0]
    assert (rs["start"], rs["end"], rs["step_s"]) == ("2013-10-25T15:00:58", "2013-10-25T15:02:27", 30)
    assert any("228 MHz: no image at or above" in m for _, m in notes)
    assert "already chosen from the images" in radio.select(world.home, config(world.home), 10)[0][1]
    user = json.loads((world.home / "config.json").read_text())
    user["radio_images"][0].pop("interval_source")
    user["radio_images"][0].update(start="2013-10-25T14:59:00", end="2013-10-25T15:03:00")
    (world.home / "config.json").write_text(json.dumps(user))
    assert "was given by hand" in radio.select(world.home, config(world.home), 10)[0][1]
    assert config(world.home)["radio_images"][0]["start"] == "2013-10-25T14:59:00"
    radio.select(world.home, config(world.home), 10, force=True)
    assert config(world.home)["radio_images"][0]["start"] == "2013-10-25T15:00:58"


def test_radio_interval_without_a_burst_leaves_the_default(world):
    make_config(world)
    user = json.loads((world.home / "config.json").read_text())
    user["radio_images"][0]["min_tb"] = 1e12
    (world.home / "config.json").write_text(json.dumps(user))
    notes = radio.select(world.home, config(world.home), 10)
    assert any("no frequency reaches min_tb" in m for _, m in notes)
    assert config(world.home)["radio_images"][0]["start"] == "2013-10-25T15:00:00"


def test_nrh_fetch_reports_missing_frequencies_and_the_manual_step(world, capsys):
    cfg = make_config(world)
    cfg["radio_images"][0]["freqs"] = [150.9]
    assert fetch.nrh(world.home, cfg) == 0
    cfg["radio_images"][0]["freqs"] = [150.9, 228.0]
    assert fetch.nrh(world.home, cfg) == 1
    out = capsys.readouterr().out
    assert "NRH 150.9 MHz: files overlap" in out and "NRH 228 MHz: missing" in out and "rsdb.obs-nancay.fr" in out


def test_ovro_lwa_is_fetched_once_for_the_search_window(world, monkeypatch):
    windows = []
    monkeypatch.setattr(fetch, "ovro_lwa", lambda h, cfg: windows.append([(e["start"], e["end"]) for e in cfg["radio_images"]]) or 0)
    a = args(time="2024-05-09T17:32", flare="S14W28", fit="skip")
    pipeline.run(world.home, a, lambda s: None)
    pipeline.run(world.home, a, lambda s: None)
    assert len(windows) == 1 and windows[0][1] == ("2024-05-09T17:22:00", "2024-05-09T17:42:00")


def test_t_ref_also_sets_wl_time_and_is_recorded_in_the_provenance(world):
    run(world, t_ref="2013-10-25T16:05")
    cfg = config(world.home)
    assert cfg["wl_time"] == "2013-10-25T16:05"
    assert cfg["gcs_fit"]["wl_time"] == "2013-10-25T16:05" and cfg["gcs_fit"]["wl_time_from"] == "--t-ref"
    assert world.fits[0][1] == "2013-10-25T16:05"
    run(world, fit="skip", t_ref="2013-10-25T16:10")                 # also without a fit
    assert config(world.home)["wl_time"] == "2013-10-25T16:10"


def test_without_t_ref_the_fit_chooses_it_and_wl_time_follows(world):
    """t_out is the time of the event (the burst); t_ref is not passed and comes from the fit, and the
    report says so. wl_time follows the chosen t_ref; the fit time is not the default time plus 60 min."""
    res = result()
    res["t_ref"]["time"] = "2013-10-25T15:20:00"
    res["t_ref_choice"] = dict(method="auto", rule="earliest developed frame", developed=["c2", "cor1b"])
    world.fit_results.append(res)
    run(world)
    assert world.fits[0][1] is None and world.fits[0][2] == "2013-10-25T15:00:00"      # t_ref auto, t_out = --time
    cfg = config(world.home)
    assert cfg["wl_time"] == "2013-10-25T15:20:00" and cfg["gcs_fit"]["wl_time_from"] == "t_ref chosen automatically"
    assert cfg["gcs_fit"]["t_out"] == "2013-10-25T15:00:00" and cfg["gcs_fit"]["t_ref_choice"]["method"] == "auto"
    text = report(world)
    assert "fit time 2013-10-25T15:20:00 was chosen from the images" in text
    assert "radio and shell time is 2013-10-25T15:00:00" in text


def test_fit_without_a_developed_cme_falls_back_and_the_report_says_so(world):
    res = result()
    res["t_ref_choice"] = dict(method="fallback", rule="r", why="no developed envelope found")
    world.fit_results.append(res)
    run(world)
    assert "No coronagraph frame showed a developed CME" in report(world)


def test_a_failed_figure_is_a_warning_and_exits_zero_unless_strict(world, monkeypatch):
    ok = reconstruction.make_figure

    def figure(out, h, half, params, **layers):
        if "wl" in out:
            raise RuntimeError("boom")
        ok(out, h, half, params, **layers)

    monkeypatch.setattr(reconstruction, "make_figure", figure)
    code, log = run(world)
    assert code == 0
    text = report(world)
    assert "## Warnings" in text and "Failed: figure wl (boom)." in text
    assert "figures/corona3d.png" in text and "figures/radio/corona3d.png" in text
    assert run(world, strict=True, force=True)[0] == 1


def test_no_figure_at_all_is_a_total_failure_with_exit_one(world, monkeypatch):
    def figure(*a, **k):
        raise RuntimeError("boom")
    monkeypatch.setattr(reconstruction, "make_figure", figure)
    code, _ = run(world)
    assert code == 1 and "No figure could be drawn" in report(world)


def test_strict_option_is_on_the_command_line():
    from corona3d import cli
    seen = []
    for extra in ([], ["--strict"]):
        pipeline_run = pipeline.run
        pipeline.run = lambda h, a, log=print: seen.append(a.strict) or 0
        try:
            with pytest.raises(SystemExit):
                cli.main(["--home", "x", "event", "--time", "2024-05-09T17:32", "--flare", "S14W28"] + extra)
        finally:
            pipeline.run = pipeline_run
    assert seen == [False, True]


def test_white_light_figure_is_drawn_without_an_instrument_that_has_a_gap(world, monkeypatch):
    """Regression (2024-05-09): LASCO C2 has no frame between 14:48 and 18:52 UT, and the whole
    white-light figure failed; it is now drawn with COR1 and COR2, and the report names the gap."""
    monkeypatch.setattr(pipeline, "wl_gaps", lambda *a: {"c2": "no frame within 30 min of wl_time 2024-05-09T17:55"})
    code, log = run(world)
    assert code == 0
    assert any(c.startswith("figure figures/wl/corona3d.png") and c.endswith("cor1,cor2") for c in world.calls)
    assert "LASCO C2: no frame within 30 min of wl_time" in report(world)


def test_wl_gaps_names_only_instruments_without_a_frame_near_wl_time(monkeypatch):
    from corona3d import whitelight
    t = datetime(2024, 5, 9, 17, 55)
    frames = {"c2": [(datetime(2024, 5, 9, 14, 48), []), (datetime(2024, 5, 9, 18, 52), [])],
              "cor2": [(datetime(2024, 5, 9, 17, 54), [])]}
    monkeypatch.setattr(whitelight, "_frames", lambda wl, k, sc="B": frames[k])
    out = pipeline.wl_gaps("x", {"wl_time": "2024-05-09T17:55", "time": "2024-05-09T17:32"},
                           {"c2": True, "cor2": True, "cor1": False}, "A")
    assert list(out) == ["c2"] and "18:52:00 (after)" in out["c2"]


def test_wl_gap_when_the_nearest_frame_is_the_base_frame(monkeypatch):
    """Regression (2024-05-09 burst 4): wl_time 18:32 is within 30 min of the first LASCO C2 frame (18:52), but
    that frame is the first of wl_window, the base frame, so the panel failed; the instrument is now a gap."""
    from corona3d import whitelight
    frames = {"c2": [(datetime(2024, 5, 9, 18, 52), []), (datetime(2024, 5, 9, 19, 4), [])],
              "cor2": [(datetime(2024, 5, 9, 17, 7), []), (datetime(2024, 5, 9, 18, 30), [])]}
    monkeypatch.setattr(whitelight, "_frames", lambda wl, k, sc="B": frames[k])
    cfg = {"wl_time": "2024-05-09T18:32", "time": "2024-05-09T17:32",
           "wl_window": ["2024-05-09T16:32", "2024-05-09T20:02"]}
    out = pipeline.wl_gaps("x", cfg, {"c2": True, "cor2": True}, "A")
    assert list(out) == ["c2"] and "base frame" in out["c2"]


def test_white_light_panel_moves_to_the_nearest_frame_when_no_instrument_has_one(world, monkeypatch):
    """No instrument has a usable frame at wl_time: wl_time moves to the nearest drawable frame, the
    coronagraph panel is still drawn, and the report says that the time was moved."""
    from corona3d import whitelight
    monkeypatch.setattr(pipeline, "wl_gaps", ORIGINAL_WL_GAPS)
    t = lambda h, m: datetime(2013, 10, 25, h, m)
    frames = {"c2": [(t(16, 40), []), (t(16, 50), [])], "cor1": [(t(14, 10), []), (t(14, 20), [])],
              "cor2": [(t(14, 30), [])]}
    monkeypatch.setattr(whitelight, "_frames", lambda wl, k, sc="B": frames[k])
    code, log = run(world)
    assert code == 0
    assert config(world.home)["wl_time"] == "2013-10-25T16:50:00"
    assert config(world.home)["gcs_fit"]["wl_time_from"].startswith("nearest frame")
    assert any(c.startswith("figure figures/wl/corona3d.png") and c.endswith("c2") for c in world.calls)
    text = report(world)
    assert "wl_time was moved to 2013-10-25T16:50:00" in text and "## Warnings" in text
    assert "LASCO C2" not in text.split("## Missing data")[1].split("## Warnings")[0]
    run(world)                                                       # a second run keeps the moved time
    assert config(world.home)["wl_time"] == "2013-10-25T16:50:00"


def test_pick_t_ref_takes_the_earliest_time_with_the_most_developed_views():
    from corona3d.gcs import pick_t_ref
    t = lambda m: datetime(2024, 5, 9, 17, m)
    r_out = {"cor1a": 4.0, "cor2a": 15.0}
    meas = {"cor1a": [(t(20), 5.0, 3.0), (t(35), 3.5, 2.8), (t(45), 6.0, 3.2), (t(55), 9.0, 3.99)],
            "cor2a": [(t(40), 4.0, 6.0), (t(52), 7.0, 6.0)]}
    got, dev = pick_t_ref(meas, r_out, t(32))                 # the large frame at 17:20 is before t_out and is ignored
    assert got == t(35) and dev == dict(cor1a=True, cor2a=True)    # COR2-A at 17:40 is within 6 min
    got, dev = pick_t_ref({"cor1a": meas["cor1a"][1:]}, r_out, t(32))
    assert got == t(35) and dev == dict(cor1a=True)
    small = {"cor1a": [(t(35), 1.0, 2.0), (t(55), 9.0, 3.95)]}   # too small, then too close to the edge of the field
    assert pick_t_ref(small, r_out, t(32)) == (None, {})
    late = {"cor2a": [(t(32) + timedelta(minutes=200), 9.0, 6.0)]}  # beyond the search window
    assert pick_t_ref(late, r_out, t(32)) == (None, {})
