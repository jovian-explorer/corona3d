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

    def fit(*a):
        calls.append("fit")
        fits.append(a)
        return fit_results.pop(0) if fit_results else result()

    monkeypatch.setattr(fetch, "model", model)
    monkeypatch.setattr(fetch, "whitelight", wl)
    monkeypatch.setattr(reconstruction, "make_figure", figure)
    monkeypatch.setattr(gcs, "fit", fit)
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
                                                  "radio_interval", "gcs", "render_base", "render_wl",
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


def test_fit_at_a_bound_falls_back_and_a_good_fit_does_not(world):
    world.fit_results += [result(kappa=0.8), result(fixed=True)]       # kappa at its upper bound
    run(world)
    assert fixed_flags(world) == [False, True]
    world.fits.clear()
    world.fit_results.append(result())
    run(world, force=True)
    assert fixed_flags(world) == [False]


def test_failed_free_fit_falls_back(world, monkeypatch):
    def fit(*a):
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


def test_when_every_fit_fails_the_shell_is_assumed_and_the_exit_code_is_one(world, monkeypatch):
    def fit(*a):
        raise SystemExit("fewer than two usable fronts")
    monkeypatch.setattr(gcs, "fit", fit)
    code, _ = run(world)
    assert code == 1
    text = report(world)
    assert "default one placed above the flare site" in text and "Failed: GCS fit." in text
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


def test_without_t_ref_wl_time_is_left_alone(world):
    run(world)
    assert config(world.home)["wl_time"] == new_config("2013-10-25T15:00", "S06E69")["wl_time"]
    assert "wl_time_from" not in config(world.home)["gcs_fit"]
