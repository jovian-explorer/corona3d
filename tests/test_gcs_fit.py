"""`corona3d gcs`: the height-time step (the limb view, the selection of the front, and the
extrapolation to the time of the reconstruction) on synthetic numbers with a known answer, and the
keys written into config.json by a mocked fit."""
import json
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pytest

from corona3d import gcs as G
from corona3d.cli import main

FIT_JSON = Path(__file__).resolve().parent.parent / "scripts" / "gcs_fit_2013-10-25.json"
needs_tree = pytest.mark.skipif(not FIT_JSON.exists(), reason="source tree not present")


def test_projection_factor_is_one_at_the_limb_and_zero_along_the_line_of_sight():
    assert G.projection_factor(-53.0, 0.0, [np.cos(np.radians(37)), np.sin(np.radians(37)), 0]) == \
        pytest.approx(1.0, abs=1e-9)                          # 90 deg from the observer
    assert G.projection_factor(0.0, 0.0, [1, 0, 0]) == pytest.approx(0.0, abs=1e-9)


def test_front_profile_and_extrapolation_recover_a_known_rise():
    """Front at 1.9 R_sun at 15:00 rising at 0.09 R_sun/min; noise frames before the eruption
    (small envelopes) and truncated fronts at the edge of a 4 R_sun field must be left out."""
    t0 = datetime(2013, 10, 25, 15, 0)
    rows = [(t0 - timedelta(minutes=10), 2.4, 0.3),          # pre-event noise at the seed
            (t0 - timedelta(minutes=2), 2.1, 0.6)]
    for m in (5, 8, 13, 18):
        rows.append((t0 + timedelta(minutes=m), 1.9 + 0.09 * m, 2.0 + m / 4))
    rows += [(t0 + timedelta(minutes=23), 3.95, 9.0),           # front at the edge: truncated
             (t0 + timedelta(minutes=28), 3.70, 9.5)]           # what remains after it left
    kept = G.front_profile(rows, r_out=4.0)
    assert [r for _, r, _ in kept] == pytest.approx([1.9 + 0.09 * m for m in (5, 8, 13, 18)])
    v, h0, err = G.extrapolate([t for t, _, _ in kept], [r for _, r, _ in kept], t0)
    assert v == pytest.approx(0.09, abs=1e-9) and h0 == pytest.approx(1.9, abs=1e-9)
    assert err == pytest.approx(0.0, abs=1e-6)


def test_front_profile_ignores_small_envelopes_before_the_eruption():
    """Regression (2013-10-25 on real data): blobs of 1.0 to 1.3 R_sun^2 at 2.5 R_sun, 20 and 55
    minutes before the rise, flattened the line from 1039 to 334 km/s. Only the run of frames
    after the last small envelope counts."""
    t0 = datetime(2013, 10, 25, 15, 0)
    rows = [(t0 - timedelta(minutes=15), 2.55, 1.23), (t0 - timedelta(minutes=12), 1.0, 0.0),
            (t0 - timedelta(minutes=5), 2.47, 1.02)]
    rows += [(t0 + timedelta(minutes=m), 1.9 + 0.09 * m, a) for m, a in ((8, 2.2), (13, 3.8), (18, 8.0))]
    kept = G.front_profile(rows, r_out=4.0)
    assert [r for _, r, _ in kept] == pytest.approx([1.9 + 0.09 * m for m in (8, 13, 18)])
    v, h0, _ = G.extrapolate([t for t, _, _ in kept], [r for _, r, _ in kept], t0)
    assert v == pytest.approx(0.09) and h0 == pytest.approx(1.9)


def fake_fit(*args, **kw):
    return json.loads(FIT_JSON.read_text())


@needs_tree
def test_gcs_command_writes_the_shell_and_its_provenance_into_config_json(tmp_path, monkeypatch):
    main(["--home", str(tmp_path), "init"])
    monkeypatch.setattr(G, "fit", fake_fit)
    main(["--home", str(tmp_path), "gcs", "--views", "c2,cor1a,cor1b", "--fix-direction"])
    cfg = json.loads((tmp_path / "config.json").read_text())
    assert cfg["gcs"] == dict(height=1.931, kappa=0.396, alpha_deg=23.46, lon_stonyhurst=-53.262,
                              lat_deg=3.146, tilt_deg=-146.654)
    assert cfg["gcs_wl"]["height"] == 5.538 and cfg["gcs_wl"]["tilt_deg"] == -146.654
    rec = cfg["gcs_fit"]
    assert rec["views"] == ["c2", "cor1a", "cor1b"] and rec["iou"] == dict(c2=0.761, cor1a=0.854, cor1b=0.744)
    assert rec["speed_km_s"] == 1039.0 and rec["height_t_out_err"] == 0.158
    assert rec["fixed"] == ["lon_stonyhurst", "lat_deg"] and rec["flare"] == [-69.0, -6.0]
    assert (tmp_path / "gcs_fit.json").exists()


@needs_tree
def test_the_fit_file_of_2013_10_25_uses_the_renamed_keys():
    d = json.loads(FIT_JSON.read_text())
    assert {"t_ref", "t_out"} <= set(d) and not {"t1524", "t1500"} & set(d)
    assert d["t_out"]["height"] == d["fit"]["height_t_out"] == 1.931


@needs_tree
def test_assumed_height_is_recorded_as_assumed_with_its_note(tmp_path, monkeypatch):
    main(["--home", str(tmp_path), "init"])
    seen = {}

    def fit(*args, **kw):
        seen["args"] = args
        d = json.loads(FIT_JSON.read_text())
        d["t_out"]["height"] = args[-2]
        d["fit"] = dict(height_t_out=args[-2], assumed=True, note=args[-1])
        return d
    monkeypatch.setattr(G, "fit", fit)
    with pytest.raises(SystemExit):                              # a height without a reason
        main(["--home", str(tmp_path), "gcs", "--height-out", "1.77"])
    main(["--home", str(tmp_path), "gcs", "--height-out", "1.77", "--height-note", "kinematic estimate",
          "--running-diff", "10"])
    cfg = json.loads((tmp_path / "config.json").read_text())
    assert cfg["gcs"]["height"] == 1.77
    assert cfg["gcs_fit"]["height_source"] == "assumed" and cfg["gcs_fit"]["height_note"] == "kinematic estimate"
    assert cfg["gcs_fit"]["running_diff_min"] == 10.0 and "speed_km_s" not in cfg["gcs_fit"]
    assert seen["args"][-3] == 10.0


def test_running_difference_takes_the_frame_nearest_to_t_minus_the_interval(monkeypatch):
    t0 = datetime(2024, 5, 9, 17, 0)
    fr = {"cor1a": [(t0 + timedelta(minutes=5 * i), [Path(f"x_s4c1A_{i}")]) for i in range(6)]}
    monkeypatch.setattr(G, "total_b", lambda files: files[0].name)
    monkeypatch.setattr(G, "seq", lambda f: "s4")
    t = t0 + timedelta(minutes=20)
    assert G.base_of(fr, {("cor1a", "s4"): "first"}, "cor1a", t, fr["cor1a"][4][1], None) == "first"
    assert G.base_of(fr, {}, "cor1a", t, fr["cor1a"][4][1], 10) == "x_s4c1A_2"
    assert G.base_of(fr, {}, "cor1a", t0, fr["cor1a"][0][1], 10) is None


@needs_tree
def test_gcs_t_ref_sets_wl_time_and_the_provenance(tmp_path, monkeypatch):
    main(["--home", str(tmp_path), "init"])
    monkeypatch.setattr(G, "fit", fake_fit)
    main(["--home", str(tmp_path), "gcs", "--t-ref", "2013-10-25T15:30"])
    cfg = json.loads((tmp_path / "config.json").read_text())
    assert cfg["wl_time"] == "2013-10-25T15:30"
    assert cfg["gcs_fit"]["wl_time"] == "2013-10-25T15:30" and cfg["gcs_fit"]["wl_time_from"] == "--t-ref"
    main(["--home", str(tmp_path), "gcs"])                          # no --t-ref: wl_time stays as it is
    assert json.loads((tmp_path / "config.json").read_text())["wl_time"] == "2013-10-25T15:30"


# ---- poor views: refit without a view whose IoU is far below the others (synthetic masks) ----

class DiscView:
    """A view whose envelope is a disc of a 100 x 100 px frame; the 'shell' q = (x, r) is a disc, and the IoU is
    the intersection over union of the two masks."""
    yy, xx = np.mgrid[:100, :100]

    def __init__(self, x, r):
        self.mask = self.disc(x, r)

    @classmethod
    def disc(cls, x, r):
        return (cls.xx - x) ** 2 + (cls.yy - 50) ** 2 < r ** 2

    def iou(self, q):
        s = self.disc(*q)
        return float((s & self.mask).sum() / max((s | self.mask).sum(), 1))


def disc_world(*specs):
    """{name: (time, view)} for discs (name, x, r), and the solve and score functions of a fit of (x, r)."""
    from scipy.optimize import differential_evolution
    vs = {n: (None, DiscView(x, r)) for n, x, r in specs}

    def score(views, q):
        return {k: v.iou(q) for k, (_, v) in views.items()}

    def solve(views):
        cost = lambda q: 1.0 - np.mean(list(score(views, q).values()))
        return list(differential_evolution(cost, [(5, 95), (3, 40)], popsize=15, maxiter=60, seed=1, polish=True).x)
    return vs, solve, score


def test_poor_view_is_the_lowest_iou_below_the_limit_and_below_a_fraction_of_the_others():
    assert G.poor_view(dict(cor1a=0.548, cor2a=0.295)) == "cor2a"        # May 9 2024
    assert G.poor_view(dict(c2=0.8, cor1a=0.85, cor1b=0.727)) is None    # a good fit
    assert G.poor_view(dict(a=0.45, b=0.46)) is None                     # both low, none far below the other
    assert G.poor_view(dict(a=0.55, b=0.2)) == "b"
    assert G.poor_view(dict(a=0.05, b=0.2, c=0.9)) == "a"                # the lowest of two poor ones
    assert G.poor_view(dict(a=0.1)) is None                              # a single view is never poor


def test_a_view_that_pulls_the_fit_is_dropped_and_both_fits_are_recorded():
    """May 9 2024 numbers (cor1a 0.548, cor2a 0.295 in the fit to both; cor1a alone 0.88)."""
    vs = {"cor1a": (None, None), "cor2a": (None, None)}
    ious = {("cor1a", "cor2a"): dict(cor1a=0.548, cor2a=0.295), ("cor1a",): dict(cor1a=0.88)}
    solve = lambda views: tuple(views)
    score = lambda views, q: ious[q]
    views, q2, iou2, rec = G.check_views(vs, tuple(vs), ious[tuple(vs)], solve, score)
    assert list(views) == ["cor1a"] and rec["kept"] == "refit" and rec["dropped"] == ["cor2a"] and q2 == ("cor1a",)
    assert rec["full"]["views"] == ["cor1a", "cor2a"] and rec["full"]["iou"]["cor2a"] == 0.295
    r = rec["refits"][0]
    assert r["kept"] and r["views"] == ["cor1a"] and r["mean_iou"] == 0.88 and r["mean_iou_before"] == 0.548
    text = " ".join(G.view_lines(rec))
    assert "Fit to all views (cor1a, cor2a)" in text and "Refit without cor2a" in text and "kept." in text


def test_a_poor_disc_view_is_flagged_from_its_synthetic_mask():
    vs, solve, score = disc_world(("a", 30, 10), ("b", 34, 10), ("c", 45, 10))
    q = solve(vs)
    assert G.poor_view(score(vs, q)) == "c"


def test_a_poor_view_that_costs_the_others_nothing_is_still_dropped():
    """May 9 2024, free fit: cor1a 0.9009 and cor2a 0 in the fit to both, cor1a 0.9017 alone. The gain is 0.001,
    but the view cannot be matched by the shell, so the refit is kept and the mean IoU is 0.90, not 0.45."""
    vs = {"cor1a": (None, None), "cor2a": (None, None)}
    ious = {("cor1a", "cor2a"): dict(cor1a=0.9009, cor2a=0.0), ("cor1a",): dict(cor1a=0.9017)}
    views, q, iou, rec = G.check_views(vs, tuple(vs), ious[tuple(vs)], lambda v: tuple(v), lambda v, q: ious[q])
    assert list(views) == ["cor1a"] and rec["dropped"] == ["cor2a"] and rec["kept"] == "refit"
    assert rec["refits"][0]["mean_iou"] == 0.902 and rec["full"]["mean_iou"] == 0.45
    assert G.quality(dict(t_ref=dict(iou=iou, height=5.4, kappa=0.4, alpha_deg=25.0, lon_stonyhurst=36.0,
                                     lat_deg=3.0, tilt_pythea_deg=-30.0)), (28.0, -14.0)) == []


def test_a_refit_that_makes_the_kept_views_worse_is_not_kept():
    vs = {"a": (None, None), "b": (None, None)}
    ious = {("a", "b"): dict(a=0.60, b=0.10), ("a",): dict(a=0.45)}        # the others lose 0.15 without the view
    views, q, iou, rec = G.check_views(vs, tuple(vs), ious[tuple(vs)], lambda v: tuple(v), lambda v, q: ious[q])
    assert list(views) == ["a", "b"] and rec["kept"] == "full" and not rec["refits"][0]["kept"]
    assert "not kept" in " ".join(G.view_lines(rec))


def test_keep_all_views_fits_nothing_again_and_good_fits_are_left_alone():
    vs, solve, score = disc_world(("a", 30, 10), ("b", 34, 10), ("c", 45, 10))
    q = solve(vs)
    calls = []
    counting = lambda v: calls.append(1) or solve(v)
    _, _, _, rec = G.check_views(vs, q, score(vs, q), counting, score, drop=False)
    assert calls == [] and rec["kept"] == "full" and not rec["enabled"]
    assert "--keep-all-views" in " ".join(G.view_lines(rec))
    good, solve, score = disc_world(("a", 50, 20), ("b", 52, 20))
    q = solve(good)
    _, _, _, rec = G.check_views(good, q, score(good, q), counting, score)
    assert calls == [] and rec["refits"] == [] and "No view was poor" in " ".join(G.view_lines(rec))


def test_drops_stop_at_one_view_and_follow_the_worst_first():
    names = ["a", "b", "c"]
    vs = {n: (None, None) for n in names}
    ious = {("a", "b", "c"): dict(a=0.9, b=0.2, c=0.05), ("a", "b"): dict(a=0.9, b=0.2), ("a",): dict(a=0.95)}
    solve = lambda views: tuple(views)
    score = lambda views, q: ious[q]
    full = ious[("a", "b", "c")]
    views, q, iou, rec = G.check_views(vs, tuple(vs), full, solve, score, improve=0.0)
    assert list(views) == ["a"] and rec["dropped"] == ["c", "b"]            # the worst first; one view stays
    assert [r["dropped"] for r in rec["refits"]] == ["c", "b"]


def test_choose_fit_follows_the_documented_criterion():
    good, poor = "good", "poor"
    f = lambda status, m: dict(status=status, mean_iou=m)
    assert G.choose_fit(f(good, 0.90), None) == "free"                                    # May 9 free fit, cor1a only
    assert G.choose_fit(f(poor, 0.45), f(good, 0.80)) == "fixed direction"                # a good fit beats a poor one
    assert G.choose_fit(f(good, 0.60), f(poor, 0.40)) == "free"
    assert G.choose_fit(f(poor, 0.45), f(poor, 0.377)) == "free"                          # May 9 with all views kept
    assert G.choose_fit(f(poor, 0.30), f(poor, 0.45)) == "fixed direction"                # higher mean IoU wins
    assert G.choose_fit(f(poor, 0.40), f(poor, 0.41)) == "free"                           # a tie goes to the free fit
    assert G.choose_fit(None, f(poor, 0.2)) == "fixed direction" and G.choose_fit(f(poor, 0.2), None) == "free"
    assert G.choose_fit(None, None) is None
    assert G.choose_fit(f(poor, None), f(poor, 0.2)) == "fixed direction"                 # no envelope in the free fit


@needs_tree
def test_the_command_records_dropped_views_in_config_and_gcs_fit_json(tmp_path, monkeypatch):
    main(["--home", str(tmp_path), "init"])
    seen = {}

    def fit(*args, **kw):
        seen.update(kw)
        d = json.loads(FIT_JSON.read_text()) if FIT_JSON.exists() else None
        d["view_check"] = dict(dropped=["cor1a"], kept="refit", enabled=True, refits=[], full={})
        d["t_ref"]["iou"] = dict(c2=0.8, cor1b=0.7)
        return d
    monkeypatch.setattr(G, "fit", fit)
    main(["--home", str(tmp_path), "gcs", "--views", "c2,cor1a,cor1b"])
    assert seen["drop_poor_views"] is True                                  # the default
    rec = json.loads((tmp_path / "config.json").read_text())["gcs_fit"]
    assert rec["views"] == ["c2", "cor1b"] and rec["views_dropped"] == ["cor1a"]
    assert json.loads((tmp_path / "gcs_fit.json").read_text())["view_check"]["dropped"] == ["cor1a"]
    main(["--home", str(tmp_path), "gcs", "--keep-all-views"])
    assert seen["drop_poor_views"] is False


def _frames(start, n=12):
    return {"cor1a": [(start + timedelta(minutes=5 * i), [Path(f"x_s4c1A_{i}")]) for i in range(n)]}


def test_base_frame_is_fixed_by_the_onset_and_not_by_the_start_of_the_download(monkeypatch):
    """The last frame at or before onset minus N minutes is the base, whatever earlier frames were downloaded."""
    monkeypatch.setattr(G, "total_b", lambda files: files[0].name)
    monkeypatch.setattr(G, "seq", lambda f: "s4")
    onset = datetime(2024, 5, 9, 17, 23)
    cut = G.base_cutoff(onset, 10)
    assert cut == datetime(2024, 5, 9, 17, 13)
    for start in (datetime(2024, 5, 9, 16, 0), datetime(2024, 5, 9, 16, 45)):   # a long download and a short one
        fr = _frames(start, 40)
        base, late = G.bases_of(fr, cut)
        assert late == {} and base[("cor1a", "s4")] == [f for f in fr["cor1a"] if f[0] == datetime(2024, 5, 9, 17, 10)][0][1][0].name
    assert G.base_cutoff("2024-05-09T17:23:00", 0) == onset


def test_base_frame_falls_back_to_the_earliest_frame_and_says_so(monkeypatch):
    monkeypatch.setattr(G, "total_b", lambda files: files[0].name)
    monkeypatch.setattr(G, "seq", lambda f: "s4")
    fr = _frames(datetime(2024, 5, 9, 17, 20))
    base, late = G.bases_of(fr, datetime(2024, 5, 9, 17, 10))
    assert base[("cor1a", "s4")] == "x_s4c1A_0" and list(late) == [("cor1a", "s4")]
    base, late = G.bases_of(fr, None)                              # no cutoff: the earliest, nothing reported
    assert base[("cor1a", "s4")] == "x_s4c1A_0" and late == {}


def test_base_is_recorded_in_the_provenance():
    out = {"base": {"onset": "2024-05-09T17:23", "minutes": 10.0, "cutoff": "2024-05-09T17:13:00",
                    "frames": {"cor1a/s4": "2024-05-09T17:10:00"}, "fallback_earliest": []},
           "t_ref": {"height": 3, "kappa": .3, "alpha_deg": 20, "lon_stonyhurst": 28, "lat_deg": -14,
                     "tilt_pythea_deg": 0, "time": "2024-05-09T17:55", "iou": {"cor1a": .8}},
           "t_out": {"height": 1.77, "kappa": .3, "alpha_deg": 20, "lon_stonyhurst": 28, "lat_deg": -14,
                     "tilt_pythea_deg": 0, "time": "2024-05-09T17:32"},
           "fit": {"height_t_out": 1.77, "assumed": True, "note": "x"}}
    rec = G.config_update(out, (28, -14), ("cor1a",))["gcs_fit"]
    assert rec["base"]["cutoff"] == "2024-05-09T17:13:00" and rec["base"]["frames"]["cor1a/s4"].endswith("17:10:00")


def test_auto_t_ref_skips_a_frame_with_another_binning_and_falls_back_when_nothing_develops(monkeypatch):
    """Regression (2024-05-09, COR2-A): three frames of 256 x 256 pixels among 512 x 512 ones stopped the fit with a
    broadcast error; such a frame is skipped. With no developed envelope the fit time is t_out plus 60 min."""
    from types import SimpleNamespace as NS
    t0 = datetime(2024, 5, 9, 17, 32)
    fr = {"cor2a": [(t0 + timedelta(minutes=m), [Path(f"{n}.fts")]) for m, n in ((10, "small"), (25, "big"))]}
    monkeypatch.setattr(G, "total_b", lambda files: NS(data=np.zeros((256, 256) if files[0].name == "small.fts" else (512, 512)),
                                                       meta=None, date=NS(isot="2024-05-09T17:00:00")))
    monkeypatch.setattr(G, "base_of", lambda *a, **k: NS(data=np.zeros((512, 512)), date=NS(isot="2024-05-09T17:00:00")))
    monkeypatch.setattr(G, "envelope", lambda dm, k, flare: (np.zeros((512, 512), bool), None))
    monkeypatch.setattr(G, "sunpy", NS(map=NS(Map=lambda d, m: NS(data=d))), raising=False)
    t, rec = G.auto_t_ref(fr, {}, ("cor2a",), (28.0, -14.0), t0)
    assert t == "2024-05-09T18:32:00" and rec["method"] == "fallback"


def _obs(*lons):
    return {n: [np.cos(np.radians(a)), np.sin(np.radians(a)), 0.0] for n, a in lons}


def test_one_view_or_close_observers_tie_the_direction_to_the_source_region():
    one = G.direction_limit(_obs(("cor1a", 10)))                                # 9 May 2024, COR1-A alone
    assert one["why"] == "single viewpoint" and one["cone_deg"] == 15.0 and one["separation_deg"] == 0.0
    assert one["text"] == "direction constrained by the source region (single viewpoint)"
    close = G.direction_limit(_obs(("c2", 0), ("cor1a", 12), ("cor2a", 12)), cone=0)   # SOHO and STEREO-A 12 deg apart
    assert close["cone_deg"] == 0 and close["separation_deg"] == pytest.approx(12.0) and "within 12 deg" in close["why"]
    assert G.direction_limit(_obs(("c2", 0), ("cor1a", 31))) is None                # at 30 deg and more the fit stays free
    assert G.direction_limit(_obs(("c2", 0), ("cor1a", 150), ("cor1b", -150))) is None   # 2013-10-25, three views


def test_the_cone_narrows_the_direction_bounds_around_the_flare():
    flare = (28.0, -14.0)
    assert G.bounds(flare)[3:5] == [(-12.0, 68.0), (-49.0, 21.0)]                  # free
    assert G.bounds(flare, cone=15)[3:5] == [(13.0, 43.0), (-29.0, 1.0)]
    assert G.bounds(flare, cone=90)[3:5] == G.bounds(flare)[3:5]                    # never wider than the free bounds
    assert len(G.bounds(flare, fix_direction=True)) == 4


def test_a_direction_at_the_edge_of_the_cone_is_not_a_poor_fit_but_is_recorded():
    shell = dict(height=5.5, kappa=0.4, alpha_deg=25.0, lon_stonyhurst=-109.0, lat_deg=-6.0, tilt_pythea_deg=-30.0)
    out = {"t_ref": dict(shell, time="2024-05-09T18:06:00", iou={"cor1a": 0.9}),
           "t_out": dict(shell, height=1.77, time="2024-05-09T17:32:00"),
           "fit": dict(height_t_out=1.77, assumed=True, note="kinematic assumption")}
    assert any("lon_stonyhurst" in w for w in G.quality(out, (-69.0, -6.0)))        # free fit at its bound is poor
    out["direction"] = dict(G.direction_limit(_obs(("cor1a", 10))), constrained=True, flare=[-69.0, -6.0])
    assert G.quality(out, (-69.0, -6.0)) == []                                      # a prior, not a limit of the data
    rec = G.config_update(out, (-69.0, -6.0), ("cor1a",))["gcs_fit"]
    assert rec["direction_constraint"]["why"] == "single viewpoint" and rec["fixed"] == []
    out["direction"]["cone_deg"] = 0.0
    assert G.config_update(out, (-69.0, -6.0), ("cor1a",))["gcs_fit"]["fixed"] == ["lon_stonyhurst", "lat_deg"]


def test_the_commands_pass_the_direction_cone(tmp_path, monkeypatch):
    main(["--home", str(tmp_path), "init"])
    seen = {}
    monkeypatch.setattr(G, "fit", lambda *a, **kw: (seen.update(kw), 1 / 0)[1])
    for args, want in ((["gcs"], 15.0), (["gcs", "--direction-cone", "0"], 0.0), (["gcs", "--direction-cone", "8"], 8.0)):
        with pytest.raises(ZeroDivisionError):
            main(["--home", str(tmp_path), *args])
        assert seen["direction_cone"] == want


def _close_world(ious, hits):
    """Synthetic world for check_views: `ious` {tuple of views: {view: IoU}}, `hits` {tuple of views: [parameters at a bound]}."""
    names = sorted({k for t in ious for k in t})
    vs = {n: (None, None) for n in names}
    return vs, (lambda v: tuple(v)), (lambda v, q: ious[q]), (lambda q, v: hits.get(q, []))


def test_same_spacecraft_pair_the_clean_single_view_beats_the_bound_hit_one_with_the_higher_iou():
    """9 May 2024: cor2a has the higher IoU in the joint fit but its own fit sits at the kappa bound; cor1a alone is clean."""
    ious = {("cor1a", "cor2a"): dict(cor1a=0.28, cor2a=0.62), ("cor2a",): dict(cor2a=0.798), ("cor1a",): dict(cor1a=0.90)}
    vs, solve, score, problems = _close_world(ious, {("cor2a",): ["kappa"]})
    close = lambda v: dict(cone_deg=15.0)                                    # both views from STEREO-A
    views, q, iou, rec = G.check_views(vs, ("cor1a", "cor2a"), ious[("cor1a", "cor2a")], solve, score, limit=close, problems=problems)
    assert list(views) == ["cor1a"] and rec["dropped"] == ["cor2a"] and rec["kept"] == "refit"
    assert len(rec["refits"]) == 2 and {r["dropped"] for r in rec["refits"]} == {"cor1a", "cor2a"}
    assert [r["bounds"] for r in rec["refits"] if r["dropped"] == "cor1a"] == [["kappa"]]
    assert sum(r["chosen"] for r in rec["refits"]) == 1
    text = " ".join(G.view_lines(rec))
    assert "each view was fitted on its own" in text and "cor1a was chosen" in text and "at a bound: kappa" in text


def test_same_spacecraft_pair_the_higher_iou_wins_when_both_or_neither_hit_a_bound():
    ious = {("a", "b"): dict(a=0.1, b=0.6), ("a",): dict(a=0.7), ("b",): dict(b=0.8)}
    vs, solve, score, _ = _close_world(ious, {})
    close = lambda v: dict(cone_deg=15.0)
    for hits in ({}, {("a",): ["kappa"], ("b",): ["alpha_deg"]}):
        views, *_ = G.check_views(vs, ("a", "b"), ious[("a", "b")], solve, score, limit=close, problems=lambda q, v: hits.get(q, []))
        assert list(views) == ["b"]


def test_a_tie_in_iou_keeps_the_view_nearer_the_sun():
    ious = {("cor1a", "cor2a"): dict(cor1a=0.2, cor2a=0.7), ("cor1a",): dict(cor1a=0.80), ("cor2a",): dict(cor2a=0.81)}
    vs, solve, score, problems = _close_world(ious, {})
    views, *_ = G.check_views(vs, ("cor1a", "cor2a"), ious[("cor1a", "cor2a")], solve, score, limit=lambda v: dict(cone_deg=15.0), problems=problems)
    assert list(views) == ["cor1a"]                                           # smaller field of view


def test_separated_views_keep_the_joint_fit_rule_exactly():
    """2013 case: observers far apart, so limit gives None and the lowest-IoU view is dropped as before; one refit only."""
    ious = {("c2", "cor1a"): dict(c2=0.9, cor1a=0.2), ("c2",): dict(c2=0.9)}
    vs, solve, score, problems = _close_world(ious, {("c2",): ["kappa"]})
    views, q, iou, rec = G.check_views(vs, ("c2", "cor1a"), ious[("c2", "cor1a")], solve, score, limit=lambda v: None, problems=problems)
    assert list(views) == ["c2"] and len(rec["refits"]) == 1 and rec["refits"][0]["candidate"] is False
    assert "Refit without cor1a" in " ".join(G.view_lines(rec)) and "each view was fitted" not in " ".join(G.view_lines(rec))


def test_bound_hits_name_the_parameters_and_leave_a_direction_at_the_cone_alone():
    t = dict(height=5.0, kappa=0.1, alpha_deg=59.9, lon_stonyhurst=43.0, lat_deg=-16.0, tilt_pythea_deg=-30.0)
    assert [h[0] for h in G.bound_hits(t, (28.0, -14.0), fix_direction=True)] == ["kappa", "alpha_deg"]
    assert G.quality(dict(t_ref=dict(t, iou=dict(a=0.8))), (28.0, -14.0), True) != []
