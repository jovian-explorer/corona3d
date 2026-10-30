"""Density from the polarised brightness, placement of the radio sources (model surface, F or H), the scattering
correction, the OVRO-LWA position check, and the lanes of dynamic spectra, against analytic answers."""
import json
from datetime import datetime, timedelta

import numpy as np
import pytest

from corona3d import heights, ovro_lwa, pb, scattering, spectra


# ------------------------------------------------------------------ Thomson scattering
def test_van_de_hulst_coefficients_tend_to_point_source():
    r = 200.0
    s2 = 1.0 / r**2
    A, B, C, D = pb.vdh_coefficients(r)
    assert A == pytest.approx(s2, rel=1e-3) and C == pytest.approx(s2, rel=1e-3)
    assert B == pytest.approx(2 / 3 * s2, rel=1e-3) and D == pytest.approx(2 / 3 * s2, rel=1e-3)
    A, B, C, D = pb.vdh_coefficients(1.0)          # at the photosphere: finite
    assert np.isfinite([A, B, C, D]).all() and C == pytest.approx(4 / 3, rel=1e-4)


def test_thomson_far_from_the_sun_is_point_source_scattering():
    n0, rho = 1e5, 20.0
    pB, tB = pb.thomson(lambda pts: np.full(len(pts), n0), [rho], [0.0], (0.0, 0.0), u=0.6, step=0.05, r_out=600.0)
    k = pb.K_THOMSON * n0
    a = np.arctan(np.sqrt(600.0**2 - rho**2) / rho)                    # the line of sight ends at r_out
    i4 = (a + np.sin(a) * np.cos(a)) / rho**3                           # Int dz / r^4
    i2 = 2 * a / rho                                                    # Int dz / r^2
    assert pB[0] == pytest.approx(k * rho**2 * i4, rel=2e-3)            # Int rho^2 / r^4 dz
    assert tB[0] == pytest.approx(k * (2 * i2 - rho**2 * i4), rel=2e-3)  # Int (2 - sin^2) / r^2 dz
    pB2, _ = pb.thomson(lambda pts: np.full(len(pts), 2 * n0), [rho], [0.0], (0.0, 0.0), u=0.6, step=0.05, r_out=600.0)
    assert pB2[0] == pytest.approx(2 * pB[0])                            # linear in the density


@pytest.mark.parametrize("k", [0.5, 1.0])
def test_stokes_from_three_polarisers(k):
    I, Q, U = 10.0, 3.0, -2.0
    ang = [0.0, 120.0, 240.0] if k == 0.5 else [60.0, 0.0, -60.0]
    ims = [k * (I + Q * np.cos(np.radians(2 * a)) + U * np.sin(np.radians(2 * a))) for a in ang]
    tb, pbv = pb.stokes_pb(ims, ang, k)
    assert tb == pytest.approx(I) and pbv == pytest.approx(np.hypot(Q, U))


def test_secchi_onboard_factors_bias_and_calibration():
    h = {"IP_00_19": " 41 76  3 50  3 50  1 97", "DIV2CORR": True, "BIASMEAN": 700.0, "IPSUM": 3, "OFFSETCR": 0}
    assert pb.sebip_factor(h) == 16.0                         # two divisions by 4; the division by 2 already corrected
    assert pb.sebip_factor(dict(h, DIV2CORR=False)) == 32.0
    assert pb.secchi_bias(h) == 700.0 * 16
    assert pb.secchi_bias(dict(h, OFFSETCR=5.0)) == 0.0
    c = pb.cor1_calfac({"OBSRVTRY": "STEREO_B", "DATE-OBS": "2008-01-17T02:20:15", "IPSUM": 3})
    assert c == pytest.approx(7.080e-11 / 16, rel=1e-6)


def test_fit_scale_and_resolve():
    pa = np.repeat(np.arange(5.0, 360.0, 10.0), 4)
    r = np.tile([1.7, 2.0, 2.5, 3.0], 36)
    f = pb.fit_scale(np.full(len(pa), 1.8), pa, r, [[1.6, 2.2], [2.2, 3.1]])
    assert f["scale"] == pytest.approx(1.8) and f["err"] == pytest.approx(0.0, abs=1e-12)
    assert [q["scale"] for q in f["ranges"]] == pytest.approx([1.8, 1.8])
    assert pb.resolve({"density_scale": 2.0})["density_scale"] == 2.0
    assert pb.resolve({"density_scale": "pb", "pb_fit": {"scale": 1.3}})["density_scale"] == 1.3
    with pytest.raises(ValueError):
        pb.resolve({"density_scale": "pb"})


def test_polar_bins_of_a_radial_image():
    y, x = np.mgrid[-4:4:200j, -4:4:200j]
    obs = dict(x=x, y=y, pb=1e-8 / np.hypot(x, y) ** 3)
    b = pb.polar_bins(obs, [1.6, 3.6], 0.5, 30.0)
    assert len(b["r"]) > 30
    assert np.allclose(b["pb"], 1e-8 / b["r"] ** 3, rtol=0.1)


# ------------------------------------------------------------------ scattering
def newkirk_r(f, fold=1, harmonic=False):
    return 4.32 / np.log10((f * 1e6 / (2 if harmonic else 1) / 8980) ** 2 / (fold * 4.2e4))


def test_scattering_shift_reproduces_chrysaphi_2018():
    ne = scattering.newkirk_profile(1.0)
    d32 = scattering.shift(32, newkirk_r(32), ne, harmonic=False)
    d40 = scattering.shift(40, newkirk_r(40), ne, harmonic=False)
    assert 0.5 < d32[0] < 0.75 and d32[1] < d32[0] < d32[2]       # about 0.6 R_sun at 32 MHz
    sep = (newkirk_r(32) + d32[0]) - (newkirk_r(40) + d40[0])
    assert sep == pytest.approx(0.3, abs=0.06)                     # apparent separation of 32 and 40 MHz
    dh = scattering.shift(32, newkirk_r(32, harmonic=True), ne, harmonic=True)
    assert dh[0] < 0.5 * d32[0]                                     # harmonic waves scatter much less


def test_scattering_correction_and_size():
    x, y = scattering.correct(1.2, -0.5, 1.5, 0.3)
    assert np.hypot(x, y) == pytest.approx(np.hypot(1.2, -0.5) * 1.5 / 1.8)
    assert scattering.fwhm_deg(32) == pytest.approx(11.78 * 32**-0.98)


# ------------------------------------------------------------------ placement of the sources
P = {"earth": [0.0, 0.0], "radio_emission": "H", "radio_deproject": "sphere"}


def test_emission_choice():
    assert heights.emission_of(P) == "H"
    assert heights.emission_of(P, {"emission": "f"}) == "F"
    e = {"emission": {"150.9": "F", "default": "H"}}
    assert heights.emission_of(P, e, 150.9) == "F" and heights.emission_of(P, e, 432.0) == "H"
    assert heights.emission_of(P, {"emission": "auto"}) == "auto"
    with pytest.raises(ValueError):
        heights.emission_of(P, {"emission": "X"})
    assert heights.emission_from_polarisation(0.45, True, 0.3)[0] == "F"
    assert heights.emission_from_polarisation(0.05, True, 0.3)[0] == "H"
    assert heights.emission_from_polarisation(0.9, False, 0.3)[0] == "H"
    assert not heights.custom(dict(P, scattering={"apply": False}), {"instrument": "nrh"})
    assert heights.custom(dict(P, radio_deproject="model"))
    assert heights.custom(P, {"emission": "F"})
    assert heights.custom(P, {"position_offset": [0.01, 0.0]})


def test_locate_on_sphere_and_on_model_surface():
    from corona3d.geometry import deproject
    q = heights.locate(-0.9, 0.2, 150.9, P, True, "sphere")
    assert q["r"] == pytest.approx(heights.newkirk_r(150.9))
    assert q["xyz"] == pytest.approx(deproject([-0.9], [0.2], q["r"], 0.0, 0.0)[0])
    n0 = 1e10
    ne_at = lambda pts: n0 * np.linalg.norm(pts, axis=1) ** -8          # plasma level of f at r = (n0 / n)^(1/8)
    for harm in (True, False):
        q = heights.locate(-0.9, 0.2, 150.9, P, harm, "model", ne_at)
        r_true = (n0 / heights.level(150.9, harm)) ** (1 / 8)
        assert q["surface"] == "model" and q["r"] == pytest.approx(r_true, abs=0.01)
        assert np.linalg.norm(q["xyz"]) == pytest.approx(q["r"])
        assert q["xyz"][2] > 0                                           # observer side
    q = heights.locate(-3.0, 0.0, 150.9, P, True, "model", ne_at)      # line of sight misses the level: sphere
    assert q["surface"] == "sphere"


def test_locate_with_scattering_moves_the_source_inward():
    q0 = heights.locate(-1.0, 0.3, 40.0, P, False, "sphere")
    q1 = heights.locate(-1.0, 0.3, 40.0, P, False, "sphere", scatter=True)
    assert q1["dr"] > 0 and np.hypot(q1["x"], q1["y"]) < np.hypot(-1.0, 0.3)
    assert q1["r"] == pytest.approx(q0["r"])                               # same sphere, corrected position on it
    assert q1["dr_lo"] <= q1["dr"] <= q1["dr_hi"]


# ------------------------------------------------------------------ OVRO-LWA position check
def disc_image(dx, dy, radius=1.3, tb=1e6, n=128, extent=6.0):
    x = np.linspace(-extent, extent, n)
    X, Y = np.meshgrid(x, x)
    img = np.where(np.hypot(X - dx, Y - dy) < radius, tb, 0.0) + 1e3 * np.random.default_rng(1).standard_normal(X.shape)
    return img, x, x


def test_quiet_disc_centroid_and_criteria():
    img, x, y = disc_image(0.08, -0.05)
    cx, cy, rad, tb, con = ovro_lwa.disc(img, x, y)
    assert cx == pytest.approx(0.08, abs=0.03) and cy == pytest.approx(-0.05, abs=0.03)
    assert rad == pytest.approx(1.3, abs=0.1) and con < 2
    X, Y = np.meshgrid(x, y)
    storm = img + np.where(np.hypot(X - 1.0, Y + 1.0) < 0.3, 1e9, 0.0)   # a noise storm: not a quiet disc
    cx, cy, rad, tb, con = ovro_lwa.disc(storm, x, y)
    assert con > ovro_lwa.QUIET["contrast"]


def test_refraction_fit_and_offset_given(tmp_path):
    res = {f: dict(dx=0.01 + 30.0 / f**2, dy=-0.02, n=3) for f in (34.0, 52.0, 70.0, 84.0)}
    fit = ovro_lwa.refraction_fit(res)
    assert fit["dx"] == pytest.approx((0.01, 30.0), rel=1e-6)
    assert ovro_lwa.offset_for(tmp_path, {"position_offset": [0.1, -0.2]}, 50.0) == (0.1, -0.2)
    assert ovro_lwa.offset_for(tmp_path, {}, 50.0) == (0.0, 0.0)
    f = tmp_path / ovro_lwa.OFFSET_FILE
    f.parent.mkdir(parents=True)
    f.write_text(json.dumps({"bands": {"52.4": dict(dx=0.05, dy=0.01, n=4), "84.6": dict(dx=np.nan, dy=np.nan, n=0)}}))
    assert ovro_lwa.offset_for(tmp_path, {"position_offset": "quiet_sun"}, 52.4) == (0.05, 0.01)
    assert ovro_lwa.offset_for(tmp_path, {"position_offset": "quiet_sun"}, 84.6) == (0.0, 0.0)


# ------------------------------------------------------------------ dynamic spectra
def test_band_split_mach_number():
    X = 1.2**2
    assert spectra.mach_bandsplit(1.2) == pytest.approx(np.sqrt(X * (X + 5) / (2 * (4 - X))))
    assert np.isnan(spectra.mach_bandsplit(2.1))
    assert spectra.mach_bandsplit(1.0) == pytest.approx(1.0)              # no compression: M_A = 1


def test_shock_speed_in_newkirk_density():
    ne = scattering.newkirk_profile(4.0)
    D = -2e-3
    v, r = spectra.shock_speed(80.0, D, ne, harmonic=True)
    assert r == pytest.approx(heights.newkirk_r(80.0, 4.0, True), abs=1e-3)
    assert v == pytest.approx(2 * D * r**2 / (-4.32 * np.log(10)) * 695700.0, rel=1e-3)


def synthetic_spectrum(drift=-2e-3, f0=80.0, split=1.25, harmonic=True, seed=0):
    t0 = datetime(2024, 5, 9, 17, 30)
    t = [t0 + timedelta(seconds=k * 0.5) for k in range(1200)]
    f = np.geomspace(20, 300, 300)[::-1]
    rng = np.random.default_rng(seed)
    d = 10 * rng.standard_normal((len(f), len(t))) + 100
    ts = np.arange(len(t)) * 0.5
    for k, tt in enumerate(ts):
        if 120 <= tt <= 480:
            fc = f0 * np.exp(drift * (tt - 300))
            for ratio, amp in [(1.0, 60)] + ([(split, 40)] if split else []) + ([(2.0, 30)] if harmonic else []):
                d[:, k] += amp * np.exp(-0.5 * (np.log(f / (ratio * fc)) / 0.03) ** 2)
    return dict(time=t, freq=f, data=d, units="digits"), t0, t[-1]


def test_lanes_found_in_a_synthetic_spectrum():
    sp, t0, t1 = synthetic_spectrum()
    s = spectra.settings({})
    res = spectra.analyse_one("test", sp, t0, t1, s, {"newkirk4": scattering.newkirk_profile(4.0)}, "H")
    kinds = {l["kind"]: l for l in res["lanes"]}
    assert "main" in kinds
    assert kinds["main"]["drift_per_s"] == pytest.approx(-2e-3, rel=0.15)
    assert "split_upper" in kinds or "split_lower" in kinds
    assert res["band_split"]["ratio"] == pytest.approx(1.25, rel=0.05)
    assert "harmonic" in kinds and kinds["main"]["emission"] == "F"       # a lane at twice its frequency: F
    assert kinds["main"]["speed_km_s"]["newkirk4"] > 0


def test_given_lane_is_refined():
    sp, t0, t1 = synthetic_spectrum(split=None, harmonic=False)
    s = dict(spectra.settings({}), lanes=[{"t0": "2024-05-09T17:32:30", "f0": 80 * np.exp(-2e-3 * (150 - 300)) * 1.03,
                                          "t1": "2024-05-09T17:37:30", "f1": 80 * np.exp(-2e-3 * (450 - 300)) * 0.97}])
    res = spectra.analyse_one("test", sp, t0, t1, s, {"newkirk4": scattering.newkirk_profile(4.0)}, "H")
    assert res["lanes"][0]["given"] and res["lanes"][0]["drift_per_s"] == pytest.approx(-2e-3, rel=0.1)


def test_callisto_listing_and_choice():
    html = ('<a href="BIR_20131025_145000_01.fit.gz">x</a> <a href="BIR_20131025_150000_01.fit.gz">x</a>'
            '<a href="GLASGOW_20131025_150000_59.fit.gz">x</a> <a href="XYZ_20131025_120000_01.fit.gz">x</a>')

    class R:
        def __init__(self, t):
            self.t = t.encode()

        def read(self):
            return self.t

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False
    lst = spectra.callisto_list(datetime(2013, 10, 25), opener=lambda url, timeout: R(html))
    assert len(lst) == 4
    ch = spectra.callisto_choose(lst, datetime(2013, 10, 25, 14, 58), datetime(2013, 10, 25, 15, 10))
    assert ch[0][0] == "BIR_01" and len(ch[0][1]) == 2
    assert all(not n.startswith("XYZ") for _, fs in ch for n in fs)


def test_read_callisto(tmp_path):
    from astropy.io import fits
    nf, nt = 5, 8
    prim = fits.PrimaryHDU(np.arange(nf * nt, dtype=np.uint8).reshape(nf, nt))
    prim.header["DATE-OBS"], prim.header["TIME-OBS"] = "2013/10/25", "15:00:00.000"
    cols = fits.ColDefs([fits.Column("TIME", f"{nt}D", array=[np.arange(nt) * 0.25]),
                         fits.Column("FREQUENCY", f"{nf}D", array=[[90.0, 80.0, 80.0, 60.0, 50.0]])])
    p = tmp_path / "BIR_20131025_150000_01.fit.gz"
    fits.HDUList([prim, fits.BinTableHDU.from_columns(cols)]).writeto(p)
    sp = spectra.read_callisto([p])
    assert list(sp["freq"]) == [50.0, 60.0, 80.0, 90.0] and sp["data"].shape == (4, nt)
    assert sp["time"][1] - sp["time"][0] == timedelta(seconds=0.25)


def test_report_section(tmp_path):
    from corona3d.report_sections import update
    assert update(tmp_path, "X", ["a"]) is None
    (tmp_path / "event_report.md").write_text("# R\n\n## Measured\n\nm\n\n## Figures\n\nf\n", encoding="utf-8")
    update(tmp_path, "Dynamic spectra", ["- one"])
    update(tmp_path, "Dynamic spectra", ["- two"])
    text = (tmp_path / "event_report.md").read_text(encoding="utf-8")
    assert text.count("## Dynamic spectra") == 1 and "- two" in text and "- one" not in text
    assert text.index("## Dynamic spectra") < text.index("## Figures")


@pytest.mark.parametrize("cmd", [["pb", "--help"], ["sources", "--help"], ["spectra", "--help"]])
def test_cli_commands_exist(cmd, capsys):
    from corona3d.cli import main
    with pytest.raises(SystemExit) as e:
        main(cmd)
    assert e.value.code == 0 and "usage" in capsys.readouterr().out


def test_render_options_reach_the_parameters(synthetic_home, monkeypatch):
    import corona3d.reconstruction as rc
    from corona3d.cli import main
    seen = {}
    monkeypatch.setattr(rc, "make_figure", lambda out, h, half, params, **k: seen.update(params))
    main(["--home", str(synthetic_home), "render", "--radio-surface", "model", "--emission", "F", "--scattering"])
    assert seen["radio_deproject"] == "model" and seen["radio_emission"] == "F" and seen["scattering"] == {"apply": True}
