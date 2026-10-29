"""The shock diagnostics (speeds, theta_Bn, Mach numbers, self-similar expansion), the heights of the
radio sources in several density models, the speeds of the radio tracks, and the degree of circular
polarisation, against analytic answers on synthetic data."""
from datetime import datetime, timedelta

import numpy as np
import pytest
from astropy.io import fits

from conftest import write_mas
from corona3d import heights, mast_model, nrh, shock

M_P, K_B = 1.67262192e-24, 1.380649e-16


def write_var(d, name, arr_fn):
    """One more MAS variable on the grid of conftest.write_mas: arr_fn(P, T, R) in code units."""
    from pyhdf.SD import SD, SDC
    phi, th, r = np.linspace(0, 2 * np.pi, 33), np.linspace(0, np.pi, 61), np.linspace(1.0, 4.0, 61)
    P, T, R = np.meshgrid(phi, th, r, indexing="ij")
    sd = SD(str(d / f"{name}.hdf"), SDC.WRITE | SDC.CREATE)
    for i, ax in enumerate((phi, th, r)):
        s = sd.create(f"fakeDim{i}", SDC.FLOAT32, (len(ax),))
        s[:] = ax.astype(np.float32)
        s.endaccess()
    s = sd.create("Data-Set-2", SDC.FLOAT32, P.shape)
    s[:] = np.broadcast_to(arr_fn(P, T, R), P.shape).astype(np.float32)
    s.endaccess()
    sd.end()


def test_mas_units_of_temperature_and_flow_follow_from_b_and_n():
    assert mast_model.V_UNIT == pytest.approx(481.37, rel=1e-4)            # km/s
    assert mast_model.T_UNIT == pytest.approx(2.8072e7, rel=1e-3)          # K
    # p0 = rho0 v0^2 = n0 k_B T0, so that p = 2 rho T holds for electrons and protons at T
    assert N0_KT0() == pytest.approx(1e8 * M_P * (mast_model.V_UNIT * 1e5) ** 2, rel=1e-6)


def N0_KT0():
    return 1e8 * K_B * mast_model.T_UNIT


def test_characteristic_speeds():
    assert shock.alfven_speed(1.0, 1e8) == pytest.approx(1.0 / np.sqrt(4 * np.pi * 1e8 * M_P) / 1e5)
    assert shock.sound_speed(1e6) == pytest.approx(np.sqrt(5 / 3 * 2 * K_B * 1e6 / M_P) / 1e5)   # 165.9 km/s
    assert shock.sound_speed(1e6) == pytest.approx(165.9, abs=0.1)
    va, cs = 300.0, 150.0
    assert shock.fast_speed(va, cs, 90.0) == pytest.approx(np.hypot(va, cs))      # perpendicular
    assert shock.fast_speed(va, cs, 0.0) == pytest.approx(va)                     # parallel: the larger
    assert shock.fast_speed(100.0, cs, 0.0) == pytest.approx(cs)
    th = np.linspace(0, 90, 7)
    assert np.all(np.diff(shock.fast_speed(va, cs, th)) >= 0)                     # grows towards perpendicular


def test_theta_bn_is_folded_to_0_90():
    n = np.array([[0, 0, 1.0]] * 4)
    b = np.array([[0, 0, 2.0], [0, 0, -3.0], [1.0, 0, 0], [1.0, 0, 1.0]])
    assert shock.theta_bn(b, n) == pytest.approx([0, 0, 90, 45])


def test_normals_of_a_tube_point_away_from_its_axis():
    a, phi = 0.3, np.linspace(0, 2 * np.pi, 25)                   # closed rings, as PyThea gives them
    zs = np.linspace(1.2, 2.0, 9)
    rings = np.array([[(2 + a * np.cos(f), a * np.sin(f), z) for f in phi] for z in zs])
    n = shock.shell_normals(rings)
    radial = rings - np.array([2.0, 0, 0]) * np.ones_like(rings)
    radial[..., 2] = 0
    radial /= np.linalg.norm(radial, axis=2, keepdims=True)
    assert np.allclose(n, radial, atol=1e-6)


def test_pythea_shell_normals_and_self_similar_speeds():
    pytest.importorskip("PyThea")
    from corona3d.geometry import pythea_shell
    h, v = 2.0, 1000.0
    rings = pythea_shell(h, 0.4, 25.0, 0.0, 0.0, 0.0, 0.0, obstime="2013-10-25T15:00",
                         nbvertsl=10, nbvertcirc=20, nbvertcircshell=40)
    n = shock.shell_normals(rings).reshape(-1, 3)
    x = rings.reshape(-1, 3)
    apex = np.argmax(np.linalg.norm(x, axis=1))
    assert np.linalg.norm(x[apex]) == pytest.approx(h, abs=0.01)
    assert n[apex] @ x[apex] / h == pytest.approx(1.0, abs=0.01)                  # radial at the apex
    vn = shock.normal_speed(x, n, v, h)
    assert vn[apex] == pytest.approx(v, rel=0.01)                                # the apex moves at V_apex
    assert np.nanmax(vn) <= v * 1.01
    axis = np.array([1.0, 0, 0])                                                 # lon 0, lat 0 with Earth at 0
    w = shock.half_width(rings, axis)
    assert shock.half_width(rings * 1.5, axis) == pytest.approx(1.5 * w)         # self-similar


def test_mach_is_empty_where_the_shell_does_not_advance():
    assert np.isnan(shock.mach(-10.0, 100.0)) and shock.mach(500.0, 100.0) == 5.0


def test_track_speed_recovers_a_uniform_motion():
    t0 = datetime(2013, 10, 25, 15)
    times = [t0 + timedelta(seconds=60 * k) for k in range(5)]
    v = 800.0 / shock.R_SUN_KM                                                   # R_sun per s
    xyz = np.array([[1.5 + v * 60 * k, 0.2, -0.1 - v * 60 * k] for k in range(5)])
    s, err, vec = shock.track_speed(times, xyz)
    assert s == pytest.approx(800.0 * np.sqrt(2), rel=1e-6) and err == pytest.approx(0, abs=1e-6)
    assert np.isnan(shock.track_speed(times[:2], xyz[:2])[1])


def test_heights_of_the_density_models():
    assert heights.newkirk_r(173.0) == pytest.approx(1.576, abs=0.001)
    assert heights.newkirk_r(173.0, 6) > heights.newkirk_r(173.0) > heights.newkirk_r(173.0, 1)
    r = heights.saito_r(173.0)
    assert heights.saito_ne(r) == pytest.approx(heights.level(173.0), rel=1e-6)
    assert 1.0 < r < heights.newkirk_r(173.0, 1)
    assert np.isnan(heights.saito_r(900.0))                                       # denser than the model base


def test_mas_height_on_the_line_of_sight_matches_a_spherical_model():
    nk = lambda pts: 4 * 4.2e4 * 10 ** (4.32 / np.linalg.norm(pts, axis=1))
    for x, y in ((-1.1, -0.2), (0.3, 0.4)):
        assert heights.mas_los_r(nk, x, y, 173.0, (357.0, 5.0)) == pytest.approx(heights.newkirk_r(173.0), abs=0.003)
    assert np.isnan(heights.mas_los_r(nk, -1.8, 0.0, 173.0, (357.0, 5.0)))      # farther out than the level


def test_plasma_sampling_of_a_synthetic_run(tmp_path):
    d = write_mas(tmp_path / "run", n_model=1.0, b0=1.0)
    write_var(d, "t002", lambda P, T, R: 1e6 / mast_model.T_UNIT)
    write_var(d, "vr002", lambda P, T, R: 100.0 / mast_model.V_UNIT)
    for v in ("vt002", "vp002"):
        write_var(d, v, lambda P, T, R: 0.0)
    m = mast_model.load(d, density_scale=2.0)
    assert set(mast_model.load_plasma(d, m)) == {"t002", "vr002", "vt002", "vp002"}
    p = dict(model_dlon=0.0, va_density="scaled", density_scale=2.0)
    pl = shock.Plasma(m, p)
    q = np.array([[0, 0, 2.0], [2.0, 0, 0]])                        # pole and equator at r = 2
    s = pl.at(q)
    assert s["ne"] == pytest.approx([2e8, 2e8], rel=1e-4)
    assert s["t"] == pytest.approx([1e6, 1e6], rel=1e-4)
    b0 = mast_model.B_UNIT / 8                                      # dipole, |B| = b0 (2 cos, sin) / r^3
    assert s["b_xyz"][0] == pytest.approx([0, 0, 2 * b0], abs=1e-4 * b0)
    assert s["b_xyz"][1] == pytest.approx([0, 0, -b0], abs=2e-3 * b0)
    assert s["u_xyz"] == pytest.approx(q / 2 * 100.0, rel=1e-3)     # radial flow of 100 km/s
    d_ = shock.diagnostics(q, q / 2, pl, v_apex=600.0, height=2.0)  # radial normals at r = 2: V_n = 600
    assert d_["vn"] == pytest.approx([500.0, 500.0], rel=1e-3)      # minus the flow
    assert d_["theta_bn"] == pytest.approx([0.0, 90.0], abs=0.5)
    assert d_["mach_a"] == pytest.approx(500.0 / d_["va"])
    assert d_["vf"][1] == pytest.approx(np.hypot(d_["va"][1], d_["cs"][1]), rel=1e-3)
    p_rot = dict(p, model_dlon=40.0)                                # a rotated model: the field turns with it
    assert np.linalg.norm(shock.Plasma(m, p_rot).at(q)["b_xyz"], axis=1) == pytest.approx(np.linalg.norm(s["b_xyz"], axis=1), rel=1e-3)


def write_nrh_iv(path, freq, n, v_frac, separate=False):
    yy, xx = np.mgrid[:64, :64]
    i_img = (1e5 + 1e9 * np.exp(-((xx - 15) ** 2 + (yy - 28) ** 2) / 4.5))[None].repeat(n, 0).astype(np.float32)
    v_img = (v_frac * (i_img - 1e5)).astype(np.float32)
    t0 = (15 * 3600) * 1000
    base = [fits.Column("TIME", "1J", array=t0 + 9900 * np.arange(n)), fits.Column("NUMIM", "1J", array=np.arange(n))]
    cols = base + [fits.Column("STOKESI", "4096E", dim="(64,64)", array=i_img)]
    if not separate:
        cols.append(fits.Column("STOKESV", "4096E", dim="(64,64)", array=v_img))
    for p, cc in ((path, cols),) + (((path.with_name(path.name.replace("_i.", "_v.")),
                                      base + [fits.Column("STOKESV", "4096E", dim="(64,64)", array=v_img)]),) if separate else ()):
        t = fits.BinTableHDU.from_columns(cc)
        t.header.update(FREQ=freq, CRPIX1=32, CRPIX2=32, SOLAR_R=16.0, EXPTIME=9.9)
        fits.HDUList([fits.PrimaryHDU(header=fits.Header({"DATE-OBS": "2013-10-25"})), t]).writeto(p)


@pytest.mark.parametrize("separate", [False, True])
def test_degree_of_circular_polarisation(tmp_path, separate):
    write_nrh_iv(tmp_path / "nrh2_1509_h60_20131025_150000c00_i.fts", 150.9, 10, -0.3, separate)
    t = datetime(2013, 10, 25, 15, 0, 20)
    img, _, _, _ = nrh.image_at(tmp_path, 150.9, t)
    v, how = nrh.stokes_v_at(tmp_path, 150.9, t)
    assert v is not None and ("_v." in how) == separate
    src, pk = nrh.polarisation(img, v, 0.5)
    assert pk == pytest.approx(-0.3, abs=1e-3) and src == pytest.approx(-0.3, abs=1e-3)
    rt = dict(instrument=np.array(["nrh", "user"]), freq_mhz=np.array([150.9, 150.9]),
              time_utc=np.array([t.isoformat(), ""]), x_rsun=np.zeros(2))
    p = dict(radio_images=[dict(instrument="nrh", folder=".", freqs=[150.9], level=0.5)])
    out = nrh.polarisation_rows(tmp_path, p, rt)
    assert out["dcp_source"][0] == pytest.approx(-0.3, abs=1e-3) and np.isnan(out["dcp_source"][1])
    assert out["dcp_valid"].tolist() == [True, False]


def test_stokes_i_only_files_give_no_polarisation(tmp_path):
    from test_nrh import write_nrh
    write_nrh(tmp_path / "nrh2_1509_h60_20131025_150000c00_i.fts", 150.9, 54000000, 5)
    with fits.open(tmp_path / "nrh2_1509_h60_20131025_150000c00_i.fts", mode="update") as h:   # drop STOKESV
        h[1] = fits.BinTableHDU.from_columns([c for c in h[1].columns if c.name != "STOKESV"], header=h[1].header)
    v, how = nrh.stokes_v_at(tmp_path, 150.9, datetime(2013, 10, 25, 15, 0, 10))
    assert v is None and "no STOKESV column" in how


def test_analysis_of_a_synthetic_event_writes_tables_and_the_report(synthetic_home):
    from corona3d.paths import config, model_dir
    d = model_dir(synthetic_home, config(synthetic_home))
    write_var(d, "t002", lambda P, T, R: 1.5e6 / mast_model.T_UNIT)
    (synthetic_home / "event_report.md").write_text("# Event report\n\n## Fit quality\n\n- x\n\n## Figures\n\n- None.\n",
                                                    encoding="utf-8")
    with pytest.raises(SystemExit, match="--speed"):
        shock.analyse(synthetic_home, radio=False, trace=False)
    summ, shell, rtab, tracks = shock.analyse(synthetic_home, speed=900.0, speed_note="test", trace=True, radio=False)
    assert summ["apex"]["vn"] == pytest.approx(900.0, rel=0.03)            # radial flow absent in this run
    assert summ["conventions"]["flow"] == "not subtracted"
    adv = shell["vn"] > 0
    assert np.allclose(shell["mach_a"][adv], shell["vn"][adv] / shell["va"][adv])
    assert np.all(np.isnan(shell["mach_a"][~adv]))
    assert set(shell["field"]) <= {"open", "closed", "undetermined"}
    out = shock.write(synthetic_home / "analysis", summ, shell, rtab, tracks)
    assert (out / "shock_shell.csv").exists() and (out / "shock_summary.json").exists()
    shock.update_report(synthetic_home, summ)
    shock.update_report(synthetic_home, summ)                                # replaced, not repeated
    text = (synthetic_home / "event_report.md").read_text(encoding="utf-8")
    assert text.count("## Shock diagnostics") == 1 and text.index("## Shock diagnostics") < text.index("## Figures")


def test_shell_colour_and_height_range_are_off_by_default(synthetic_home):
    from corona3d.paths import config
    from corona3d.reconstruction import reconstruct
    p = config(synthetic_home)
    on = dict(surface=False, open=False, closed=False, streamer=False, gcs=True)
    sc = reconstruct(p, on, synthetic_home)
    assert "gcs_scalars" not in sc and "radio_ranges" not in sc
    write_var(synthetic_home / "data" / p["run"].replace("/", "_"), "t002", lambda P, T, R: 1e6 / mast_model.T_UNIT)
    sc = reconstruct(dict(p, shell_color="mach_f", shock_speed_km_s=800.0), on, synthetic_home)
    label, vals, cmap, lim = sc["gcs_scalars"]
    assert vals.shape == (sc["gcs"].shape[0] * sc["gcs"].shape[1],) and np.nanmax(vals) > 0
    with pytest.raises(ValueError, match="shell_color"):
        reconstruct(dict(p, shell_color="mach_x", shock_speed_km_s=800.0), on, synthetic_home)
