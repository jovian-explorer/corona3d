"""Shock physics beyond the shell (corona3d.shockphys): the stand-off relation, the curvature of a
surface, the Rankine-Hugoniot relations and the first critical Mach number, crossings of field lines
with a surface, collapsing traps, threshold crossings in time, and the analysis of a synthetic event
with every option."""
from argparse import Namespace
from datetime import datetime, timedelta

import numpy as np
import pytest

from corona3d import mast_model, shock, shockphys as sp


def test_standoff_relation_of_russell_and_mulligan():
    g = 5 / 3
    assert sp.standoff_ratio(1e4, g) == pytest.approx(0.81 * (g - 1) / (g + 1), rel=1e-6)   # strong shock: 0.2025
    assert sp.standoff_ratio(2.0, g) == pytest.approx(0.81 * ((g - 1) * 4 + 2) / ((g + 1) * 3))
    assert np.isnan(sp.standoff_ratio(1.0, g)) and np.isnan(sp.standoff_ratio(0.5, g))
    m = np.linspace(1.1, 10, 50)
    assert np.all(np.diff(sp.standoff_ratio(m, g)) < 0)                    # the shock moves in as M grows
    assert sp.standoff_ratio(3.0, 4 / 3) < sp.standoff_ratio(3.0, 5 / 3)    # a softer gas compresses more: thinner sheath


def grid_sphere(R, n_th=40, n_ph=48):
    th = np.linspace(0.2, np.pi - 0.2, n_th)
    ph = np.linspace(0, 2 * np.pi, n_ph, endpoint=False)
    T, P = np.meshgrid(th, ph, indexing="ij")
    n = np.stack([np.sin(T) * np.cos(P), np.sin(T) * np.sin(P), np.cos(T)], axis=-1)
    return R * n, n


def test_principal_curvatures_of_a_sphere_and_a_tube():
    x, n = grid_sphere(2.0)
    k1, H = sp.principal_curvatures(x, n)
    assert k1[2:-2] == pytest.approx(0.5, rel=2e-3) and H[2:-2] == pytest.approx(0.5, rel=2e-3)
    a, phi = 0.3, np.linspace(0, 2 * np.pi, 40, endpoint=False)
    z = np.linspace(0, 2, 30)
    tube = np.array([[(2 + a * np.cos(f), a * np.sin(f), zz) for f in phi] for zz in z])
    nt = shock.shell_normals(tube)
    k1, H = sp.principal_curvatures(tube, nt)
    assert k1[1:-1] == pytest.approx(1 / a, rel=5e-3)                     # across the tube
    assert H[1:-1] == pytest.approx(0.5 / a, rel=5e-3)                    # mean of 1/a and 0


def fluxes(d, ma, gamma):
    """Normal fluxes of mass, normal and tangential momentum, and energy on both sides (upstream units)."""
    u2, X, bn = d["u2"], d["X"], d["bn"]
    up = (ma, ma**2 + d["p1"] + d["bt1"] ** 2 / 2, -bn * d["bt1"],
          ma * (ma**2 / 2 + gamma * d["p1"] / (gamma - 1)) + ma * d["bt1"] ** 2)
    dn = (X * u2, X * u2**2 + d["p2"] + d["bt2"] ** 2 / 2, X * u2 * d["v2t"] - bn * d["bt2"],
          ma * ((u2**2 + d["v2t"] ** 2) / 2 + gamma * d["p2"] / ((gamma - 1) * X)) + u2 * d["bt2"] ** 2 - bn * d["v2t"] * d["bt2"])
    return np.array(up), np.array(dn)


@pytest.mark.parametrize("theta,beta,ma", [(90, 0.1, 3.0), (60, 0.5, 2.2), (30, 0.01, 4.0), (5, 1.0, 2.5)])
def test_fast_shock_conserves_mass_momentum_and_energy(theta, beta, ma):
    g = 5 / 3
    d = sp.fast_shock(ma, theta, beta, g)
    assert d is not None and 1 < d["X"] < 4
    up, dn = fluxes(d, ma, g)
    assert dn == pytest.approx(up, rel=1e-8, abs=1e-10)
    assert abs(d["bt2"]) > abs(d["bt1"]) or theta < 1                     # a fast shock bends the field away from the normal


def test_fast_shock_limits():
    g = 5 / 3
    assert sp.fast_shock(300.0, 50.0, 0.1, g)["X"] == pytest.approx(4.0, rel=1e-3)        # strong shock
    ms = 2.0                                                                               # parallel, beta >> 1: gas dynamics
    beta = 20.0
    cs = np.sqrt(g * beta / 2)
    assert sp.fast_shock(ms * cs, 0.0, beta, g)["X"] == pytest.approx((g + 1) * ms**2 / (2 + (g - 1) * ms**2), rel=1e-6)
    assert sp.fast_shock(0.5, 60.0, 0.1, g) is None                                        # below the fast speed


def test_first_critical_mach_number():
    # Edmiston & Kennel (1984): M_c = 2.76 for a perpendicular shock at beta -> 0, 1.53 for a parallel one
    assert sp.critical_mach_exact(90, 1e-3) == pytest.approx(2.76, abs=0.01)
    assert sp.critical_mach_exact(0, 1e-3) == pytest.approx(1.53, abs=0.02)
    th = np.arange(0, 91, 15)
    mc = [sp.critical_mach_exact(t, 0.01) for t in th]
    assert np.all(np.diff(mc) > 0)                                                        # grows towards perpendicular
    b = [1e-3, 0.1, 0.5, 1.0, 3.0]
    mc = [sp.critical_mach_exact(60, x) for x in b]
    assert np.all(np.diff(mc) < 0)                                                        # falls with beta
    assert sp.critical_mach_exact(10, 50.0) == pytest.approx(1.0, abs=0.01)                    # high beta: M_c close to 1
    # at M_c the downstream normal speed equals the downstream sound speed
    m, t, beta = sp.critical_mach_exact(45, 0.2), 45, 0.2
    vf = sp._fast_speed(beta, t, 5 / 3)
    d = sp.fast_shock(m * vf, t, beta)
    assert d["u2"] ** 2 == pytest.approx(5 / 3 * d["p2"] / d["X"], rel=1e-4)
    # the table agrees with the exact solution
    rng = np.random.default_rng(1)
    tt, bb = rng.uniform(1, 89, 12), 10 ** rng.uniform(-2.5, 1.5, 12)
    assert sp.critical_mach(tt, bb) == pytest.approx([sp.critical_mach_exact(a, c) for a, c in zip(tt, bb)], rel=0.01)
    assert np.isnan(sp.critical_mach(np.nan, 0.1)[0])


def test_plasma_beta_and_scaled_plasma():
    assert sp.plasma_beta(100.0, 100.0) == pytest.approx(2 / (5 / 3))

    class P:
        def at(self, xyz):
            k = len(np.atleast_2d(xyz))
            return dict(ne=np.full(k, 1e8), n=np.full(k, 1e8), b=np.ones(k), b_xyz=np.ones((k, 3)), t=np.full(k, 1e6),
                        u_xyz=np.zeros((k, 3)))
    a, b = P().at([[1, 0, 0]]), sp.ScaledPlasma(P(), 4.0).at([[1, 0, 0]])
    assert b["ne"] == pytest.approx(4 * a["ne"]) and b["n"] == pytest.approx(4 * a["n"])
    assert shock.alfven_speed(b["b"], b["n"]) == pytest.approx(shock.alfven_speed(a["b"], a["n"]) / 2)
    assert sp.ScaledPlasma(P(), 4.0, scaled_va=False).at([[1, 0, 0]])["n"] == pytest.approx(a["n"])


class Uniform:
    """A uniform upstream plasma: B along z, n, T, no flow."""

    def __init__(self, b=(0, 0, 1.0), n=1e8, t=1e6):
        self.b, self.n, self.t = np.asarray(b, float), n, t

    def at(self, xyz):
        k = len(np.atleast_2d(xyz))
        return dict(ne=np.full(k, self.n), n=np.full(k, self.n), b_xyz=np.tile(self.b, (k, 1)),
                    b=np.full(k, np.linalg.norm(self.b)), t=np.full(k, self.t), u_xyz=np.zeros((k, 3)))


def sphere_surface(R=2.0, v=1000.0):
    x, n = grid_sphere(R)
    x, n = x.reshape(-1, 3), n.reshape(-1, 3)
    pl = Uniform()
    dg = sp.add_critical(shock.diagnostics(x, n, pl, v, R), 5 / 3)
    return sp.Surface(x, n, dg, np.zeros(len(x)), np.full(len(x), R), R, v, 5 / 3, "shell"), pl


def test_crossings_and_a_collapsing_trap_on_a_sphere():
    S, pl = sphere_surface()
    tol = sp.mesh_tolerance(S)
    f = np.linspace(-0.75 * np.pi, 0.75 * np.pi, 2001)            # a loop from inside the sphere out to 3 R_sun and back
    loop = np.c_[1.5 + 1.5 * np.cos(f), 1.5 * np.sin(f), 0.3 * np.ones_like(f)]
    c = sp.crossings(loop, S, tol)
    assert len(c) == 2 and c[0]["outward"] and not c[1]["outward"]
    assert [np.linalg.norm(x["point"]) for x in c] == pytest.approx([2.0, 2.0], abs=0.01)
    out = np.linalg.norm(loop, axis=1) > 2
    expect = np.sum(np.linalg.norm(np.diff(loop, axis=0), axis=1)[out[:-1] & out[1:]])
    s = dict(shock.DEFAULTS)
    tr, dist = sp.traps([(loop, 0, "closed"), (loop, 0, "open"), (loop + 1e-4, 0, "closed")], S, pl, s,
                        radio_pts=np.array([[3.0, 0.0, 0.3], [0.0, 0.0, 2.5]]))
    assert len(tr) == 1 and tr[0]["n_lines"] == 2                  # the open line is no trap; the copy is merged
    assert tr[0]["length_rsun"] == pytest.approx(expect, rel=0.01)
    assert tr[0]["top_r"] == pytest.approx(np.hypot(3.0, 0.3), abs=1e-3)
    assert tr[0]["theta_bn_a"] == pytest.approx(90 - np.degrees(np.arcsin(0.3 / 2)), abs=1.5)   # B along z
    assert dist[0] == pytest.approx(0.0, abs=1e-3) and dist[1] > 1 and tr[0]["n_radio_near"] == 1
    inside = np.c_[np.linspace(0, 1.5, 50), np.zeros(50), np.zeros(50)]
    assert sp.crossings(inside, S, tol) == []


def test_stand_off_surface_lies_ahead_of_a_sphere():
    x, n = grid_sphere(2.0)
    keep = np.ones(x.shape[:2], bool)
    s = dict(shock.DEFAULTS)
    pl = Uniform(b=(0, 0, 1.0), n=1e8, t=1e6)
    S0 = sp.surface_state(x, n, keep, 1.0, pl, 1500.0, 2.0, s, standoff=False)
    S = sp.surface_state(x, n, keep, 1.0, pl, 1500.0, 2.0, s, standoff=True)
    assert np.allclose(S0.xyz, x[keep]) and np.all(S0.delta == 0)
    adv = S0.dg["mach_f"] > 1
    assert np.all(S.delta[adv] > 0) and np.all(S.delta[~adv] == 0)
    expect = 2.0 * np.minimum(sp.standoff_ratio(S.dg["vn"] * 0 + S0.dg["mach_f"]), 2.0)  # uniform plasma: M at the shock = M at the shell
    assert S.delta[adv] == pytest.approx(expect[adv], rel=1e-3)
    assert S.rc == pytest.approx(2.0, rel=5e-3)
    r = np.linalg.norm(S.xyz, axis=1)
    assert np.all(r[adv] > 2.0)
    assert S.dg["vn"][adv] == pytest.approx(1500.0 * r[adv] / 2.0, rel=1e-6)          # the shock moves self-similarly
    S2 = sp.surface_state(x, n, keep, 1.5, pl, 1500.0, 2.0, s, standoff=True)
    assert S2.rc == pytest.approx(3.0, rel=5e-3) and S2.ht == 3.0


def test_threshold_crossings_and_quantiles():
    t0 = datetime(2013, 10, 25, 15)
    ts = [t0 + timedelta(minutes=i) for i in range(5)]
    c = sp.threshold_crossings(ts, [0.5, 0.8, 1.2, 1.4, 0.9])
    assert [x["direction"] for x in c] == ["up", "down"]
    assert c[0]["time"] == t0 + timedelta(minutes=1.5)
    q = sp.quantiles(np.arange(101.0))
    assert (q["median"], q["p16"], q["p84"]) == (50.0, 16.0, 84.0)
    assert sp.quantiles([np.nan])["n"] == 0


def test_shock_physics_of_a_synthetic_event(synthetic_home):
    pytest.importorskip("PyThea")
    from corona3d.paths import config, model_dir
    from test_shock import write_var
    d = model_dir(synthetic_home, config(synthetic_home))
    write_var(d, "t002", lambda P, T, R: 1.5e6 / mast_model.T_UNIT)
    p = dict(radio_sources=[dict(freq=150.9, x=[-1.3, -1.25], y=[0.1, 0.12])],
             shock=dict(standoff=True, mc_samples=3, evolution_step_s=120.0))
    tables = {}
    summ, shell, rtab, tracks = shock.analyse(synthetic_home, p=p, speed=900.0, speed_note="test", tables=tables)
    adv = shell["vn"] > 0
    assert np.allclose(shell["mach_f_mc"][adv], shell["mach_f"][adv] / shell["mach_c"][adv], equal_nan=True)
    assert set(tables) == {"shock_surface", "shock_connectivity", "shock_traps", "shock_evolution", "shock_mc", "shock_mc_radio"}
    so = summ["standoff"]
    assert so["apex_delta_rsun"] > 0 and so["apex"]["r"] == pytest.approx(summ["apex"]["r"] + so["apex_delta_rsun"], abs=5e-3)
    assert len(tables["shock_connectivity"]) == 2 * len(rtab["x"]) and "behind_shock_front" in rtab
    assert summ["evolution"]["n_times"] == len(tables["shock_evolution"]) > 3
    assert summ["uncertainty"]["n"] == 3 and len(tables["shock_mc"]) == 3 and len(tables["shock_mc_radio"]) == len(rtab["x"])
    out = shock.write(synthetic_home / "analysis_phys", summ, shell, rtab, tracks, tables)
    for f in ("shock_surface", "shock_connectivity", "shock_evolution", "shock_mc", "shock_mc_radio"):
        assert (out / f"{f}.csv").exists()
    text = "\n".join(shock.report_lines(summ))
    for key in ("Stand-off shock surface", "Magnetic connectivity", "Time evolution", "Ensemble of 3 members", "M_c"):
        assert key in text
    a = Namespace(speed=900.0, speed_note="test", no_trace=True, no_radio=True, out="analysis_cli", standoff=True,
                  standoff_radius="mean", mc=None, mc_seed=None, no_connectivity=False, no_traps=False, no_evolution=True,
                  evolution_step=None)
    shock.run(synthetic_home, config(synthetic_home), a)
    assert (synthetic_home / "analysis_cli" / "shock_surface.csv").exists()
    assert not (synthetic_home / "analysis_cli" / "shock_evolution.csv").exists()
