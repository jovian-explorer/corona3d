"""Shock diagnostics of the CME shell in the MAS/MAST corona (`corona3d shock`).

The steady MAS/MAST run holds no CME, so that the plasma it gives at a point of the GCS shell is the
undisturbed (upstream) corona into which the shell advances. At every point of the shell above
1.02 R_sun the module samples the electron density, the magnetic field, the temperature, and the
flow of the run, and derives

    v_A = |B| / sqrt(4 pi n m_p)                        Alfven speed
    c_s = sqrt(gamma 2 k_B T / m_p)                     sound speed (p = 2 n k_B T, gamma = 5/3)
    v_f = sqrt{(v_A^2 + c_s^2 + sqrt[(v_A^2 + c_s^2)^2 - 4 v_A^2 c_s^2 cos^2 theta_Bn]) / 2}
                                                        fast-magnetosonic speed along the normal
    theta_Bn = arccos(|B . n| / |B|)                    angle between B and the shell normal, 0 to 90 deg
    V_n = (V_apex / h) (x . n) - u . n                  normal speed of the self-similar shell
                                                        relative to the upstream flow u
    M_A = V_n / v_A,  M_f = V_n / v_f,  M_s = V_n / c_s

The GCS shell with fixed aspect ratio and half-angle is self-similar about Sun centre (every length
scales with the apex height h), so that a point x of the shell moves with the velocity
(dh/dt / h) x. The apex speed dh/dt is the speed of the height-time fit of `corona3d gcs`.
Points whose normal speed is not positive (the side of the shell that faces the Sun) do not
advance into the corona, and their Mach numbers are left empty.

At the radio sources, the same quantities are sampled at the deprojected source position, and the
shell point nearest to the source (on the shell scaled to the time of the source) gives the normal
and the normal speed. The heights of the sources are also computed in several density models
(`heights`), and the degree of circular polarisation is measured where NRH Stokes V is available.
"""
import csv
import json
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial.transform import Rotation

from . import mast_model
from .mast_model import K_B, M_P

GAMMA = 5.0 / 3.0                         # ratio of specific heats of the MAST energy equation
R_SUN_KM = 695700.0
DEFAULTS = dict(
    mesh=dict(nbvertsl=10, nbvertcirc=20, nbvertcircshell=40),   # analysis shell (finer than the drawn one)
    r_min=1.02,                           # R_sun: shell points below this are left out (model boundary)
    gamma=GAMMA,
    flow=True,                            # subtract the flow of the run when vr002, vt002, vp002 are present
    temperature_k=None,                   # a constant temperature [K] when t002 is not used (then assumed)
    nose_deg=30.0,                        # the nose: normal within this angle of the CME axis
    flank_deg=60.0,                       # the flanks: normal at least this angle from the CME axis
    quasi_perp_deg=45.0,                  # theta_Bn above this is quasi-perpendicular
    height_models=["newkirk1", "newkirk2", "newkirk4", "newkirk6", "saito", "mas"],
    # shock physics beyond the shell (corona3d.shockphys, METHODS Section 16)
    standoff=False,                       # also place a shock surface ahead of the shell at the stand-off distance
    standoff_radius="min",                # radius of curvature: min (1 / largest principal curvature) or mean
    standoff_gamma=None,                  # gamma of the stand-off relation (None: gamma)
    standoff_max_rc=2.0,                  # the stand-off is at most this many radii of curvature
    standoff_iterations=8,
    connectivity=True,                    # trace the field line of every radio source to the shock surface
    traps=True,                           # closed loops crossed twice by the shock (collapsing traps)
    trap_near_rsun=0.1,                   # a radio source this close to a trap lies on it
    trap_merge_rsun=0.05,                 # traps whose crossings lie this close are one trap
    evolution=True,                       # diagnostics at a series of times through the burst
    evolution_step_s=60.0,
    evolution_pad_min=10.0,               # minutes before the first and after the last radio source
    evolution_min_height=1.1,             # R_sun: earlier shells are left out
    mc_samples=0,                         # ensemble size of the uncertainty (0: off)
    mc_seed=0,
    mc_height_err_rsun=None,              # 1-sigma apex height (None: the error of the fit, else 0.15, assumed)
    mc_speed_frac=0.1,                    # 1-sigma of the apex speed as a fraction of it
    mc_density_scale=[1.5, 3.0],          # density_scale drawn uniformly in this range
    mc_runs="auto",                       # other model runs: "auto" (folders in data/), a list of run names, or []
)


# ------------------------------------------------------------------ pure physics (NumPy only)
def alfven_speed(b_gauss, n_cm3):
    """v_A [km/s] for |B| [G] and the number density n [cm^-3] (one proton mass per particle)."""
    return np.asarray(b_gauss, float) / np.sqrt(4 * np.pi * np.asarray(n_cm3, float) * M_P) / 1e5


def sound_speed(t_k, gamma=GAMMA):
    """c_s [km/s] = sqrt(gamma p / rho) with p = 2 n k_B T and rho = n m_p (electrons and protons at
    the temperature T [K])."""
    return np.sqrt(gamma * 2 * K_B * np.asarray(t_k, float) / M_P) / 1e5


def fast_speed(va, cs, theta_deg):
    """Phase speed of the fast-magnetosonic wave [km/s] at the angle theta [deg] between the wave
    normal and B: v_f^2 = (v_A^2 + c_s^2 + sqrt((v_A^2 + c_s^2)^2 - 4 v_A^2 c_s^2 cos^2 theta)) / 2."""
    va, cs = np.asarray(va, float), np.asarray(cs, float)
    c2 = np.cos(np.radians(theta_deg)) ** 2
    s = va**2 + cs**2
    return np.sqrt(0.5 * (s + np.sqrt(np.clip(s**2 - 4 * va**2 * cs**2 * c2, 0, None))))


def theta_bn(b, n):
    """Angle [deg] between the field b (k, 3) and the normal n (k, 3), folded to 0..90 deg."""
    b, n = np.atleast_2d(b), np.atleast_2d(n)
    c = np.abs(np.sum(b * n, axis=1)) / np.maximum(np.linalg.norm(b, axis=1) * np.linalg.norm(n, axis=1), 1e-30)
    return np.degrees(np.arccos(np.clip(c, 0, 1)))


def spherical_to_cartesian(vr, vt, vp, xyz):
    """Cartesian components of a vector given by its (r, theta, phi) components at the points xyz."""
    x, y, z = np.atleast_2d(xyz).T
    r = np.linalg.norm(np.atleast_2d(xyz), axis=1)
    th = np.arccos(np.clip(z / r, -1, 1))
    ph = np.arctan2(y, x)
    st, ct, sp, cp = np.sin(th), np.cos(th), np.sin(ph), np.cos(ph)
    return np.c_[vr * st * cp + vt * ct * cp - vp * sp, vr * st * sp + vt * ct * sp + vp * cp, vr * ct - vt * st]


def shell_normals(rings):
    """Outward unit normals (n_axis, n_circle, 3) of a GCS shell given as rings around its axis. The
    rings of PyThea are closed (the last point repeats the first), and the normal points away from the
    centre of its ring."""
    rings = np.asarray(rings, float)
    m = rings.shape[1] - 1 if np.allclose(rings[:, 0], rings[:, -1]) else rings.shape[1]
    x = rings[:, :m]
    t_circ = np.roll(x, -1, axis=1) - np.roll(x, 1, axis=1)
    t_axis = np.gradient(x, axis=0)
    n = np.cross(t_axis, t_circ)
    centre = x.mean(axis=1, keepdims=True)
    n *= np.sign(np.sum(n * (x - centre), axis=2, keepdims=True) + 1e-300)
    n /= np.maximum(np.linalg.norm(n, axis=2, keepdims=True), 1e-30)
    if m < rings.shape[1]:
        n = np.concatenate([n, n[:, :1]], axis=1)
    return n


def normal_speed(xyz, normals, v_apex, height):
    """Normal speed [km/s] of a self-similar shell: (v_apex / height) (x . n)."""
    return v_apex / height * np.sum(np.atleast_2d(xyz) * np.atleast_2d(normals), axis=1)


def mach(vn, speed):
    """vn / speed where the shell advances (vn > 0) and the speed is positive, otherwise NaN."""
    vn, speed = np.asarray(vn, float), np.asarray(speed, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where((vn > 0) & (speed > 0), vn / speed, np.nan)


def half_width(rings, axis_dir):
    """Largest distance of the shell from its axis (the line from Sun centre along axis_dir) [R_sun]."""
    p = np.asarray(rings, float).reshape(-1, 3)
    a = np.asarray(axis_dir, float) / np.linalg.norm(axis_dir)
    return float(np.linalg.norm(p - np.outer(p @ a, a), axis=1).max())


def track_speed(times, xyz):
    """Speed [km/s] of a straight-line fit to positions xyz (k, 3) [R_sun] against time, with its
    standard error (NaN for fewer than three points); returns (speed, error, velocity vector)."""
    t = np.array([(x - times[0]).total_seconds() for x in times], float)
    xyz = np.atleast_2d(np.asarray(xyz, float))
    if len(t) < 2 or np.ptp(t) == 0:
        return np.nan, np.nan, np.full(xyz.shape[1], np.nan)
    if len(t) < 3:
        v = (xyz[-1] - xyz[0]) / (t[-1] - t[0]) * R_SUN_KM
        return float(np.linalg.norm(v)), np.nan, v
    v, var = np.empty(xyz.shape[1]), np.empty(xyz.shape[1])
    for i in range(xyz.shape[1]):
        (v[i], _), cov = np.polyfit(t, xyz[:, i], 1, cov=True)
        var[i] = cov[0, 0]
    v *= R_SUN_KM
    s = float(np.linalg.norm(v))
    err = float(np.sqrt(np.sum((v / max(s, 1e-30)) ** 2 * var)) * R_SUN_KM) if s > 0 else np.nan
    return s, err, v


# ------------------------------------------------------------------ sampling of the model
class Plasma:
    """Interpolates the run at points of the reconstruction (Carrington Cartesian, R_sun), with the
    rotation model_dlon of the configuration. Quantities: n_e (scaled density), n (density of the
    Alfven speed, after va_density), B (Cartesian, G), T (K), u (flow, Cartesian, km/s)."""

    def __init__(self, m, p, temperature_k=None, flow=True):
        self.r_min, self.r_max = m["r"][0], m["r"][-1]
        axes = (m["phi"], m["theta"], m["r"])
        self.R = Rotation.from_euler("z", p["model_dlon"], degrees=True).as_matrix()
        n_mass = m["ne"] if p["va_density"] == "scaled" else m["ne"] / p["density_scale"]
        self.f = {k: RegularGridInterpolator(axes, a, bounds_error=False, fill_value=None)
                  for k, a in (("ne", m["ne"]), ("n", n_mass), ("br", m["br"]), ("bt", m["bt"]), ("bp", m["bp"]))}
        self.has_t = "T" in m and temperature_k is None
        if self.has_t:
            self.f["T"] = RegularGridInterpolator(axes, m["T"], bounds_error=False, fill_value=None)
        self.t_const = temperature_k
        self.has_flow = flow and all(k in m for k in ("vr", "vt", "vp"))
        if self.has_flow:
            for k in ("vr", "vt", "vp"):
                self.f[k] = RegularGridInterpolator(axes, m[k], bounds_error=False, fill_value=None)

    def _q(self, xyz_model):
        x, y, z = xyz_model.T
        r = np.linalg.norm(xyz_model, axis=1)
        return np.c_[np.mod(np.arctan2(y, x), 2 * np.pi), np.arccos(np.clip(z / r, -1, 1)),
                     np.clip(r, self.r_min, self.r_max)]

    def at(self, xyz):
        """Dictionary of arrays at the points xyz (k, 3) of the reconstruction."""
        xyz = np.atleast_2d(np.asarray(xyz, float))
        xm = xyz @ self.R                                  # into the frame of the model
        q = self._q(xm)
        out = {k: self.f[k](q) for k in ("ne", "n")}
        b = spherical_to_cartesian(self.f["br"](q), self.f["bt"](q), self.f["bp"](q), xm) @ self.R.T
        out["b_xyz"], out["b"] = b, np.linalg.norm(b, axis=1)
        if self.has_t:
            out["t"] = self.f["T"](q)
        else:
            out["t"] = np.full(len(xyz), np.nan if self.t_const is None else float(self.t_const))
        if self.has_flow:
            out["u_xyz"] = spherical_to_cartesian(self.f["vr"](q), self.f["vt"](q), self.f["vp"](q), xm) @ self.R.T
        else:
            out["u_xyz"] = np.zeros_like(xyz)
        return out


def diagnostics(xyz, normals, plasma, v_apex, height, gamma=GAMMA, b_at=None):
    """Upstream speeds, theta_Bn, normal speed, and Mach numbers at the points xyz with outward
    normals; b_at: plasma samples at other points (the radio sources) to use instead of those at xyz.
    Returns a dictionary of arrays (speeds in km/s, angles in deg)."""
    s = b_at or plasma.at(xyz)
    va = alfven_speed(s["b"], s["n"])
    cs = sound_speed(s["t"], gamma)
    th = theta_bn(s["b_xyz"], normals)
    vf = fast_speed(va, np.nan_to_num(cs, nan=0.0), th)
    vf = np.where(np.isfinite(cs), vf, np.nan)
    u_n = np.sum(s["u_xyz"] * np.atleast_2d(normals), axis=1)
    vn = normal_speed(xyz, normals, v_apex, height) - u_n
    return dict(ne=s["ne"], b=s["b"], t=s["t"], va=va, cs=cs, vf=vf, theta_bn=th, u_n=u_n, vn=vn,
                mach_a=mach(vn, va), mach_f=mach(vn, vf), mach_s=mach(vn, cs))


# ------------------------------------------------------------------ the analysis of a working directory
def settings(p):
    """DEFAULTS updated with the key `shock` of the configuration."""
    s = dict(DEFAULTS, **(p.get("shock") or {}))
    s["mesh"] = dict(DEFAULTS["mesh"], **((p.get("shock") or {}).get("mesh") or {}))
    return s


def apex_speed(home, p, speed=None, note=None):
    """(speed [km/s], source text, height-time frames) of the apex: --speed when given (assumed), the
    speed of the gcs_fit record of the configuration, or that of gcs_fit.json in the working directory."""
    if speed is not None:
        return float(speed), f"given with --speed ({note or 'no reason given'}); assumed", []
    g = p.get("gcs_fit") or {}
    frames = []
    f = Path(home) / "gcs_fit.json"
    rec = json.loads(f.read_text(encoding="utf-8")) if f.exists() else {}
    frames = (rec.get("height_time") or {}).get("frames", [])
    if g.get("speed_km_s"):
        return float(g["speed_km_s"]), f"height-time fit of the front in {g.get('front')} (config.json gcs_fit)", frames
    if (rec.get("fit") or {}).get("speed_km_s"):
        h = rec.get("height_time") or {}
        return float(rec["fit"]["speed_km_s"]), f"height-time fit of the front in {h.get('instrument')} (gcs_fit.json)", frames
    raise SystemExit("no apex speed: fit the shell with `corona3d gcs` (its height-time fit gives the speed), or give "
                     "--speed KM_S with --speed-note")


def axis_direction(g, p):
    """Unit vector of the CME axis (Carrington Cartesian) for a shell of the configuration."""
    from .geometry import unit_vec
    return unit_vec(g["lon_stonyhurst"] + p["earth"][0], g["lat_deg"])


def field_class(m, p, xyz):
    """open, closed, or undetermined for the field line through each point (the tracer of corona3d)."""
    from .geometry import FieldTracer
    R = Rotation.from_euler("z", p["model_dlon"], degrees=True).as_matrix()
    return np.array([k for _, k in FieldTracer(m).field_lines(np.atleast_2d(xyz) @ R)])


def _summ(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return dict(n=int(len(x)), median=float(np.median(x)), min=float(x.min()), max=float(x.max())) if len(x) else dict(n=0)


def analyse(home, p=None, speed=None, speed_note=None, trace=True, radio=True, log=print, tables=None):
    """Shock diagnostics, kinematics, radio heights, and polarisation of a working directory, with the
    shock physics of corona3d.shockphys (stand-off surface, critical Mach number, connectivity, traps,
    time evolution, and ensemble) as set in the key `shock` of the configuration.
    Returns (summary dict, shell table dict, radio table dict, track list); the further tables
    (shock_surface, shock_connectivity, shock_traps, shock_evolution, shock_mc, shock_mc_radio) are put
    into the dictionary `tables` when one is given."""
    from . import heights, nrh
    from . import shockphys as sp
    from .geometry import pythea_shell
    from .paths import config, merge, model_dir
    from .whitelight import utc
    from .pb import resolve
    p = resolve(merge(config(home), p))                    # density_scale "pb": the scale fitted by corona3d pb
    s = settings(p)
    tables = {} if tables is None else tables
    g = p["gcs"]
    v_apex, v_src, frames = apex_speed(home, p, speed, speed_note)
    d = model_dir(home, p)
    m = mast_model.load(d, density_scale=p["density_scale"], va_density=p["va_density"])
    found = mast_model.load_plasma(d, m)
    if "t002" not in found and s["temperature_k"] is None:
        raise SystemExit(f"t002.hdf (temperature) missing in {d}: run `corona3d fetch --model-extra`, or set "
                         "shock.temperature_k in config.json (then assumed)")
    plasma = Plasma(m, p, s["temperature_k"], s["flow"])
    rings = pythea_shell(g["height"], g["kappa"], g["alpha_deg"], g["lon_stonyhurst"], g["lat_deg"], g["tilt_deg"] + 90.0,
                         p["earth"][0], obstime=p["time"], **s["mesh"])
    normals = shell_normals(rings)
    if np.allclose(rings[:, 0], rings[:, -1]):                         # drop the repeated point of each ring
        rings, normals = rings[:, :-1], normals[:, :-1]
    keep2 = np.linalg.norm(rings, axis=2) > s["r_min"]
    xyz, nrm = rings[keep2], normals[keep2]
    h = float(g["height"])
    dg = sp.add_critical(diagnostics(xyz, nrm, plasma, v_apex, h, s["gamma"]), s["gamma"])
    axis = axis_direction(g, p)
    ang = np.degrees(np.arccos(np.clip(nrm @ axis, -1, 1)))
    region = np.where(ang <= s["nose_deg"], "nose", np.where(ang >= s["flank_deg"], "flank", "shoulder"))
    lines = sp.trace_lines(m, p, xyz) if trace else None
    kind = np.array([k for _, _, k in lines]) if trace else np.full(len(xyz), "")
    from .api import carrington_angles
    lon, lat, r = carrington_angles(xyz)
    shell = dict(x=xyz[:, 0], y=xyz[:, 1], z=xyz[:, 2], r=r, lon_carrington_deg=lon, lat_deg=lat,
                 nx=nrm[:, 0], ny=nrm[:, 1], nz=nrm[:, 2], normal_axis_deg=ang, region=region, field=kind, **dg)
    apex = int(np.argmax(r))
    w = half_width(rings, axis)
    adv = dg["vn"] > 0
    qperp = dg["theta_bn"] > s["quasi_perp_deg"]

    def stats(sel, dg=dg, adv=adv, qperp=qperp):
        sup = np.nan_to_num(dg["mach_f_mc"]) > 1
        a = sel & adv
        return {k: _summ(dg[k][sel]) for k in ("vn", "va", "cs", "vf", "theta_bn", "mach_a", "mach_f", "mach_c", "mach_f_mc")} | dict(
            quasi_perpendicular_fraction=float(np.mean(qperp[a])) if a.any() else None,
            supercritical_fraction=float(np.mean(sup[a])) if a.any() else None,
            qperp_supercritical_fraction=float(np.mean((sup & qperp)[a])) if a.any() else None)

    summ = dict(
        time=p["time"], height_rsun=h, apex_speed_km_s=v_apex, apex_speed_source=v_src,
        conventions=dict(va_density=p["va_density"], density_scale=p["density_scale"], gamma=s["gamma"],
                         temperature="MAS t002" if plasma.has_t else f"constant {s['temperature_k']} K (assumed)",
                         flow="MAS vr002, vt002, vp002 subtracted" if plasma.has_flow else "not subtracted",
                         expansion="self-similar GCS shell about Sun centre: V(x) = (V_apex / h) x",
                         nose_deg=s["nose_deg"], flank_deg=s["flank_deg"], quasi_perp_deg=s["quasi_perp_deg"],
                         critical_mach="first critical Mach number from the Rankine-Hugoniot relations (Edmiston & Kennel 1984)"),
        kinematics=dict(half_width_rsun=w, lateral_speed_km_s=v_apex * w / h,
                        flank_normal_speed_km_s=_summ(dg["vn"][region == "flank"]),
                        shells=[dict(time=f["time"], height_rsun=f["height"], half_width_rsun=f["height"] * w / h) for f in frames]),
        apex={k: (float(dg[k][apex]) if np.ndim(dg[k]) == 1 else None)
              for k in ("vn", "ne", "b", "t", "va", "cs", "vf", "theta_bn", "mach_a", "mach_f", "mach_s", "beta", "mach_c", "mach_f_mc")}
        | dict(r=float(r[apex]), lon_carrington_deg=float(lon[apex]), lat_deg=float(lat[apex]), field=str(kind[apex])),
        nose=stats(region == "nose"), flank=stats(region == "flank"), all=stats(np.ones(len(r), bool)),
        n_points=int(len(r)), n_advancing=int(adv.sum()))
    if trace:
        summ["closed_field"] = stats(kind == "closed")
        summ["open_field"] = stats(kind == "open")
    log(f"shock: apex speed {v_apex:g} km/s ({v_src}); height {h:g} R_sun at {p['time']}; {adv.sum()} of {len(r)} shell points advance")

    # the shock surface: the shell itself, or the surface at the stand-off distance ahead of it
    standoff = bool(s["standoff"])
    surf_at = {}

    def surface(k):
        k = round(float(k), 9)
        if k not in surf_at:
            surf_at[k] = (sp.Surface(xyz * k, nrm, dict(dg), np.zeros(len(xyz)), np.full(len(xyz), np.nan), h * k, v_apex,
                                     s["gamma"], "shell") if (k == 1.0 and not standoff) else
                          sp.surface_state(rings, normals, keep2, k, plasma, v_apex, h, s, standoff))
        return surf_at[k]
    S0 = surface(1.0)
    s_lines = lines
    if standoff:
        s_lines = sp.trace_lines(m, p, S0.xyz) if trace else None
        s_kind = np.array([k for _, _, k in s_lines]) if trace else np.full(len(S0.xyz), "")
        sang = np.degrees(np.arccos(np.clip(S0.nrm @ axis, -1, 1)))
        sreg = np.where(sang <= s["nose_deg"], "nose", np.where(sang >= s["flank_deg"], "flank", "shoulder"))
        slon, slat, sr = carrington_angles(S0.xyz)
        tables["shock_surface"] = dict(x=S0.xyz[:, 0], y=S0.xyz[:, 1], z=S0.xyz[:, 2], r=sr, lon_carrington_deg=slon, lat_deg=slat,
                                       nx=S0.nrm[:, 0], ny=S0.nrm[:, 1], nz=S0.nrm[:, 2], normal_axis_deg=sang, region=sreg,
                                       field=s_kind, delta_rsun=S0.delta, rc_rsun=S0.rc, driver_mach_f=dg["mach_f"], **S0.dg)
        sadv, sq = S0.dg["vn"] > 0, S0.dg["theta_bn"] > s["quasi_perp_deg"]
        sa = apex                                           # the point in front of the apex of the shell

        def st(sel):
            return stats(sel, S0.dg, sadv, sq)
        capped = (S0.delta > 0) & np.isclose(S0.delta, s["standoff_max_rc"] * np.nan_to_num(S0.rc))
        summ["standoff"] = dict(
            relation="Delta/R_c = 0.81 [(gamma-1) M^2 + 2] / [(gamma+1) (M^2 - 1)] (Russell & Mulligan 2002; Farris & Russell 1994)",
            radius=s["standoff_radius"], gamma=s.get("standoff_gamma") or s["gamma"], max_rc=s["standoff_max_rc"],
            residual_rc=S0.residual, apex_delta_rsun=float(S0.delta[apex]), apex_rc_rsun=float(S0.rc[apex]),
            apex_driver_mach_f=float(dg["mach_f"][apex]),
            delta=_summ(np.where(S0.delta > 0, S0.delta, np.nan)), n_shifted=int((S0.delta > 0).sum()), n_capped=int(capped.sum()),
            apex={k: float(S0.dg[k][sa]) for k in ("vn", "va", "vf", "theta_bn", "mach_a", "mach_f", "mach_c", "mach_f_mc")}
            | dict(r=float(sr[sa]), field=str(s_kind[sa])),
            nose=st(sreg == "nose"), flank=st(sreg == "flank"), all=st(np.ones(len(sr), bool)))
        if trace:
            summ["standoff"]["closed_field"] = st(s_kind == "closed")
            summ["standoff"]["open_field"] = st(s_kind == "open")
        log(f"shock: stand-off surface {summ['standoff']['apex_delta_rsun']:.3f} R_sun ahead of the apex "
            f"(R_c {summ['standoff']['apex_rc_rsun']:.2f} R_sun, driver M_f {summ['standoff']['apex_driver_mach_f']:.2f})")
    summ["surface"] = "standoff" if standoff else "shell"

    # radio sources
    rtab, tracks = {}, []
    pts, times, scale = np.empty((0, 3)), [], np.empty(0)
    t0 = utc(p["time"])
    if radio:
        from .api import build
        rec = build(home, p, surface=False, open=False, closed=False, streamer=False, gcs=False, radio=True)
        rt = rec.radio
        if rt and len(rt["x"]):
            pts = np.c_[rt["x"], rt["y"], rt["z"]]
            times = [utc(t) if t else t0 for t in rt["time_utc"]]
            scale = np.array([(h + v_apex * (t - t0).total_seconds() / R_SUN_KM) / h for t in times])
            near_i, near_d = [], []
            for q, k in zip(pts, scale):                      # shell scaled to the time of the source
                dd = np.linalg.norm(xyz * k - q, axis=1)
                near_i.append(int(np.argmin(dd)))
                near_d.append(float(dd.min()))
            near_i = np.array(near_i)
            # the normal speed of a self-similar shell does not change with its scale, so that the shell at
            # t_out gives it; the upstream plasma is that at the source itself
            rd = sp.add_critical(diagnostics(xyz[near_i], nrm[near_i], plasma, v_apex, h, s["gamma"], b_at=plasma.at(pts)),
                                 s["gamma"])
            inside = np.array([float(nrm[i] @ (q - xyz[i] * k)) < 0 for i, q, k in zip(near_i, pts, scale)])
            hr = heights.ranges(rt, p, m, s["height_models"])
            pol = nrh.polarisation_rows(home, p, rt)
            rtab = {c: rt[c] for c in ("instrument", "freq_mhz", "time_utc", "x_rsun", "y_rsun", "r_rsun", "x", "y", "z",
                                       "lon_carrington_deg", "lat_deg")}
            rtab.update(hr)
            rtab.update(shell_distance_rsun=np.array(near_d), behind_shell_front=inside,
                        shell_field=kind[near_i] if trace else np.full(len(pts), ""),
                        **{f"shell_{k}": v for k, v in rd.items() if k in ("vn", "theta_bn")},
                        **{k: v for k, v in rd.items() if k not in ("vn", "theta_bn", "u_n")})
            rtab.update(pol)
            for key in sorted({(i, f) for i, f in zip(rt["instrument"], rt["freq_mhz"])}):
                sel = (rt["instrument"] == key[0]) & (rt["freq_mhz"] == key[1]) & (rt["time_utc"] != "")
                if sel.sum() < 2:
                    continue
                tt = [utc(t) for t in rt["time_utc"][sel]]
                v3, e3, _ = track_speed(tt, pts[sel])
                vs, es, _ = track_speed(tt, np.c_[rt["x_rsun"][sel], rt["y_rsun"][sel]])
                tracks.append(dict(instrument=key[0], freq_mhz=float(key[1]), n=int(sel.sum()), start=str(rt["time_utc"][sel][0]),
                                   end=str(rt["time_utc"][sel][-1]), speed_3d_km_s=v3, speed_3d_err_km_s=e3,
                                   speed_sky_km_s=vs, speed_sky_err_km_s=es))
            sup = np.nan_to_num(rd["mach_f_mc"]) > 1
            summ["radio"] = dict(n=int(len(pts)), mach_a=_summ(rd["mach_a"]), mach_f=_summ(rd["mach_f"]),
                                 va=_summ(rd["va"]), vf=_summ(rd["vf"]), theta_bn=_summ(rd["theta_bn"]),
                                 mach_c=_summ(rd["mach_c"]), mach_f_mc=_summ(rd["mach_f_mc"]),
                                 supercritical=int(sup.sum()), behind_shell=int(inside.sum()),
                                 shell_vn=_summ(rd["vn"]), distance=_summ(near_d),
                                 height_range=dict(min=float(np.nanmin(hr["r_min"])), max=float(np.nanmax(hr["r_max"]))),
                                 tracks=tracks, polarisation=_summ(np.where(pol["dcp_valid"], pol["dcp_source"], np.nan)),
                                 polarisation_invalid=int((np.isfinite(pol["dcp_source"]) & ~pol["dcp_valid"]).sum()))
            surfs = [surface(k) for k in scale]
            if standoff:                                       # the ahead/behind test against the shock surface
                cols = dict(shock_distance_rsun=[], behind_shock_front=[], shock_delta_rsun=[], shock_vn=[], shock_theta_bn=[],
                            shock_mach_f=[], shock_mach_c=[], shock_mach_f_mc=[])
                for q, S in zip(pts, surfs):
                    _, sd, i = S.nearest(q)
                    dd = S.at_points(q[None], i, plasma)
                    for c, val in (("shock_distance_rsun", abs(sd[0])), ("behind_shock_front", bool(sd[0] < 0)),
                                   ("shock_delta_rsun", S.delta[i[0]]), ("shock_vn", dd["vn"][0]), ("shock_theta_bn", dd["theta_bn"][0]),
                                   ("shock_mach_f", dd["mach_f"][0]), ("shock_mach_c", dd["mach_c"][0]),
                                   ("shock_mach_f_mc", dd["mach_f_mc"][0])):
                        cols[c].append(val)
                rtab.update({c: np.array(val) for c, val in cols.items()})
                summ["radio"]["shock_surface"] = dict(behind=int(np.sum(rtab["behind_shock_front"])),
                                                      distance=_summ(rtab["shock_distance_rsun"]),
                                                      mach_f=_summ(rtab["shock_mach_f"]), mach_f_mc=_summ(rtab["shock_mach_f_mc"]),
                                                      theta_bn=_summ(rtab["shock_theta_bn"]))
            if trace and s["connectivity"]:                    # field lines of the sources to the shock surface
                rows = sp.connectivity(m, p, pts, surfs, plasma, s)
                for j, row in enumerate(rows):
                    i = j // 2
                    rows[j] = dict(source=i, instrument=str(rt["instrument"][i]), freq_mhz=float(rt["freq_mhz"][i]),
                                   time_utc=str(rt["time_utc"][i]), **row)
                tables["shock_connectivity"] = rows

                def hit(i, ok):
                    return any(ok(x) for x in rows[2 * i:2 * i + 2] if x["meets"])
                rtab["connected"] = np.array([hit(i, lambda x: True) for i in range(len(pts))])
                rtab["connected_qperp_supercritical"] = np.array([hit(i, lambda x: x["quasi_perp"] and x["supercritical"])
                                                                  for i in range(len(pts))])
                met = [x for x in rows if x["meets"]]
                freqs = np.asarray(rt["freq_mhz"], float)
                summ["connectivity"] = dict(
                    n_sources=int(len(pts)), connected=int(rtab["connected"].sum()),
                    qperp_supercritical=int(rtab["connected_qperp_supercritical"].sum()),
                    path_rsun=_summ([x["path_rsun"] for x in met]), theta_bn=_summ([x["theta_bn"] for x in met]),
                    mach_f=_summ([x["mach_f"] for x in met]), mach_f_mc=_summ([x["mach_f_mc"] for x in met]),
                    by_frequency={f"{i} {f:g}": dict(n=int(sel.sum()), connected=int(rtab["connected"][sel].sum()),
                                                     qperp_supercritical=int(rtab["connected_qperp_supercritical"][sel].sum()))
                                  for i, f in sorted({(str(a), float(b)) for a, b in zip(rt["instrument"], freqs)})
                                  for sel in [(rt["instrument"] == i) & (freqs == f)]})

    # collapsing magnetic traps at the time of the reconstruction
    if trace and s["traps"]:
        tl = [ln for ln, a in zip(s_lines, S0.dg["vn"] > 0) if a]
        tr, tdist = sp.traps(tl, S0, plasma, s, pts if len(pts) else None)
        for t in tr:
            t.setdefault("n_radio_near", 0)
            t.setdefault("radio_min_distance_rsun", np.nan)
        tables["shock_traps"] = [dict(trap=i, xa=t["a"][0], ya=t["a"][1], za=t["a"][2], xb=t["b"][0], yb=t["b"][1], zb=t["b"][2],
                                      **{k: v for k, v in t.items() if k not in ("a", "b")}) for i, t in enumerate(tr)]
        if tdist is not None and rtab:
            rtab["trap_distance_rsun"] = tdist
        both = [t for t in tr if np.nan_to_num(t["mach_f_a"]) > 1 and np.nan_to_num(t["mach_f_b"]) > 1]
        summ["traps"] = dict(n=len(tr), shocked_both=len(both), length_rsun=_summ([t["length_rsun"] for t in tr]),
                             top_r=_summ([t["top_r"] for t in tr]),
                             theta_bn=_summ([x for t in tr for x in (t["theta_bn_a"], t["theta_bn_b"])]),
                             radio_on_trap=int(np.sum(tdist <= s["trap_near_rsun"])) if tdist is not None else 0,
                             near_rsun=s["trap_near_rsun"])

    # time evolution through the burst
    if s["evolution"]:
        tb = sorted(times) if len(times) else [t0]
        burst = dict(radio_start=tb[0].isoformat(), radio_end=tb[-1].isoformat())
        burst.update(swpc_burst(home, t0) or {})
        ts = sp.times_through(tb[0], tb[-1], t0, h, v_apex, s, extra=set(times))
        rows = [dict(time=t.isoformat(), dt_s=(t - t0).total_seconds(),
                     **sp.evolution_row(surface((h + v_apex * (t - t0).total_seconds() / R_SUN_KM) / h), axis, s)) for t in ts]
        tables["shock_evolution"] = rows
        cols = ("apex_mach_f", "apex_mach_f_mc", "nose_mach_f", "nose_mach_f_mc", "flank_mach_f", "flank_mach_f_mc")
        side = lambda x: None if not np.isfinite(x) else ("above" if x >= 1 else "below")
        summ["evolution"] = dict(
            n_times=len(rows), start=ts[0].isoformat() if ts else None, end=ts[-1].isoformat() if ts else None,
            step_s=s["evolution_step_s"], burst=burst,
            crossings={c: [dict(time=x["time"].isoformat(), direction=x["direction"])
                           for x in sp.threshold_crossings(ts, [r[c] for r in rows])] for c in cols},
            state_at_start={c: side(rows[0][c]) for c in cols} if rows else {},
            state_at_end={c: side(rows[-1][c]) for c in cols} if rows else {})

    # uncertainty: ensemble over apex height, speed, density scale, and model runs
    if int(s["mc_samples"]) > 0:
        summ["uncertainty"], tables["shock_mc"], tables["shock_mc_radio"] = ensemble(
            home, p, s, rings, normals, keep2, plasma, m, h, v_apex, speed is not None, axis, t0, pts, times, rtab, log)
    return summ, shell, rtab, tracks


def swpc_burst(home, t0):
    """Begin and end of the type II burst of the SWPC list of the day (cached in data/bursts) that
    contains t0, or else begins nearest to it; None without the cached list."""
    from . import bursts
    f = Path(home) / "data" / "bursts" / f"{t0:%Y%m%d}events.txt"
    if not f.exists():
        return None
    rows = [r for r in bursts.parse(f.read_text(encoding="utf-8", errors="replace"), t0.date()) if r["type"] == "II"]
    if not rows:
        return None
    inside = [r for r in rows if r["begin"] <= t0 <= (r["end"] or r["begin"])]
    r = inside[0] if inside else min(rows, key=lambda r: abs((r["begin"] - t0).total_seconds()))
    return dict(type_ii_begin=r["begin"].isoformat(), type_ii_end=r["end"].isoformat() if r["end"] else None,
                type_ii_observatory=r["obs"], type_ii_speed_km_s=r["speed_km_s"])


def neighbour_runs(home, p, s, plasma):
    """[(run name, Plasma)] of the run of the configuration and of the other runs of the ensemble:
    shock.mc_runs "auto" takes every other folder of data/ that holds the files of a run and its
    temperature (at most two), a list takes those runs, and [] none."""
    from .paths import model_dir
    out = [(p["run"], plasma)]
    here = model_dir(home, p)
    if s["mc_runs"] == "auto":
        cands = [d for d in sorted((Path(home) / "data").glob("*")) if d.is_dir() and d != here
                 and all((d / f"{f}.hdf").exists() for f in mast_model.FILES)][:2]
    else:
        cands = [model_dir(home, dict(p, run=r)) for r in (s["mc_runs"] or [])]
    for d in cands:
        try:
            m2 = mast_model.load(d, density_scale=p["density_scale"], va_density=p["va_density"])
            if "t002" not in mast_model.load_plasma(d, m2) and s["temperature_k"] is None:
                print(f"ensemble: run {d.name} left out (no t002)")
                continue
            out.append((d.name, Plasma(m2, p, s["temperature_k"], s["flow"])))
        except Exception as e:                                # a run that cannot be read is left out
            print(f"ensemble: run {d.name} left out ({e})")
    return out


def ensemble(home, p, s, rings, normals, keep, plasma, m, h, v, speed_given, axis, t0, pts, times, rtab, log=print):
    """Ensemble of shock.mc_samples members over the apex height (normal, 1 sigma from the fit or
    mc_height_err_rsun), the apex speed (normal, mc_speed_frac), the density scale (uniform in
    mc_density_scale), and the model runs of neighbour_runs. Returns (summary, member rows, radio rows)."""
    from . import shockphys as sp
    rng = np.random.default_rng(int(s["mc_seed"]))
    fit = p.get("gcs_fit") or {}
    f = Path(home) / "gcs_fit.json"
    rec = (json.loads(f.read_text(encoding="utf-8")).get("fit") or {}) if f.exists() else {}
    err = fit.get("height_t_out_err") if fit else rec.get("height_t_out_err")
    if s["mc_height_err_rsun"] is not None:
        sh, sh_src = float(s["mc_height_err_rsun"]), "shock.mc_height_err_rsun"
    elif err:
        sh, sh_src = float(err), "error of the apex height of the fit"
    else:
        sh, sh_src = 0.15, "0.15 R_sun (assumed; the fit gives no error)"
    sv = float(s["mc_speed_frac"]) * v
    runs = neighbour_runs(home, p, s, plasma)
    lo, hi = s["mc_density_scale"]
    n = int(s["mc_samples"])
    standoff = bool(s["standoff"])
    rows, rad = [], []
    for i in range(n):
        h_i = max(float(rng.normal(h, sh)), float(s["evolution_min_height"]))
        v_i = max(float(rng.normal(v, sv)), 1.0)
        ds = float(rng.uniform(lo, hi))
        ri = int(rng.integers(len(runs)))
        pl = sp.ScaledPlasma(runs[ri][1], ds / p["density_scale"], p["va_density"] == "scaled")
        S = sp.surface_state(rings, normals, keep, h_i / h, pl, v_i, h, s, standoff)
        rows.append(dict(member=i, speed_km_s=v_i, density_scale=ds, run=runs[ri][0], **sp.evolution_row(S, axis, s)))
        if len(pts):
            ks = np.array([(h_i + v_i * (t - t0).total_seconds() / R_SUN_KM) / h for t in times])
            for k in np.unique(ks):
                sel = np.flatnonzero(ks == k)
                St = sp.surface_state(rings, normals, keep, k, pl, v_i, h, s, standoff)
                _, sd, idx = St.nearest(pts[sel])
                dd = St.at_points(pts[sel], idx, pl)
                for j, a in enumerate(sel):
                    rad.append(dict(member=i, source=int(a), signed_distance_rsun=float(sd[j]), ahead=bool(sd[j] >= 0),
                                    theta_bn=float(dd["theta_bn"][j]), mach_f=float(dd["mach_f"][j]),
                                    mach_c=float(dd["mach_c"][j]), mach_f_mc=float(dd["mach_f_mc"][j])))
        if (i + 1) % 10 == 0:
            log(f"shock: ensemble member {i + 1} of {n}")
    keys = [k for k in rows[0] if k not in ("member", "run")] if rows else []
    summ = dict(n=n, seed=int(s["mc_seed"]), height_sigma_rsun=sh, height_sigma_source=sh_src, speed_sigma_km_s=sv,
                speed_sigma_source=f"{s['mc_speed_frac']:g} of the apex speed" + (" (given with --speed)" if speed_given else ""),
                density_scale=[lo, hi], runs=[r for r, _ in runs], surface="standoff" if standoff else "shell",
                quantities={k: sp.quantiles([r[k] for r in rows]) for k in keys})
    src = []
    for a in range(len(pts) if rad else 0):
        rr = [x for x in rad if x["source"] == a]
        src.append(dict(source=a, instrument=str(rtab["instrument"][a]), freq_mhz=float(rtab["freq_mhz"][a]),
                        time_utc=str(rtab["time_utc"][a]), p_ahead=float(np.mean([x["ahead"] for x in rr])),
                        **{f"{k}_{q}": sp.quantiles([x[k] for x in rr]).get(q, np.nan)
                           for k in ("signed_distance_rsun", "theta_bn", "mach_f", "mach_f_mc") for q in ("median", "p16", "p84")}))
    if src:
        summ["radio"] = src
    return summ, rows, src


# ------------------------------------------------------------------ output
SHELL_COLUMNS = ("x", "y", "z", "r", "lon_carrington_deg", "lat_deg", "nx", "ny", "nz", "normal_axis_deg", "region", "field",
                 "ne", "b", "t", "va", "cs", "vf", "theta_bn", "u_n", "vn", "mach_a", "mach_f", "mach_s",
                 "beta", "mach_c", "mach_f_mc")
UNITS = dict(x="R_sun", y="R_sun", z="R_sun", r="R_sun", ne="cm^-3", b="G", t="K", va="km/s", cs="km/s", vf="km/s",
             theta_bn="deg", u_n="km/s", vn="km/s", normal_axis_deg="deg", delta_rsun="R_sun", rc_rsun="R_sun",
             path_rsun="R_sun", length_rsun="R_sun", top_r="R_sun", height_rsun="R_sun", dt_s="s")
TABLES = ("shock_surface", "shock_connectivity", "shock_traps", "shock_evolution", "shock_mc", "shock_mc_radio")


def write(out, summ, shell, rtab, tracks, tables=None):
    """analysis/shock_shell.csv, shock_radio.csv, radio_tracks.csv, shock.npz, shock_summary.json, and
    the tables of `tables` (TABLES; a dictionary of columns or a list of rows) as <name>.csv."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)

    def table(path, data, cols):
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(cols)
            w.writerows(zip(*(np.asarray(data[c]).tolist() for c in cols)))
    table(out / "shock_shell.csv", shell, SHELL_COLUMNS)
    if rtab:
        table(out / "shock_radio.csv", rtab, list(rtab))
    if tracks:
        with open(out / "radio_tracks.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, list(tracks[0]))
            w.writeheader()
            w.writerows(tracks)
    for name in TABLES:
        data = (tables or {}).get(name)
        if isinstance(data, dict) and data:
            table(out / f"{name}.csv", data, list(data))
        elif data:
            with open(out / f"{name}.csv", "w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, list(dict.fromkeys(k for row in data for k in row)))
                w.writeheader()
                w.writerows([{k: (v.item() if isinstance(v, np.generic) else v) for k, v in row.items()} for row in data])
    arrays = {f"shell_{k}": np.asarray(v) for k, v in shell.items()}
    arrays.update({f"radio_{k}": np.asarray(v) for k, v in rtab.items()})
    if isinstance((tables or {}).get("shock_surface"), dict):
        arrays.update({f"surface_{k}": np.asarray(v) for k, v in tables["shock_surface"].items()})
    np.savez_compressed(out / "shock.npz", **arrays)
    (out / "shock_summary.json").write_text(json.dumps(dict(summ, units=UNITS), indent=1, default=_json) + "\n", encoding="utf-8")
    return out


def _json(x):
    if isinstance(x, (np.floating, np.integer)):
        return x.item()
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, datetime):
        return x.isoformat()
    return str(x)


def _f(x, fmt=".2f"):
    return "n/a" if x is None or (isinstance(x, float) and not np.isfinite(x)) else format(x, fmt)


def report_lines(summ):
    """Lines of the section 'Shock diagnostics' of event_report.md."""
    c, k, a = summ["conventions"], summ["kinematics"], summ["apex"]
    L = [f"- Shell at {summ['time']} with apex height {summ['height_rsun']:g} R_sun and apex speed {summ['apex_speed_km_s']:g} km/s "
         f"({summ['apex_speed_source']}). The shell expands self-similarly about Sun centre, so that its half-width of "
         f"{k['half_width_rsun']:.2f} R_sun grows at {k['lateral_speed_km_s']:.0f} km/s.",
         f"- Upstream plasma from the MAS/MAST run (no CME in the run): Alfven speed from the {c['va_density']} density "
         f"(density_scale {c['density_scale']:g}), temperature {c['temperature']}, flow {c['flow']}, gamma {c['gamma']:.3f}.",
         f"- Apex: v_A {_f(a['va'], '.0f')} km/s, c_s {_f(a['cs'], '.0f')} km/s, v_f {_f(a['vf'], '.0f')} km/s, theta_Bn "
         f"{_f(a['theta_bn'], '.0f')} deg, M_A {_f(a['mach_a'])}, M_f {_f(a['mach_f'])} (field line {a['field'] or 'not traced'})."]
    for name in ("nose", "flank", "closed_field", "open_field"):
        st = summ.get(name)
        if not st or not st["mach_a"].get("n"):
            continue
        L.append(f"- {name.replace('_', ' ').capitalize()} points ({st['mach_a']['n']} advancing): normal speed median "
                 f"{st['vn']['median']:.0f} km/s, v_A median {st['va']['median']:.0f} km/s (min {st['va']['min']:.0f}), "
                 f"M_A median {st['mach_a']['median']:.2f} (max {st['mach_a']['max']:.2f}), M_f median {st['mach_f'].get('median', float('nan')):.2f} "
                 f"(max {st['mach_f'].get('max', float('nan')):.2f}), theta_Bn median {st['theta_bn']['median']:.0f} deg, "
                 f"quasi-perpendicular fraction {_f(st['quasi_perpendicular_fraction'])}"
                 + (f", supercritical fraction (M_f > M_c) {_f(st['supercritical_fraction'])}, quasi-perpendicular and supercritical "
                    f"{_f(st['qperp_supercritical_fraction'])}" if "supercritical_fraction" in st else "") + ".")
    if a.get("mach_c") is not None:
        L.append(f"- First critical Mach number at the apex (Edmiston & Kennel 1984, from the Rankine-Hugoniot relations): "
                 f"M_c {_f(a['mach_c'])} at beta {_f(a['beta'], '.3g')}, so that M_f / M_c is {_f(a['mach_f_mc'])}.")
    L += physics_lines(summ)
    if summ.get("figure"):
        L.append(f"- Figure: {summ['figure']}, the shell coloured by the fast-magnetosonic Mach number (grey: not advancing).")
    r = summ.get("radio")
    if r:
        L.append(f"- Radio sources ({r['n']}): upstream v_A median {_f(r['va'].get('median'), '.0f')} km/s, M_A median "
                 f"{_f(r['mach_a'].get('median'))} (range {_f(r['mach_a'].get('min'))} to {_f(r['mach_a'].get('max'))}), M_f median "
                 f"{_f(r['mach_f'].get('median'))}, theta_Bn median {_f(r['theta_bn'].get('median'), '.0f')} deg, distance to the shell "
                 f"median {_f(r['distance'].get('median'))} R_sun.")
        if "mach_f_mc" in r:
            L.append(f"- At the radio sources M_c has a median of {_f(r['mach_c'].get('median'))} and M_f / M_c a median of "
                     f"{_f(r['mach_f_mc'].get('median'))}; {r['supercritical']} of {r['n']} sources face a supercritical front, and "
                     f"{r['behind_shell']} lie behind the shell at their time.")
        ss = r.get("shock_surface")
        if ss:
            L.append(f"- Against the stand-off shock surface, {ss['behind']} of {r['n']} sources lie behind the shock and "
                     f"{r['n'] - ss['behind']} ahead of it (distance median {_f(ss['distance'].get('median'))} R_sun); M_f median "
                     f"{_f(ss['mach_f'].get('median'))}, M_f / M_c median {_f(ss['mach_f_mc'].get('median'))}.")
        L.append(f"- Radio heights over the density models: {r['height_range']['min']:.2f} to {r['height_range']['max']:.2f} R_sun.")
        for t in r["tracks"]:
            L.append(f"- Radio track {t['instrument']} {t['freq_mhz']:g} MHz ({t['n']} positions, {t['start'][11:19]} to "
                     f"{t['end'][11:19]} UT): {_f(t['speed_3d_km_s'], '.0f')} +/- {_f(t['speed_3d_err_km_s'], '.0f')} km/s in three "
                     f"dimensions, {_f(t['speed_sky_km_s'], '.0f')} +/- {_f(t['speed_sky_err_km_s'], '.0f')} km/s in the plane of the sky.")
        if r["polarisation"].get("n"):
            L.append(f"- Degree of circular polarisation V/I at the sources: median {r['polarisation']['median']:+.2f} "
                     f"(range {r['polarisation']['min']:+.2f} to {r['polarisation']['max']:+.2f}, {r['polarisation']['n']} sources).")
        if r.get("polarisation_invalid"):
            L.append(f"- {r['polarisation_invalid']} sources have |V/I| > 1, which is not physical (Stokes V artefacts larger "
                     "than the source); they are flagged dcp_valid false in shock_radio.csv and left out above.")
    return L


def physics_lines(summ):
    """Report lines of the stand-off surface, connectivity, traps, evolution, and ensemble."""
    L = []
    so = summ.get("standoff")
    if so:
        ap = so["apex"]
        L.append(f"- Stand-off shock surface ({so['relation']}; R_c from the {so['radius']} curvature of the shell): "
                 f"{so['apex_delta_rsun']:.3f} R_sun ahead of the apex (R_c {_f(so['apex_rc_rsun'])} R_sun, driver M_f "
                 f"{_f(so['apex_driver_mach_f'])}); {so['n_shifted']} points moved (median {_f(so['delta'].get('median'), '.3f')} R_sun, "
                 f"{so['n_capped']} at the limit of {so['max_rc']:g} R_c). On the shock surface the apex lies at {ap['r']:.2f} R_sun with "
                 f"M_f {_f(ap['mach_f'])}, theta_Bn {_f(ap['theta_bn'], '.0f')} deg, and M_f / M_c {_f(ap['mach_f_mc'])}.")
        for name in ("nose", "flank", "closed_field", "open_field"):
            st = so.get(name)
            if st and st["mach_f"].get("n"):
                L.append(f"- Shock surface, {name.replace('_', ' ')} ({st['mach_f']['n']} advancing): M_f median {st['mach_f']['median']:.2f} "
                         f"(max {st['mach_f']['max']:.2f}), theta_Bn median {st['theta_bn']['median']:.0f} deg, supercritical fraction "
                         f"{_f(st['supercritical_fraction'])} (quasi-perpendicular and supercritical {_f(st['qperp_supercritical_fraction'])}).")
    c = summ.get("connectivity")
    if c:
        L.append(f"- Magnetic connectivity: the field lines of {c['connected']} of {c['n_sources']} radio sources meet the shock "
                 f"surface ({summ.get('surface', 'shell')}) within a path of median {_f(c['path_rsun'].get('median'))} R_sun, where theta_Bn has "
                 f"a median of {_f(c['theta_bn'].get('median'), '.0f')} deg and M_f / M_c a median of {_f(c['mach_f_mc'].get('median'))}; "
                 f"{c['qperp_supercritical']} sources connect to a quasi-perpendicular supercritical part of the shock.")
        for k, v in c["by_frequency"].items():
            L.append(f"  - {k} MHz: {v['connected']} of {v['n']} connected, {v['qperp_supercritical']} to a quasi-perpendicular "
                     "supercritical shock.")
    t = summ.get("traps")
    if t:
        L.append(f"- Collapsing magnetic traps (closed loops that the shock crosses twice, loop top ahead of the shock): {t['n']}, "
                 f"{t['shocked_both']} with M_f > 1 at both crossings; loop length between the crossings median "
                 f"{_f(t['length_rsun'].get('median'))} R_sun (range {_f(t['length_rsun'].get('min'))} to {_f(t['length_rsun'].get('max'))}), "
                 f"loop top median {_f(t['top_r'].get('median'))} R_sun; {t['radio_on_trap']} radio sources lie within "
                 f"{t['near_rsun']:g} R_sun of a trap.")
    e = summ.get("evolution")
    if e and e["n_times"]:
        b = e["burst"]
        end = b.get("type_ii_end") or "unknown"
        L.append(f"- Time evolution ({e['n_times']} times from {e['start'][11:19]} to {e['end'][11:19]} UT; radio sources "
                 f"{b['radio_start'][11:19]} to {b['radio_end'][11:19]} UT"
                 + (f", SWPC type II {b['type_ii_begin'][11:16]} to {end[11:16]} UT" if b.get("type_ii_begin") else "") + "):")
        names = dict(apex_mach_f="apex M_f = 1", apex_mach_f_mc="apex M_f = M_c", nose_mach_f="nose median M_f = 1",
                     nose_mach_f_mc="nose median M_f = M_c", flank_mach_f="flank median M_f = 1", flank_mach_f_mc="flank median M_f = M_c")
        for k, lab in names.items():
            cr = e["crossings"].get(k) or []
            when = ", ".join(f"{x['time'][11:19]} UT ({x['direction']})" for x in cr) or (
                f"not crossed (always {e['state_at_start'].get(k)})")
            L.append(f"  - {lab}: {when}.")
    u = summ.get("uncertainty")
    if u:
        q = u["quantities"]

        def rng(k, fmt=".2f"):
            x = q.get(k) or {}
            return f"{_f(x.get('median'), fmt)} ({_f(x.get('p16'), fmt)} to {_f(x.get('p84'), fmt)})" if x.get("n") else "n/a"
        L.append(f"- Ensemble of {u['n']} members (apex height sigma {u['height_sigma_rsun']:.3f} R_sun, {u['height_sigma_source']}; speed "
                 f"sigma {u['speed_sigma_km_s']:.0f} km/s; density_scale {u['density_scale'][0]:g} to {u['density_scale'][1]:g}; runs "
                 f"{', '.join(u['runs'])}), median and 16 to 84 % range: apex M_f {rng('apex_mach_f')}, apex M_f / M_c "
                 f"{rng('apex_mach_f_mc')}, nose M_f {rng('nose_mach_f')}, flank M_f {rng('flank_mach_f')}, fraction shocked "
                 f"{rng('frac_shock')}, supercritical {rng('frac_supercritical')}, quasi-perpendicular supercritical "
                 f"{rng('frac_qperp_supercritical')}" + (f", stand-off at the apex {rng('apex_standoff_rsun', '.3f')} R_sun"
                                                        if u["surface"] == "standoff" else "") + ".")
        for x in u.get("radio", []):
            L.append(f"  - {x['instrument']} {x['freq_mhz']:g} MHz {x['time_utc'][11:19]}: ahead of the shock in "
                     f"{100 * x['p_ahead']:.0f} % of members, M_f {_f(x['mach_f_median'])} ({_f(x['mach_f_p16'])} to {_f(x['mach_f_p84'])}), "
                     f"theta_Bn {_f(x['theta_bn_median'], '.0f')} deg, M_f / M_c {_f(x['mach_f_mc_median'])} "
                     f"({_f(x['mach_f_mc_p16'])} to {_f(x['mach_f_mc_p84'])}).")
    return L


def update_report(home, summ):
    """Write (or replace) the section 'Shock diagnostics' of event_report.md, when the report exists."""
    f = Path(home) / "event_report.md"
    if not f.exists():
        return None
    text = f.read_text(encoding="utf-8")
    sec = "\n".join(["## Shock diagnostics", ""] + report_lines(summ)) + "\n"
    if "## Shock diagnostics" in text:
        head, rest = text.split("## Shock diagnostics", 1)
        nxt = rest.find("\n## ")
        text = head + sec + (rest[nxt + 1:] if nxt >= 0 else "")
    else:
        marker = "## Figures"
        text = text.replace(marker, sec + "\n" + marker, 1) if marker in text else text + "\n" + sec
    f.write_text(text, encoding="utf-8")
    return f


def run(home, cfg, a):
    """The `corona3d shock` command."""
    over = {k: v for k, v in (("standoff", True if getattr(a, "standoff", False) else None),
                              ("standoff_radius", getattr(a, "standoff_radius", None)),
                              ("mc_samples", getattr(a, "mc", None)), ("mc_seed", getattr(a, "mc_seed", None)),
                              ("connectivity", False if getattr(a, "no_connectivity", False) else None),
                              ("traps", False if getattr(a, "no_traps", False) else None),
                              ("evolution", False if getattr(a, "no_evolution", False) else None),
                              ("evolution_step_s", getattr(a, "evolution_step", None))) if v is not None}
    tables = {}
    summ, shell, rtab, tracks = analyse(home, p={"shock": over} if over else None, speed=a.speed, speed_note=a.speed_note,
                                        trace=not a.no_trace, radio=not a.no_radio, tables=tables)
    out = write(Path(home) / a.out, summ, shell, rtab, tracks, tables)
    lines = report_lines(summ)
    (out / "shock_summary.md").write_text("# Shock diagnostics\n\n" + "\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    rep = update_report(home, summ)
    print(f"wrote {out}" + (f" and the section Shock diagnostics of {rep}" if rep else ""))
