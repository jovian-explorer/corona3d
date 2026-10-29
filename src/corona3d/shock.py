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


def analyse(home, p=None, speed=None, speed_note=None, trace=True, radio=True, log=print):
    """Shock diagnostics, kinematics, radio heights, and polarisation of a working directory.
    Returns (summary dict, shell table dict, radio table dict, track list)."""
    from . import heights, nrh
    from .geometry import pythea_shell
    from .paths import config, merge, model_dir
    from .whitelight import utc
    p = merge(config(home), p)
    s = settings(p)
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
    xyz, nrm = rings.reshape(-1, 3), normals.reshape(-1, 3)
    keep = np.linalg.norm(xyz, axis=1) > s["r_min"]
    xyz, nrm = xyz[keep], nrm[keep]
    h = float(g["height"])
    dg = diagnostics(xyz, nrm, plasma, v_apex, h, s["gamma"])
    axis = axis_direction(g, p)
    ang = np.degrees(np.arccos(np.clip(nrm @ axis, -1, 1)))
    region = np.where(ang <= s["nose_deg"], "nose", np.where(ang >= s["flank_deg"], "flank", "shoulder"))
    kind = field_class(m, p, xyz) if trace else np.full(len(xyz), "")
    from .api import carrington_angles
    lon, lat, r = carrington_angles(xyz)
    shell = dict(x=xyz[:, 0], y=xyz[:, 1], z=xyz[:, 2], r=r, lon_carrington_deg=lon, lat_deg=lat,
                 nx=nrm[:, 0], ny=nrm[:, 1], nz=nrm[:, 2], normal_axis_deg=ang, region=region, field=kind, **dg)
    apex = int(np.argmax(r))
    w = half_width(rings, axis)
    adv = dg["vn"] > 0
    qperp = dg["theta_bn"] > s["quasi_perp_deg"]

    def stats(sel):
        return {k: _summ(dg[k][sel]) for k in ("vn", "va", "cs", "vf", "theta_bn", "mach_a", "mach_f")} | dict(
            quasi_perpendicular_fraction=float(np.mean(qperp[sel & adv])) if (sel & adv).any() else None)

    summ = dict(
        time=p["time"], height_rsun=h, apex_speed_km_s=v_apex, apex_speed_source=v_src,
        conventions=dict(va_density=p["va_density"], density_scale=p["density_scale"], gamma=s["gamma"],
                         temperature="MAS t002" if plasma.has_t else f"constant {s['temperature_k']} K (assumed)",
                         flow="MAS vr002, vt002, vp002 subtracted" if plasma.has_flow else "not subtracted",
                         expansion="self-similar GCS shell about Sun centre: V(x) = (V_apex / h) x",
                         nose_deg=s["nose_deg"], flank_deg=s["flank_deg"], quasi_perp_deg=s["quasi_perp_deg"]),
        kinematics=dict(half_width_rsun=w, lateral_speed_km_s=v_apex * w / h,
                        flank_normal_speed_km_s=_summ(dg["vn"][region == "flank"]),
                        shells=[dict(time=f["time"], height_rsun=f["height"], half_width_rsun=f["height"] * w / h) for f in frames]),
        apex={k: (float(dg[k][apex]) if np.ndim(dg[k]) == 1 else None) for k in ("vn", "ne", "b", "t", "va", "cs", "vf", "theta_bn", "mach_a", "mach_f", "mach_s")}
        | dict(r=float(r[apex]), lon_carrington_deg=float(lon[apex]), lat_deg=float(lat[apex]), field=str(kind[apex])),
        nose=stats(region == "nose"), flank=stats(region == "flank"), all=stats(np.ones(len(r), bool)),
        n_points=int(len(r)), n_advancing=int(adv.sum()))
    if trace:
        summ["closed_field"] = stats(kind == "closed")
        summ["open_field"] = stats(kind == "open")
    log(f"shock: apex speed {v_apex:g} km/s ({v_src}); height {h:g} R_sun at {p['time']}; {adv.sum()} of {len(r)} shell points advance")

    # radio sources
    rtab, tracks = {}, []
    if radio:
        from .api import build
        rec = build(home, p, surface=False, open=False, closed=False, streamer=False, gcs=False, radio=True)
        rt = rec.radio
        if rt and len(rt["x"]):
            pts = np.c_[rt["x"], rt["y"], rt["z"]]
            t0 = utc(p["time"])
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
            rd = diagnostics(xyz[near_i], nrm[near_i], plasma, v_apex, h, s["gamma"], b_at=plasma.at(pts))
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
            summ["radio"] = dict(n=int(len(pts)), mach_a=_summ(rd["mach_a"]), mach_f=_summ(rd["mach_f"]),
                                 va=_summ(rd["va"]), vf=_summ(rd["vf"]), theta_bn=_summ(rd["theta_bn"]),
                                 shell_vn=_summ(rd["vn"]), distance=_summ(near_d),
                                 height_range=dict(min=float(np.nanmin(hr["r_min"])), max=float(np.nanmax(hr["r_max"]))),
                                 tracks=tracks, polarisation=_summ(np.where(pol["dcp_valid"], pol["dcp_source"], np.nan)),
                                 polarisation_invalid=int((np.isfinite(pol["dcp_source"]) & ~pol["dcp_valid"]).sum()))
    return summ, shell, rtab, tracks


# ------------------------------------------------------------------ output
SHELL_COLUMNS = ("x", "y", "z", "r", "lon_carrington_deg", "lat_deg", "nx", "ny", "nz", "normal_axis_deg", "region", "field",
                 "ne", "b", "t", "va", "cs", "vf", "theta_bn", "u_n", "vn", "mach_a", "mach_f", "mach_s")
UNITS = dict(x="R_sun", y="R_sun", z="R_sun", r="R_sun", ne="cm^-3", b="G", t="K", va="km/s", cs="km/s", vf="km/s",
             theta_bn="deg", u_n="km/s", vn="km/s", normal_axis_deg="deg")


def write(out, summ, shell, rtab, tracks):
    """analysis/shock_shell.csv, shock_radio.csv, radio_tracks.csv, shock.npz, shock_summary.json."""
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
    arrays = {f"shell_{k}": np.asarray(v) for k, v in shell.items()}
    arrays.update({f"radio_{k}": np.asarray(v) for k, v in rtab.items()})
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
                 f"quasi-perpendicular fraction {_f(st['quasi_perpendicular_fraction'])}.")
    if summ.get("figure"):
        L.append(f"- Figure: {summ['figure']}, the shell coloured by the fast-magnetosonic Mach number (grey: not advancing).")
    r = summ.get("radio")
    if r:
        L.append(f"- Radio sources ({r['n']}): upstream v_A median {_f(r['va'].get('median'), '.0f')} km/s, M_A median "
                 f"{_f(r['mach_a'].get('median'))} (range {_f(r['mach_a'].get('min'))} to {_f(r['mach_a'].get('max'))}), M_f median "
                 f"{_f(r['mach_f'].get('median'))}, theta_Bn median {_f(r['theta_bn'].get('median'), '.0f')} deg, distance to the shell "
                 f"median {_f(r['distance'].get('median'))} R_sun.")
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
    summ, shell, rtab, tracks = analyse(home, speed=a.speed, speed_note=a.speed_note, trace=not a.no_trace,
                                        radio=not a.no_radio)
    out = write(Path(home) / a.out, summ, shell, rtab, tracks)
    lines = report_lines(summ)
    (out / "shock_summary.md").write_text("# Shock diagnostics\n\n" + "\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    rep = update_report(home, summ)
    print(f"wrote {out}" + (f" and the section Shock diagnostics of {rep}" if rep else ""))
