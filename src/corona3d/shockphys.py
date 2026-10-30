"""Shock physics of the CME front beyond the shell itself (`corona3d shock`, METHODS Section 16).

The functions of this module extend the diagnostics of `corona3d.shock`, which treat the GCS shell as
the shock, in six ways:

1. stand-off: a shock surface ahead of the shell (the driver), offset along the normal by the stand-off
   distance of Russell & Mulligan (2002), which follows the relation of Farris & Russell (1994) with the
   radius of curvature R_c of the driver,
       Delta / R_c = 0.81 [(gamma - 1) M^2 + 2] / [(gamma + 1) (M^2 - 1)],
   with M the local fast-magnetosonic Mach number of the driver and R_c from the principal curvatures of
   the GCS shell;
2. first critical Mach number M_c(theta_Bn, beta) of a fast MHD shock (Edmiston & Kennel 1984), the Mach
   number at which the downstream flow speed along the normal equals the downstream sound speed, solved
   from the Rankine-Hugoniot relations, and the ratio M_f / M_c (a shock with M_f > M_c is supercritical);
3. magnetic connectivity of the radio sources: the field line of the model through each source is traced
   both ways to where it meets the shock surface;
4. collapsing magnetic traps (Kong et al. 2015): closed field lines that the shock crosses twice, with the
   loop between the two crossings still ahead of the shock;
5. time evolution of the diagnostics through the burst, with the apex height h(t) = h0 + V (t - t0);
6. uncertainty from an ensemble over apex height, speed, density scale, and neighbouring model runs.

All positions are Carrington Cartesian in R_sun, speeds in km/s, and angles in degrees.
"""
from datetime import timedelta
from functools import lru_cache

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import brentq
from scipy.spatial import cKDTree

R_SUN_KM = 695700.0


# ------------------------------------------------------------------ 1. stand-off distance
def standoff_ratio(mach, gamma=5.0 / 3.0):
    """Delta / R_c of Russell & Mulligan (2002) after Farris & Russell (1994); NaN where M <= 1."""
    m2 = np.asarray(mach, float) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(m2 > 1, 0.81 * ((gamma - 1) * m2 + 2) / ((gamma + 1) * (m2 - 1)), np.nan)


def principal_curvatures(x, n):
    """Largest principal curvature and mean curvature [1/R_sun] of a surface given on a grid (n_axis,
    n_circle, 3), periodic around the circle, with outward unit normals n. Both are positive where the
    surface is convex towards the normal (1/R on a sphere of radius R)."""
    x, n = np.asarray(x, float), np.asarray(n, float)
    xi, ni = np.gradient(x, axis=0), np.gradient(n, axis=0)
    xj = 0.5 * (np.roll(x, -1, axis=1) - np.roll(x, 1, axis=1))
    nj = 0.5 * (np.roll(n, -1, axis=1) - np.roll(n, 1, axis=1))
    dot = lambda a, b: np.sum(a * b, axis=-1)
    E, F, G = dot(xi, xi), dot(xi, xj), dot(xj, xj)
    L, M, N = dot(xi, ni), 0.5 * (dot(xi, nj) + dot(xj, ni)), dot(xj, nj)
    det = np.maximum(E * G - F * F, 1e-30)
    H = (E * N + G * L - 2 * F * M) / (2 * det)
    K = (L * N - M * M) / det
    return H + np.sqrt(np.clip(H * H - K, 0, None)), H


# ------------------------------------------------------------------ 2. Rankine-Hugoniot and M_c
def fast_shock(ma, theta_deg, beta, gamma=5.0 / 3.0):
    """Downstream state of the fast MHD shock with the upstream normal Alfven Mach number ma, angle
    theta_deg between B and the normal, and plasma beta, in units of the upstream density and Alfven
    speed (B / sqrt(4 pi) = 1). Returns dict(X, u2, v2t, bt2, p2, p1, bt1, bn) or None when no fast
    shock exists. X is the root of the shock adiabatic (Priest 2014, Section 5.4) with 1 < X <
    (gamma + 1) / (gamma - 1) and a downstream normal speed above the normal Alfven speed of the field
    (u1^2 > X v_A^2 cos^2 theta), the fast branch."""
    c, s = np.cos(np.radians(theta_deg)), np.sin(np.radians(theta_deg))
    A, C, S = ma * ma, c * c, s * s
    cs2, p1 = gamma * beta / 2, beta / 2
    t1 = np.polymul([C * C, -2 * A * C, A * A], [cs2 + 0.5 * A * (gamma - 1), -0.5 * A * (gamma + 1)])
    t2 = np.polymul([0.5 * A * S, 0.0], [C * (gamma - 1), (2 - gamma) * A - C * (gamma + 1), gamma * A])
    r = np.roots(np.polyadd(t1, t2))
    r = r[np.abs(r.imag) < 1e-9].real
    r = r[(r > 1 + 1e-12) & (r < (gamma + 1) / (gamma - 1)) & (A > r * C * (1 + 1e-12))]
    if not len(r):
        return None
    X = float(r.max())
    bt2 = s * X * (A - C) / (A - X * C)
    p2 = p1 + A * (1 - 1 / X) + (S - bt2**2) / 2
    return dict(X=X, u2=ma / X, v2t=c * (bt2 - s) / ma, bt2=bt2, p2=p2, p1=p1, bt1=s, bn=c)


def _fast_speed(beta, theta_deg, gamma):
    from .shock import fast_speed
    return float(fast_speed(1.0, np.sqrt(gamma * beta / 2), theta_deg))


def critical_mach_exact(theta_deg, beta, gamma=5.0 / 3.0, m_max=20.0):
    """First critical fast-magnetosonic Mach number: the M_f at which the downstream normal flow speed
    equals the downstream sound speed, u_2n^2 = gamma p_2 / rho_2 (Edmiston & Kennel 1984). Returns 1
    when the shock is already subsonic downstream at M_f -> 1, and NaN when no root lies below m_max."""
    th = max(float(theta_deg), 0.5)              # exactly parallel: the fast branch is degenerate
    vf = _fast_speed(beta, th, gamma)

    def f(m):
        d = fast_shock(m * vf, th, beta, gamma)
        return np.nan if d is None else d["u2"] ** 2 - gamma * d["p2"] / d["X"]
    grid = np.r_[1 + 1e-4, np.geomspace(1.01, m_max, 40)]
    vals = np.array([f(m) for m in grid])
    ok = np.isfinite(vals)
    if ok.any() and vals[ok][0] <= 0:
        return 1.0
    i = np.flatnonzero(ok[:-1] & ok[1:] & (np.sign(vals[:-1]) != np.sign(vals[1:])))
    if not len(i):
        return np.nan
    return float(brentq(f, grid[i[0]], grid[i[0] + 1], xtol=1e-6))


def _critical_function(m, th, beta, gamma):
    """u_2n^2 - gamma p_2 / rho_2 (upstream units) of the fast shock at the fast Mach numbers m, angles
    th (0.5 to 89.5 deg), and beta (arrays of one shape), with the cubic of fast_shock solved for all
    at once through the eigenvalues of its companion matrices; NaN where no fast shock exists."""
    from .shock import fast_speed
    m, th, beta = np.broadcast_arrays(*(np.asarray(a, float) for a in (m, th, beta)))
    c, s = np.cos(np.radians(th)), np.sin(np.radians(th))
    C, S = c * c, s * s
    cs2, p1 = gamma * beta / 2, beta / 2
    ma = m * fast_speed(1.0, np.sqrt(cs2), th)
    A = ma * ma
    l1, l0 = cs2 + 0.5 * A * (gamma - 1), -0.5 * A * (gamma + 1)
    q2, q1, q0 = C * C, -2 * A * C, A * A                       # (A - C X)^2
    w = 0.5 * A * S                                             # 0.5 A S X (r2 X^2 + r1 X + r0)
    r2, r1, r0 = C * (gamma - 1), (2 - gamma) * A - C * (gamma + 1), gamma * A
    a3, a2, a1, a0 = q2 * l1 + w * r2, q2 * l0 + q1 * l1 + w * r1, q1 * l0 + q0 * l1 + w * r0, q0 * l0
    comp = np.zeros(m.shape + (3, 3))
    comp[..., 0, 0], comp[..., 0, 1], comp[..., 0, 2] = -a2 / a3, -a1 / a3, -a0 / a3
    comp[..., 1, 0] = comp[..., 2, 1] = 1.0
    ev = np.linalg.eigvals(comp)
    X = np.where((np.abs(ev.imag) < 1e-7 * np.maximum(1, np.abs(ev.real))), ev.real, np.nan)
    X = np.where((X > 1 + 1e-12) & (X < (gamma + 1) / (gamma - 1)) & (A[..., None] > X * C[..., None] * (1 + 1e-12)), X, np.nan)
    with np.errstate(all="ignore"):
        X = np.nanmax(np.where(np.isfinite(X), X, -np.inf), axis=-1)
        X = np.where(np.isfinite(X), X, np.nan)
        bt2 = s * X * (A - C) / (A - X * C)
        p2 = p1 + A * (1 - 1 / X) + (S - bt2**2) / 2
        return (ma / X) ** 2 - gamma * p2 / X


TH_GRID = np.arange(0.0, 90.01, 2.0)
LB_GRID = np.arange(-3.0, 2.001, 0.125)                 # log10 beta


@lru_cache(maxsize=4)
def _critical_table(gamma):
    """M_c on TH_GRID x LB_GRID: the first sign change of _critical_function on 41 Mach numbers from
    1 to 20, refined by 40 bisections (1 where the shock is subsonic downstream at once)."""
    th, lb = np.meshgrid(np.clip(TH_GRID, 0.5, 89.5), LB_GRID, indexing="ij")
    th, beta = th.ravel(), 10**lb.ravel()
    grid = np.r_[1 + 1e-4, np.geomspace(1.01, 20.0, 40)]
    F = _critical_function(grid[None, :], th[:, None], beta[:, None], gamma)
    out = np.full(len(th), np.nan)
    lo, hi = np.full(len(th), np.nan), np.full(len(th), np.nan)
    for n in range(len(th)):
        ok = np.flatnonzero(np.isfinite(F[n]))
        if len(ok) and F[n, ok[0]] <= 0:
            out[n] = 1.0
            continue
        sc = [i for i in range(len(grid) - 1) if np.isfinite(F[n, i]) and np.isfinite(F[n, i + 1]) and F[n, i] * F[n, i + 1] < 0]
        if sc:
            lo[n], hi[n] = grid[sc[0]], grid[sc[0] + 1]
    b = np.isfinite(lo)
    for _ in range(40):
        mid = 0.5 * (lo[b] + hi[b])
        fm = _critical_function(mid, th[b], beta[b], gamma)
        pos = np.nan_to_num(fm, nan=-1.0) > 0
        lo[b] = np.where(pos, mid, lo[b])
        hi[b] = np.where(pos, hi[b], mid)
    out[b] = 0.5 * (lo[b] + hi[b])
    return RegularGridInterpolator((TH_GRID, LB_GRID), out.reshape(len(TH_GRID), len(LB_GRID)),
                                   bounds_error=False, fill_value=None)


def critical_mach(theta_deg, beta, gamma=5.0 / 3.0):
    """M_c at the angles theta_deg and plasma beta (arrays), interpolated in a table of
    critical_mach_exact over 0 to 90 deg and beta 1e-3 to 100 (beta outside is clipped); NaN where
    theta or beta is not finite."""
    th = np.atleast_1d(np.asarray(theta_deg, float))
    b = np.atleast_1d(np.asarray(beta, float))
    th, b = np.broadcast_arrays(th, b)
    ok = np.isfinite(th) & np.isfinite(b) & (b > 0)
    out = np.full(th.shape, np.nan)
    if ok.any():
        q = np.c_[np.clip(th[ok], 0, 90), np.clip(np.log10(b[ok]), LB_GRID[0], LB_GRID[-1])]
        out[ok] = _critical_table(float(gamma))(q)
    return out


def plasma_beta(va, cs, gamma=5.0 / 3.0):
    """beta = 8 pi p / B^2 = (2 / gamma) (c_s / v_A)^2."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return 2.0 / gamma * (np.asarray(cs, float) / np.asarray(va, float)) ** 2


def add_critical(d, gamma=5.0 / 3.0):
    """Add beta, mach_c (M_c), and mach_f_mc (M_f / M_c) to a dictionary of shock.diagnostics."""
    d["beta"] = plasma_beta(d["va"], d["cs"], gamma)
    d["mach_c"] = critical_mach(d["theta_bn"], d["beta"], gamma)
    with np.errstate(divide="ignore", invalid="ignore"):
        d["mach_f_mc"] = d["mach_f"] / d["mach_c"]
    return d


# ------------------------------------------------------------------ the shock surface at one time
class ScaledPlasma:
    """A Plasma whose density is multiplied by `factor` (the density of the Alfven speed too when
    va_density is scaled), for the ensemble over the density scale."""

    def __init__(self, plasma, factor=1.0, scaled_va=True):
        self.p, self.f, self.sv = plasma, float(factor), scaled_va

    def at(self, xyz):
        s = dict(self.p.at(xyz))
        s["ne"] = s["ne"] * self.f
        if self.sv:
            s["n"] = s["n"] * self.f
        return s


class Surface:
    """The shock surface at one time: points xyz (k, 3) with normals nrm, the diagnostics dg of
    shock.diagnostics with M_c, the stand-off delta and radius of curvature rc [R_sun] (zero and NaN for
    the shell), and the apex height ht and speed v of the driver at that time."""

    def __init__(self, xyz, nrm, dg, delta, rc, ht, v, gamma, kind, residual=0.0, apex=None):
        self.xyz, self.nrm, self.dg, self.delta, self.rc = xyz, nrm, dg, delta, rc
        self.apex = int(np.argmax(np.linalg.norm(xyz, axis=1))) if apex is None else apex
        self.ht, self.v, self.gamma, self.kind, self.residual = ht, v, gamma, kind, residual
        self._tree = None

    @property
    def tree(self):
        if self._tree is None:
            self._tree = cKDTree(self.xyz)
        return self._tree

    def nearest(self, pts):
        """(distance to the nearest point, signed distance along its normal, index) of points; the
        signed distance is negative behind the surface (on the side of the driver)."""
        pts = np.atleast_2d(pts)
        d, i = self.tree.query(pts)
        return d, np.sum(self.nrm[i] * (pts - self.xyz[i]), axis=1), i

    def at_points(self, pts, idx, plasma):
        """Diagnostics at the points pts with the normal and normal speed of the surface points idx
        and the upstream plasma sampled at pts."""
        from .shock import diagnostics
        d = diagnostics(self.xyz[idx], self.nrm[idx], plasma, self.v, self.ht, self.gamma, b_at=plasma.at(pts))
        return add_critical(d, self.gamma)


def surface_state(rings, normals, keep, k, plasma, v, h, s, standoff=False):
    """The shock surface when the shell (rings, normals: grids (n_axis, n_circle, 3) without the repeated
    point; keep: points above shock.r_min) is scaled by k about Sun centre, its apex then at h k and
    moving at v. Without stand-off it is the shell; with it, every advancing point whose local
    fast-magnetosonic Mach number exceeds one is moved along its normal by Delta = R_c Delta/R_c(M),
    with M sampled at the moved point (fixed-point iteration with damping), R_c from the curvature of
    the shell (shock.standoff_radius: min = 1 / largest principal curvature, mean = 1 / mean
    curvature), and Delta at most standoff_max_rc R_c. A point of the shock surface keeps the normal of
    the point of the shell below it (the surface parallel to the driver, which holds where Delta varies
    slowly along the front), and it moves self-similarly with (v / (h k)) x."""
    from .shock import diagnostics
    g = s["gamma"]
    ht = h * k
    x0 = np.asarray(rings, float) * k
    delta = np.zeros(keep.shape)
    rc = np.full(keep.shape, np.nan)
    residual = 0.0
    if standoff:
        k1, kh = principal_curvatures(rings, normals)
        kap = k1 if s["standoff_radius"] == "min" else kh
        with np.errstate(divide="ignore"):
            rc = np.where(kap > 1e-6, 1.0 / kap, np.nan) * k
        xs, ns, rk = x0[keep], normals[keep], rc[keep]
        gs = s.get("standoff_gamma") or g
        d = np.zeros(len(xs))
        for it in range(int(s["standoff_iterations"])):
            mf = diagnostics(xs, ns, plasma, v, ht, g, b_at=plasma.at(xs + d[:, None] * ns))["mach_f"]
            ok = np.isfinite(mf) & (mf > 1) & np.isfinite(rk)
            new = np.where(ok, np.minimum(standoff_ratio(np.where(ok, mf, 2.0), gs), s["standoff_max_rc"]) * np.nan_to_num(rk), 0.0)
            residual = float(np.nanmax(np.abs(new - d) / np.where(np.isfinite(rk), rk, np.inf))) if len(d) else 0.0
            d = new if it == 0 else 0.5 * (d + new)
        delta[keep] = d
    xsurf = x0 + delta[..., None] * normals
    xf, nf = xsurf[keep], np.asarray(normals, float)[keep]
    dg = add_critical(diagnostics(xf, nf, plasma, v, ht, g), g)
    return Surface(xf, nf, dg, delta[keep], rc[keep], ht, v, g, "standoff" if standoff else "shell", residual,
                   apex=int(np.argmax(np.linalg.norm(x0[keep], axis=1))))


def summary_masks(S, axis, s):
    """Masks of the advancing points, the nose, and the flanks of a surface (angle of the normal to
    the CME axis), and the index of the apex (the point in front of the apex of the shell)."""
    ang = np.degrees(np.arccos(np.clip(S.nrm @ axis, -1, 1)))
    adv = S.dg["vn"] > 0
    return dict(adv=adv, nose=ang <= s["nose_deg"], flank=ang >= s["flank_deg"], ang=ang,
                apex=S.apex)


def fractions(S, sel, s):
    """Fractions of the advancing points in sel that are shocked (M_f > 1), supercritical (M_f > M_c),
    and quasi-perpendicular and supercritical."""
    sel = sel & (S.dg["vn"] > 0)
    if not sel.any():
        return dict(shock=None, supercritical=None, qperp_supercritical=None)
    mf, r, th = S.dg["mach_f"][sel], S.dg["mach_f_mc"][sel], S.dg["theta_bn"][sel]
    sup = np.nan_to_num(r) > 1
    return dict(shock=float(np.mean(np.nan_to_num(mf) > 1)), supercritical=float(np.mean(sup)),
                qperp_supercritical=float(np.mean(sup & (th > s["quasi_perp_deg"]))))


# ------------------------------------------------------------------ 3. and 4. field lines and the surface
def trace_lines(m, p, seeds, n_steps=1500):
    """Field lines of the model through the seeds (reconstruction frame), traced both ways with the
    tracer of corona3d: list of (points (n, 3), index of the seed on the line, kind)."""
    from scipy.spatial.transform import Rotation
    from .geometry import FieldTracer
    R = Rotation.from_euler("z", p["model_dlon"], degrees=True).as_matrix()
    tr = FieldTracer(m)
    seeds = np.atleast_2d(seeds) @ R
    fw, bw = tr.trace(seeds, sign=+1.0, n_steps=n_steps), tr.trace(seeds, sign=-1.0, n_steps=n_steps)
    out = []
    for a, b in zip(fw, bw):
        pts = np.vstack([b[::-1], a[1:]])
        r0, r1 = np.linalg.norm(pts[[0, -1]], axis=1)
        kind = "open" if (r0 > tr.r_max or r1 > tr.r_max) else "closed" if (r0 < tr.r_min and r1 < tr.r_min) else "undetermined"
        out.append((pts @ R.T, len(b) - 1, kind))
    return out


def mesh_tolerance(S):
    """Distance within which a point counts as near the surface: twice the 90th percentile of the
    distance from a point of the surface to its fourth nearest neighbour."""
    d, _ = S.tree.query(S.xyz, k=min(5, len(S.xyz)))
    return 2.0 * float(np.percentile(d[:, -1], 90))


def crossings(line, S, tol):
    """Crossings of a polyline with the surface: list of dict(j, point, idx, outward) where the signed
    distance changes sign between points j and j + 1 within tol of the surface; outward is True when
    the line goes from behind the surface to ahead of it."""
    d, sd, i = S.nearest(line)
    near = d < tol
    j = np.flatnonzero(near[:-1] & near[1:] & ((sd[:-1] >= 0) != (sd[1:] >= 0)) & (sd[:-1] != sd[1:]))
    out = []
    for a in j:
        t = sd[a] / (sd[a] - sd[a + 1])
        out.append(dict(j=int(a), point=line[a] + t * (line[a + 1] - line[a]), idx=int(i[a] if t < 0.5 else i[a + 1]),
                        outward=bool(sd[a] < 0)))
    return out


def arc(line):
    """Cumulative arc length along a polyline [R_sun]."""
    return np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(line, axis=0), axis=1))]


def connectivity(m, p, pts, surfaces, plasma, s, n_steps=1500):
    """Field line of every radio source traced both ways (along B and against B) to its first crossing of
    the shock surface at the time of the source (surfaces: one Surface per source). Returns rows with
    the path length from the source, the crossing point, and theta_Bn, M_f, M_c, and M_f / M_c there."""
    lines = trace_lines(m, p, pts, n_steps)
    rows = []
    for q, (line, i0, kind), S in zip(pts, lines, surfaces):
        tol = mesh_tolerance(S)
        L = arc(line)
        _, sd0, _ = S.nearest(q)
        for name, seg, sl in (("along_B", line[i0:], L[i0:] - L[i0]), ("against_B", line[i0::-1], L[i0] - L[i0::-1])):
            c = crossings(seg, S, tol) if len(seg) > 1 else []
            row = dict(direction=name, line=kind, source_ahead=bool(sd0[0] >= 0), meets=bool(c))
            if c:
                c0 = c[0]
                dd = S.at_points(c0["point"][None], np.array([c0["idx"]]), plasma)
                th, mf, mc = float(dd["theta_bn"][0]), float(dd["mach_f"][0]), float(dd["mach_c"][0])
                row.update(path_rsun=float(sl[c0["j"]] + np.linalg.norm(c0["point"] - seg[c0["j"]])),
                           x=c0["point"][0], y=c0["point"][1], z=c0["point"][2], r=float(np.linalg.norm(c0["point"])),
                           vn=float(dd["vn"][0]), theta_bn=th, mach_f=mf, mach_c=mc, mach_f_mc=float(dd["mach_f_mc"][0]),
                           quasi_perp=bool(th > s["quasi_perp_deg"]), supercritical=bool(np.nan_to_num(mf / mc) > 1),
                           shocked=bool(np.nan_to_num(mf) > 1), delta_rsun=float(S.delta[c0["idx"]]))
            rows.append(row)
    return rows


def traps(lines, S, plasma, s, radio_pts=None):
    """Collapsing magnetic traps: closed field lines (lines: output of trace_lines) on which a stretch of
    the loop lies ahead of the shock surface between two crossings at advancing points, so that the
    shock reaches the loop at both legs; a stretch that stays within half the mesh tolerance of the
    surface only grazes it and is left out. Returns (traps, distance of every radio point to the nearest
    trap [R_sun]); each trap gives the two crossings, the length of the loop between them, the height
    of its top, and theta_Bn, M_f, and M_f / M_c at both crossings. Traps whose crossings lie within
    shock.trap_merge_rsun of those of a trap already found are merged into it."""
    tol = mesh_tolerance(S)
    adv = S.dg["vn"] > 0
    out, segs = [], []
    for line, _, kind in lines:
        if kind != "closed":
            continue
        c = crossings(line, S, tol)
        for a, b in zip(c[:-1], c[1:]):
            if not (a["outward"] and not b["outward"] and adv[a["idx"]] and adv[b["idx"]]):
                continue
            seg = np.vstack([a["point"], line[a["j"] + 1:b["j"] + 1], b["point"]])
            if S.tree.query(seg)[0].max() < 0.5 * tol:          # a line that only grazes the surface
                continue
            ends = np.array([a["point"], b["point"]])
            dup = False
            for t in out:
                e = np.array([t["a"], t["b"]])
                if (np.linalg.norm(e - ends, axis=1).max() < s["trap_merge_rsun"]
                        or np.linalg.norm(e - ends[::-1], axis=1).max() < s["trap_merge_rsun"]):
                    t["n_lines"] += 1
                    dup = True
                    break
            if dup:
                continue
            dd = S.at_points(ends, np.array([a["idx"], b["idx"]]), plasma)
            out.append(dict(a=a["point"], b=b["point"], length_rsun=float(arc(seg)[-1]),
                            top_r=float(np.linalg.norm(seg, axis=1).max()),
                            r_a=float(np.linalg.norm(a["point"])), r_b=float(np.linalg.norm(b["point"])),
                            theta_bn_a=float(dd["theta_bn"][0]), theta_bn_b=float(dd["theta_bn"][1]),
                            mach_f_a=float(dd["mach_f"][0]), mach_f_b=float(dd["mach_f"][1]),
                            mach_f_mc_a=float(dd["mach_f_mc"][0]), mach_f_mc_b=float(dd["mach_f_mc"][1]), n_lines=1))
            segs.append(seg)
    dist = None
    if radio_pts is not None and len(radio_pts):
        dist = np.full(len(radio_pts), np.inf)
        for t, seg in zip(out, segs):
            dd = cKDTree(seg).query(np.atleast_2d(radio_pts))[0]
            t["radio_min_distance_rsun"] = float(dd.min())
            t["n_radio_near"] = int((dd <= s["trap_near_rsun"]).sum())
            dist = np.minimum(dist, dd)
    return out, dist


# ------------------------------------------------------------------ 5. time evolution
def times_through(t_start, t_end, t0, h, v, s, extra=()):
    """Times from t_start - pad to t_end + pad every evolution_step_s seconds, and the times `extra`,
    at which the apex height h + v (t - t0) is at least evolution_min_height."""
    pad = timedelta(minutes=float(s["evolution_pad_min"]))
    a, b = t_start - pad, t_end + pad
    n = int((b - a).total_seconds() // float(s["evolution_step_s"])) + 1
    ts = sorted({a + timedelta(seconds=float(s["evolution_step_s"]) * i) for i in range(n)} | set(extra))
    return [t for t in ts if h + v * (t - t0).total_seconds() / R_SUN_KM >= s["evolution_min_height"]]


def _med(x, sel):
    x = np.asarray(x, float)[sel]
    x = x[np.isfinite(x)]
    return float(np.median(x)) if len(x) else np.nan


def evolution_row(S, axis, s):
    """Diagnostics of one surface for the evolution table."""
    mk = summary_masks(S, axis, s)
    a, adv = mk["apex"], mk["adv"]
    dg = S.dg
    row = dict(height_rsun=S.ht, apex_r=float(np.linalg.norm(S.xyz[a])), apex_mach_f=float(dg["mach_f"][a]),
               apex_mach_c=float(dg["mach_c"][a]), apex_mach_f_mc=float(dg["mach_f_mc"][a]), apex_theta_bn=float(dg["theta_bn"][a]),
               apex_standoff_rsun=float(S.delta[a]))
    for name in ("nose", "flank"):
        sel = mk[name] & adv
        row.update({f"{name}_mach_f": _med(dg["mach_f"], sel), f"{name}_mach_f_mc": _med(dg["mach_f_mc"], sel),
                    f"{name}_theta_bn": _med(dg["theta_bn"], sel)})
        row[f"{name}_mach_f_max"] = float(np.nanmax(dg["mach_f"][sel])) if sel.any() and np.isfinite(dg["mach_f"][sel]).any() else np.nan
    fr = fractions(S, np.ones(len(S.xyz), bool), s)
    row.update(frac_shock=fr["shock"], frac_supercritical=fr["supercritical"], frac_qperp_supercritical=fr["qperp_supercritical"])
    return row


def threshold_crossings(times, values, level=1.0):
    """Times at which a series crosses `level` (linear interpolation), with the direction up or down."""
    out = []
    v = np.asarray(values, float)
    for i in range(len(v) - 1):
        a, b = v[i], v[i + 1]
        if np.isfinite(a) and np.isfinite(b) and (a - level) * (b - level) < 0:
            f = (level - a) / (b - a)
            out.append(dict(time=times[i] + (times[i + 1] - times[i]) * f, direction="up" if b > a else "down"))
    return out


# ------------------------------------------------------------------ 6. uncertainty
def quantiles(x):
    """Median and 16th and 84th percentiles of the finite values."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if not len(x):
        return dict(n=0)
    q = np.percentile(x, [16, 50, 84])
    return dict(n=int(len(x)), median=float(q[1]), p16=float(q[0]), p84=float(q[2]))
