"""Field-line tracing on the MAS grid, the GCS shell, and observer frames."""
import warnings

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial.transform import Rotation


def unit_vec(lon_deg, lat_deg):
    lo, la = np.radians(lon_deg), np.radians(lat_deg)
    return np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])


def observer_basis(lon_deg, lat_deg):
    """Heliocentric-Cartesian axes (X west, Y north, Z to observer) in the Carrington frame."""
    ez = unit_vec(lon_deg, lat_deg)
    north = np.array([0.0, 0.0, 1.0])
    ey = north - ez * (north @ ez)
    ey /= np.linalg.norm(ey)
    ex = np.cross(ey, ez)
    return ex, ey, ez


def deproject(x_pos, y_pos, r, obs_lon, obs_lat):
    """Plane-of-sky (x west, y north, R_sun) at radial distance r -> Carrington Cartesian,
    taking the solution on the observer side of the plane of sky. A position farther from Sun
    centre than r cannot lie at r; it is kept in the plane of the sky, with a warning."""
    ex, ey, ez = observer_basis(obs_lon, obs_lat)
    rho2 = np.asarray(x_pos) ** 2 + np.asarray(y_pos) ** 2
    if (rho2 > r**2).any():
        warnings.warn(f"{int((rho2 > r**2).sum())} radio position(s) lie farther than {r:.3f} R_sun "
                      "from Sun centre in the plane of the sky; kept in the plane of the sky")
    z = np.sqrt(np.clip(r**2 - rho2, 0, None))
    return np.outer(x_pos, ex) + np.outer(y_pos, ey) + np.outer(z, ez)


# ---------------------------------------------------------------- field lines
OPEN, CLOSED, UNDETERMINED = "open", "closed", "undetermined"


class FieldTracer:
    """RK4 integration along the unit vector of B, with a step proportional to r."""

    def __init__(self, m):
        self.r_min, self.r_max = m["r"][0], m["r"][-1]
        axes = (m["phi"], m["theta"], m["r"])
        self.f = [RegularGridInterpolator(axes, m[k], bounds_error=False, fill_value=0.0)
                  for k in ("br", "bt", "bp")]

    def direction(self, xyz):
        x, y, z = xyz.T
        r = np.linalg.norm(xyz, axis=1)
        th = np.arccos(np.clip(z / r, -1, 1))
        ph = np.mod(np.arctan2(y, x), 2 * np.pi)
        q = np.c_[ph, th, np.clip(r, self.r_min, self.r_max)]
        br, bt, bp = (f(q) for f in self.f)
        st, ct, sp, cp = np.sin(th), np.cos(th), np.sin(ph), np.cos(ph)
        bx = br * st * cp + bt * ct * cp - bp * sp
        by = br * st * sp + bt * ct * sp + bp * cp
        bz = br * ct - bt * st
        b = np.c_[bx, by, bz]
        return b / np.maximum(np.linalg.norm(b, axis=1, keepdims=True), 1e-30)

    def trace(self, seeds, ds=0.01, n_steps=1500, sign=1.0):
        """Trace all seeds in one direction with steps of ds * r [R_sun]; returns a list of
        (n, 3) arrays, each ending at the first point outside [r_min, r_max] or after n_steps."""
        p = np.array(seeds, float)
        alive = np.ones(len(p), bool)
        P = np.empty((len(p), n_steps + 1, 3), np.float32)
        P[:, 0] = p
        last = np.full(len(p), n_steps)
        for k in range(n_steps):
            idx = np.nonzero(alive)[0]
            if idx.size == 0:
                last[:] = np.minimum(last, k)
                break
            q = p[idx]
            h = (ds * sign * np.linalg.norm(q, axis=1))[:, None]
            k1 = self.direction(q)
            k2 = self.direction(q + 0.5 * h * k1)
            k3 = self.direction(q + 0.5 * h * k2)
            k4 = self.direction(q + h * k3)
            p[idx] = q + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
            P[:, k + 1] = p
            r = np.linalg.norm(p[idx], axis=1)
            out = idx[(r < self.r_min) | (r > self.r_max)]
            alive[out] = False
            last[out] = k + 1
        return [P[i, : last[i] + 1].astype(float) for i in range(len(p))]

    def field_lines(self, seeds, **kw):
        """Trace both ways; return a list of (points, kind): kind is OPEN when one end leaves
        the outer boundary, CLOSED when both ends return to the inner boundary, and
        UNDETERMINED otherwise (a line that ran out of steps or stalled where B vanishes)."""
        fw = self.trace(seeds, sign=+1.0, **kw)
        bw = self.trace(seeds, sign=-1.0, **kw)
        lines = []
        for a, b in zip(fw, bw):
            pts = np.vstack([b[::-1], a[1:]])
            r0, r1 = np.linalg.norm(pts[[0, -1]], axis=1)
            if r0 > self.r_max or r1 > self.r_max:
                kind = OPEN
            elif r0 < self.r_min and r1 < self.r_min:
                kind = CLOSED
            else:
                kind = UNDETERMINED
            lines.append((pts, kind))
        return lines


def fibonacci_sphere(n, r=1.0):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5**0.5) * i
    return r * np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]


# ---------------------------------------------------------------- GCS shell
def pythea_shell(height, kappa, alpha_deg, lon_stonyhurst, lat_deg, tilt_pythea_deg,
                 earth_carr_lon, obstime, **n):
    """GCS shell from PyThea (Thernisien skeleton, PyThea conventions), returned as
    rings (n_axis, n_circle, 3) in Carrington-Cartesian R_sun."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from sunpy.coordinates import frames
    from PyThea.geometrical_models import gcs

    rc = gcs.rcenter_(height * u.R_sun, alpha_deg * u.deg, kappa)
    c = SkyCoord(lon_stonyhurst * u.deg, lat_deg * u.deg, rc,
                 frame=frames.HeliographicStonyhurst, obstime=obstime)
    m = gcs(c, height * u.R_sun, alpha_deg * u.deg, kappa, tilt_pythea_deg * u.deg, **n)
    x, y, z = m.cloud()
    X, Y, Z = m.rotate(x * u.R_sun, y * u.R_sun, z * u.R_sun)
    xyz = np.c_[X.value, Y.value, Z.value] @ Rotation.from_euler(
        "z", earth_carr_lon, degrees=True).as_matrix().T          # Stonyhurst -> Carrington
    n_axis = 2 * (m.nbvertcirc + m.nbvertsl) - 3
    return xyz.reshape(n_axis, m.nbvertcircshell, 3)
