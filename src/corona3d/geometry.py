"""Field-line tracing on the MAS grid, the GCS shell, and observer frames."""
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
    taking the solution on the observer side of the plane of sky."""
    ex, ey, ez = observer_basis(obs_lon, obs_lat)
    z = np.sqrt(np.clip(r**2 - x_pos**2 - y_pos**2, 0, None))
    return np.outer(x_pos, ex) + np.outer(y_pos, ey) + np.outer(z, ez)


# ---------------------------------------------------------------- field lines
class FieldTracer:
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

    def trace(self, seeds, ds=0.01, n_steps=800, sign=1.0):
        """RK4 trace of all seeds in one direction; returns list of (n,3) arrays."""
        p = seeds.copy()
        alive = np.ones(len(p), bool)
        paths = [p.copy()]
        for _ in range(n_steps):
            idx = np.nonzero(alive)[0]
            if idx.size == 0:
                break
            q = p[idx]
            k1 = self.direction(q)
            k2 = self.direction(q + 0.5 * ds * sign * k1)
            k3 = self.direction(q + 0.5 * ds * sign * k2)
            k4 = self.direction(q + ds * sign * k3)
            q = q + ds * sign * (k1 + 2 * k2 + 2 * k3 + k4) / 6
            p[idx] = q
            r = np.linalg.norm(q, axis=1)
            alive[idx[(r < self.r_min) | (r > self.r_max)]] = False
            paths.append(p.copy())
        P = np.stack(paths, axis=1)                    # (nseed, nstep, 3)
        out = []
        for i in range(len(seeds)):
            r = np.linalg.norm(P[i], axis=1)
            # keep until the first exit, then drop the frozen repeats
            stop = np.argmax((r < self.r_min) | (r > self.r_max)) or len(r) - 1
            out.append(P[i, : stop + 1])
        return out

    def field_lines(self, seeds, **kw):
        """Trace both ways; return list of (points, is_open)."""
        fw = self.trace(seeds, sign=+1.0, **kw)
        bw = self.trace(seeds, sign=-1.0, **kw)
        lines = []
        for a, b in zip(fw, bw):
            pts = np.vstack([b[::-1], a[1:]])
            r_end = np.linalg.norm(pts[[0, -1]], axis=1)
            is_open = bool((r_end > self.r_max - 0.05).any())
            lines.append((pts, is_open))
        return lines


def fibonacci_sphere(n, r=1.0):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5**0.5) * i
    return r * np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]


# ---------------------------------------------------------------- GCS shell
def pythea_shell(height, kappa, alpha_deg, lon_stonyhurst, lat_deg, tilt_pythea_deg,
                 earth_carr_lon, obstime="2013-10-25 15:00", **n):
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
