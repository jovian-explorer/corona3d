"""Heights of the radio sources in several coronal density models, and their range.

A radio source at the frequency f is placed where the plasma frequency of the corona equals f/2
(second-harmonic emission, the convention of the radio layer), that is, at the electron density

    n_e(f) = (f / (2 x 8980 Hz))^2  cm^-3.

The height of that density differs between density models, and the spread over the models is the
uncertainty of the height that the observation alone cannot remove. The models are

    newkirkK   K-fold Newkirk (1961):  n_e(r) = K 4.2e4 10^(4.32 / r)
    saito      Saito et al. (1977), equatorial background corona:
               n_e(r) = 1.36e6 r^-2.14 + 1.68e8 r^-6.13
    mas        the scaled density of the MAS/MAST run of the event along the line of sight of the
               source, at the first point from the observer at which it reaches n_e(f)

with r in R_sun and n_e in cm^-3. The range of a source is the smallest and largest of these
heights, and its three-dimensional end points lie on the line of sight of the observed position.
"""
import numpy as np
from scipy.optimize import brentq

FP = 8980.0                                   # Hz per sqrt(cm^-3)


def level(f_mhz, harmonic=True):
    """Electron density [cm^-3] whose plasma frequency is f (or f/2 for harmonic emission)."""
    return (np.asarray(f_mhz, float) * 1e6 / (2 if harmonic else 1) / FP) ** 2


def newkirk_r(f_mhz, fold=4.0, harmonic=True):
    """Height [R_sun] of the plasma level of f in the fold-times Newkirk (1961) model."""
    return 4.32 / np.log10(level(f_mhz, harmonic) / (fold * 4.2e4))


def saito_ne(r):
    """Electron density [cm^-3] of the equatorial background corona of Saito et al. (1977)."""
    r = np.asarray(r, float)
    return 1.36e6 * r**-2.14 + 1.68e8 * r**-6.13


def saito_r(f_mhz, harmonic=True):
    """Height [R_sun] of the plasma level of f in the model of Saito et al. (1977); NaN when the level
    is denser than the model at 1 R_sun."""
    out = []
    for n in np.atleast_1d(level(f_mhz, harmonic)):
        out.append(np.nan if saito_ne(1.0) < n else brentq(lambda r: saito_ne(r) - n, 1.0, 215.0))
    return np.array(out) if np.ndim(f_mhz) else out[0]


def mas_los_r(ne_at, x, y, f_mhz, observer, harmonic=True, r_out=10.0, step=0.005):
    """Height [R_sun] along the line of sight of the sky-plane position (x, y) [R_sun] of an observer
    (Carrington longitude, latitude) at which the density ne_at(points) first reaches the plasma level
    of f, coming from the observer, on the observer side of the plane of the sky (as the deprojection of
    the radio layer); NaN when it does not reach it there."""
    from .geometry import observer_basis
    ex, ey, ez = observer_basis(*observer)
    rho2 = x * x + y * y
    if rho2 >= r_out**2:
        return np.nan
    z_top = np.sqrt(r_out**2 - rho2)
    z_bot = np.sqrt(1.0 - rho2) if rho2 < 1.0 else 0.0
    z = np.arange(z_top, z_bot - step / 2, -step)
    if len(z) < 2:
        return np.nan
    pts = np.outer(np.full(len(z), x), ex) + np.outer(np.full(len(z), y), ey) + np.outer(z, ez)
    ne = ne_at(pts)
    n = float(level(f_mhz, harmonic))
    hit = np.flatnonzero(ne >= n)
    if not len(hit):
        return np.nan
    i = hit[0]
    if i == 0:
        zz = z[0]
    else:                                         # interpolation in log n_e between the two samples
        a, b = np.log(ne[i - 1]), np.log(ne[i])
        zz = z[i - 1] + (z[i] - z[i - 1]) * (np.log(n) - a) / (b - a) if b != a else z[i]
    return float(np.sqrt(rho2 + zz * zz))


def model_r(name, f_mhz, x=None, y=None, ne_at=None, observer=None):
    """Height [R_sun] of one source in the model `name` (newkirkK, saito, or mas)."""
    if name.startswith("newkirk"):
        return float(newkirk_r(f_mhz, float(name[7:] or 4)))
    if name == "saito":
        return float(saito_r(f_mhz))
    if name == "mas":
        return mas_los_r(ne_at, x, y, f_mhz, observer) if ne_at is not None else np.nan
    raise ValueError(f"height model {name!r}: use newkirkK (e.g. newkirk4), saito, or mas")


def ne_sampler(m, p):
    """ne_at(points) for the scaled density of the run, with the rotation model_dlon."""
    from .shock import Plasma
    pl = Plasma(m, p, temperature_k=0.0, flow=False)
    return lambda pts: pl.at(pts)["ne"]


def ranges(rt, p, m=None, models=("newkirk1", "newkirk2", "newkirk4", "newkirk6", "saito", "mas")):
    """Heights of every row of a radio table (corona3d.api RADIO_COLUMNS) in each model, their range,
    and the three-dimensional end points of the range on the line of sight from Earth. Returns columns
    r_<model>, r_min, r_max, x_min, y_min, z_min, x_max, y_max, z_max."""
    from .geometry import deproject
    ne_at = ne_sampler(m, p) if (m is not None and "mas" in models) else None
    out = {f"r_{k}": [] for k in models}
    for f, x, y in zip(rt["freq_mhz"], rt["x_rsun"], rt["y_rsun"]):
        for k in models:
            out[f"r_{k}"].append(model_r(k, f, x, y, ne_at, p["earth"]))
    out = {k: np.array(v, float) for k, v in out.items()}
    allr = np.vstack([out[f"r_{k}"] for k in models]) if models else np.empty((0, len(rt["x_rsun"])))
    with np.errstate(all="ignore"):
        out["r_min"], out["r_max"] = np.nanmin(allr, axis=0), np.nanmax(allr, axis=0)
    for end in ("min", "max"):
        q = np.array([deproject([x], [y], r, *p["earth"], warn=False)[0] if np.isfinite(r) else [np.nan] * 3
                      for x, y, r in zip(rt["x_rsun"], rt["y_rsun"], out[f"r_{end}"])]).reshape(-1, 3)
        out[f"x_{end}"], out[f"y_{end}"], out[f"z_{end}"] = q.T
    return out
