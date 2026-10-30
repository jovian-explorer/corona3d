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
from pathlib import Path

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
    f = lambda pts: pl.at(pts)["ne"]
    f.r_max = float(m["r"][-1])                             # outer radius of the model, for the radial profiles
    return f


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


# ------------------------------------------------------------------ emission (F or H) and placement of the sources
def emission_of(p, entry=None, freq=None):
    """"F", "H", or "auto" for a radio source: the key `emission` of its radio_images entry (a letter, or a mapping
    from the frequency in MHz, written as in the entry, to a letter, with an optional "default"), else
    `radio_emission` of the configuration (default "H", the convention of the radio layer)."""
    base = str(p.get("radio_emission", "H"))
    e = (entry or {}).get("emission", base)
    if isinstance(e, dict):
        e = next((v for k, v in e.items() if k != "default" and freq is not None and abs(float(k) - freq) < 0.05),
                 e.get("default", base))
    e = str(e).strip()
    if e.lower() == "auto":
        return "auto"
    if e.upper() not in ("F", "H"):
        raise ValueError(f"emission {e!r}: give F (fundamental), H (harmonic), or auto")
    return e.upper()


def emission_from_polarisation(dcp, valid, threshold):
    """("F" or "H", reason) from the degree of circular polarisation V/I of a source: fundamental emission is
    usually strongly polarised and harmonic emission weakly, so that |V/I| at or above `threshold` is taken as F."""
    if dcp is None or not np.isfinite(dcp) or not valid:
        return "H", "no valid V/I; harmonic assumed"
    if abs(dcp) >= threshold:
        return "F", f"|V/I| = {abs(dcp):.2f} >= {threshold:g}"
    return "H", f"|V/I| = {abs(dcp):.2f} < {threshold:g}"


def custom(p, entry=None):
    """True when a source is placed otherwise than on the default sphere of the four-fold Newkirk height for
    harmonic emission (another surface, F emission, the scattering correction, or a position offset)."""
    return (p.get("radio_deproject", "sphere") != "sphere" or emission_of(p, entry) != "H"
            or (entry or {}).get("emission") not in (None, "H", "h")
            or bool((p.get("scattering") or {}).get("apply")) or bool((entry or {}).get("position_offset")))


def locate(x, y, f, p, harmonic=True, mode="sphere", ne_at=None, r_max_model=30.0, scatter=False, eps2_h=None):
    """Three-dimensional position of one source at the sky-plane position (x, y) [R_sun, seen from Earth].

    mode "sphere": on the sphere of the four-fold Newkirk height of f; "model": where the line of sight first meets
    the plasma level of f in the scaled MAS density (falling back to the sphere when it does not reach it). With
    `scatter` the observed position is first moved inward by the radial shift of scattering (corona3d.scattering),
    computed in the density of the chosen surface (four-fold Newkirk, or the MAS density along the radius through
    the source). Returns dict(r, xyz, surface, dr, dr_lo, dr_hi, x, y) with x, y the position used."""
    from . import scattering
    from .geometry import deproject

    def place(xx, yy):
        if mode == "model" and ne_at is not None:
            r = mas_los_r(ne_at, xx, yy, f, p["earth"], harmonic)
            if np.isfinite(r):
                return r, deproject([xx], [yy], r, *p["earth"], warn=False)[0], "model"
        r = float(newkirk_r(f, 4.0, harmonic))
        return r, deproject([xx], [yy], r, *p["earth"], warn=False)[0], "sphere"
    r, xyz, surf = place(x, y)
    out = dict(r=r, xyz=xyz, surface=surf, dr=0.0, dr_lo=0.0, dr_hi=0.0, x=x, y=y)
    if scatter:
        prof = (scattering.radial_profile(ne_at, xyz, r_max_model) if surf == "model"
                else scattering.newkirk_profile(4.0))
        dr, lo, hi = scattering.shift(f, r, prof, harmonic, tuple(eps2_h or scattering.EPS2_H))
        xc, yc = scattering.correct(x, y, r, dr)
        r2, xyz2, surf2 = place(float(xc), float(yc))
        out.update(r=r2, xyz=xyz2, surface=surf2, dr=dr, dr_lo=lo, dr_hi=hi, x=float(xc), y=float(yc))
    return out


def place_blobs(entry, f, xy, times, p, data_root, ne_loader):
    """Positions of the centroids of one frequency of a radio_images entry when `custom` is true: the position
    offset of the entry (OVRO-LWA, corona3d.ovro_lwa.offset_for) is subtracted, the emission is F or H (or chosen
    from the NRH polarisation of the frame with "auto"), and every centroid is placed with `locate`.
    ne_loader() returns the sampler of the scaled MAS density. Returns (xyz (n, 3), r (n,), emission list)."""
    from . import ovro_lwa
    mode = p.get("radio_deproject", "sphere")
    if mode not in ("sphere", "model"):
        raise ValueError(f'radio_deproject {mode!r}: give "sphere" or "model"')
    sc = p.get("scattering") or {}
    dx, dy = ovro_lwa.offset_for(data_root, entry, f) if entry.get("instrument") == "ovro_lwa" else (0.0, 0.0)
    em0 = emission_of(p, entry, f)
    ne_at = ne_loader() if (mode == "model" or sc.get("apply")) else None
    r_max = float(getattr(ne_at, "r_max", 30.0))         # the model's outer radius, as in sources.analyse
    xyz, rr, ems = [], [], []
    for (x, y), t in zip(xy, times):
        em = em0
        if em0 == "auto":
            em = "H"
            if entry.get("instrument") == "nrh":
                from . import nrh
                folder = Path(data_root) / entry["folder"]
                img, _, _, _ = nrh.image_at(folder, f, t)
                v, _ = nrh.stokes_v_at(folder, f, t)
                if v is not None:
                    src, pk = nrh.polarisation(img, v, entry["level"])
                    em, _ = emission_from_polarisation(src, abs(src) <= 1 and abs(pk) <= 1, p.get("dcp_fundamental", 0.3))
        q = locate(x - dx, y - dy, f, p, em == "H", mode, ne_at, r_max, bool(sc.get("apply")), sc.get("eps2_h_km"))
        xyz.append(q["xyz"])
        rr.append(q["r"])
        ems.append(em)
    return np.array(xyz), np.array(rr), ems
