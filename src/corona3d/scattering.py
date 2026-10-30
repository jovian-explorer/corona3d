"""Apparent displacement and broadening of radio sources by the scattering of radio waves on density
turbulence (Chrysaphi et al. 2018; Kontar et al. 2019).

Radio waves emitted near the plasma frequency are scattered by small-scale density fluctuations, and the
observer sees them arrive from the distance at which the scattering becomes weak. Following Chrysaphi et al.
(2018, their Eqs. 8 and 9), the angular scattering of a wave of frequency f in fluctuations
<dn^2> = eps^2 n^2 with correlation length h accumulates as

    d<dtheta^2>/dr = (sqrt(pi) / 2) f_pe(r)^4 / (f^2 - f_pe(r)^2)^2 (eps^2 / h),

and the optical depth to scattering from r outward is tau(r) = Int_r^{1 AU} d<dtheta^2>/dr dr. The apparent
source lies where tau = 1, outward of the true source at the plasma level, so that the radial shift is

    dr = r(tau = 1) - r_true.

With eps^2/h from 4.5e-5 to 7e-5 km^-1 (the range of Chrysaphi et al. 2018) the shift has a range, which is its
uncertainty here. Since the shift is radial, a sky-plane position (x, y) observed at the apparent distance
r_true + dr belongs to a true position (x, y) r_true / (r_true + dr), whatever its angle from the line of sight
(the projected shift is proportional to sin theta_s, Kontar et al. 2019). The emission frequency is f = f_pe at the
true source for fundamental and f = 2 f_pe for harmonic emission; harmonic waves are scattered much less.

Kontar et al. (2019) fitted the observed FWHM sizes of type III sources from 0.05 to 500 MHz with

    FWHM = 11.78 f^-0.98 deg   (f in MHz),

which is given beside every source as the size that scattering alone produces.
"""
import numpy as np

R_SUN_KM = 695700.0
FP = 8980.0
EPS2_H = (4.5e-5, 7.0e-5)          # km^-1, Chrysaphi et al. (2018)
AU_RSUN = 215.0


def fwhm_deg(f_mhz):
    """FWHM size [deg] of a scattered source at f [MHz] (Kontar et al. 2019, Eq. 50)."""
    return 11.78 * np.asarray(f_mhz, float) ** -0.98


def optical_depth(r, f_hz, ne_of_r, eps2_h, r_out=AU_RSUN, n=4000):
    """tau(r) for the radial density profile ne_of_r(r) [cm^-3] (vectorised in r)."""
    r = np.atleast_1d(np.asarray(r, float))
    out = np.empty(len(r))
    for i, r0 in enumerate(r):
        rr = np.geomspace(r0, r_out, n)
        fpe = FP * np.sqrt(ne_of_r(rr))
        with np.errstate(divide="ignore", invalid="ignore"):
            g = np.where(fpe < f_hz, np.sqrt(np.pi) / 2 * fpe**4 / (f_hz**2 - fpe**2) ** 2, np.inf)
        out[i] = np.trapezoid(g, rr) * eps2_h * R_SUN_KM
    return out


def _tau_unit(f_hz, r_true, ne_of_r, harmonic, n=20000):
    """Grid r from just above r_true to 1 AU and tau(r) per unit eps^2/h [km], integrated from outside."""
    r0 = float(r_true) * (1 + (1e-6 if harmonic else 1e-5))
    rr = np.geomspace(r0, AU_RSUN, n)
    fpe = FP * np.sqrt(ne_of_r(rr))
    with np.errstate(divide="ignore", invalid="ignore"):
        g = np.where(fpe < f_hz, np.sqrt(np.pi) / 2 * fpe**4 / (f_hz**2 - fpe**2) ** 2, np.nan)
    g = np.nan_to_num(g, nan=np.nanmax(g) if np.isfinite(g).any() else 0.0)
    seg = 0.5 * (g[1:] + g[:-1]) * np.diff(rr) * R_SUN_KM
    return rr, np.r_[np.cumsum(seg[::-1])[::-1], 0.0]


def apparent_r(f_mhz, r_true, ne_of_r, eps2_h, harmonic=True):
    """Distance [R_sun] where tau = 1 for a source at r_true emitting at f (f = 2 f_pe for harmonic); r_true when
    the scattering is weak already there."""
    rr, tu = _tau_unit(float(f_mhz) * 1e6, r_true, ne_of_r, harmonic)
    tau = tu * eps2_h
    above = np.flatnonzero(tau >= 1.0)
    if not len(above):
        return float(r_true)
    i = above[-1]
    if i + 1 >= len(rr):
        return float(rr[-1])
    a, b = np.log(tau[i]), np.log(max(tau[i + 1], 1e-300))       # interpolation in log tau
    return float(rr[i] + (rr[i + 1] - rr[i]) * (a / (a - b) if a != b else 0.0))


def shift(f_mhz, r_true, ne_of_r, harmonic=True, eps2_h=EPS2_H):
    """Radial shift [R_sun] by scattering: (value at the mean eps^2/h, low, high) over the range eps2_h."""
    rr, tu = _tau_unit(float(f_mhz) * 1e6, r_true, ne_of_r, harmonic)

    def one(e):
        tau = tu * e
        above = np.flatnonzero(tau >= 1.0)
        if not len(above):
            return 0.0
        i = min(above[-1], len(rr) - 2)
        a, b = np.log(tau[i]), np.log(max(tau[i + 1], 1e-300))
        return float(rr[i] + (rr[i + 1] - rr[i]) * (a / (a - b) if a != b else 0.0)) - float(r_true)
    return one(float(np.mean(eps2_h))), one(min(eps2_h)), one(max(eps2_h))


def correct(x, y, r_true, dr):
    """True sky-plane position of a source observed at (x, y) whose distance r_true was shifted outward by dr."""
    k = r_true / (r_true + dr)
    return np.asarray(x) * k, np.asarray(y) * k


R_JOIN = 4.0                       # R_sun: beyond this the Newkirk density is continued as r^-2


def newkirk_profile(fold, r_join=R_JOIN):
    """fold-times Newkirk (1961) density [cm^-3] up to r_join, continued as r^-2 beyond it (the Newkirk model tends
    to a constant far from the Sun, which would make the integral of the scattering to 1 AU grow without limit)."""
    nj = fold * 4.2e4 * 10 ** (4.32 / r_join)

    def f(r):
        r = np.asarray(r, float)
        return np.where(r <= r_join, fold * 4.2e4 * 10 ** (4.32 / np.minimum(r, r_join)), nj * (r_join / np.maximum(r, r_join)) ** 2)
    return f


def radial_profile(ne_at, direction, r_max_model):
    """ne(r) along the radial line through a unit vector, from a three-dimensional sampler; beyond the model
    the density falls as r^-2 from its outer boundary."""
    d = np.asarray(direction, float) / np.linalg.norm(direction)

    def f(r):
        r = np.asarray(r, float)
        inside = np.minimum(r, r_max_model)
        ne = ne_at(np.outer(inside, d))
        return np.where(r > r_max_model, ne * (r_max_model / r) ** 2, ne)
    return f
