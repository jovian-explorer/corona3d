"""Read PSI MAS/MAST coronal output and derive electron density, |B| and Alfven speed.

MAS code units (Predictive Science Inc.):
    length  1 R_sun
    rho     1e8 cm^-3  (number density n, proton mass per particle)
    B       2.2068914 G
    v       B / sqrt(4 pi n m_p) = 481.37 km/s
    T       m_p v^2 / k_B = B^2 / (4 pi n k_B) = 2.807e7 K   (the run has p = 2 rho T exactly:
            electrons and protons at one temperature)
"""
from pathlib import Path

import numpy as np
from pyhdf.SD import SD
from scipy.interpolate import RegularGridInterpolator

B_UNIT = 2.2068914          # G
N_UNIT = 1.0e8              # cm^-3
M_P = 1.67262192e-24        # g
K_B = 1.380649e-16          # erg/K
V_UNIT = B_UNIT / np.sqrt(4 * np.pi * N_UNIT * M_P) / 1e5      # km/s (481.37)
T_UNIT = B_UNIT**2 / (4 * np.pi * N_UNIT * K_B)                # K (2.807e7)


def read_mas(name, data_dir):
    """Return (phi, theta, r, array[phi, theta, r]) for one MAS HDF4 file."""
    sd = SD(str(Path(data_dir) / f"{name}.hdf"))
    ds = sd.datasets()
    dims = sorted(k for k in ds if k.startswith("fakeDim"))
    phi, theta, r = (sd.select(k)[:].astype(float) for k in dims)
    arr = sd.select("Data-Set-2")[:].astype(float)
    return phi, theta, r, arr


def to_grid(name, phi, theta, r, data_dir):
    """Interpolate a (possibly staggered) MAS variable onto the given grid."""
    p0, t0, r0, a = read_mas(name, data_dir)
    f = RegularGridInterpolator((p0, t0, r0), a, bounds_error=False, fill_value=None)
    P, T, R = np.meshgrid(phi, theta, r, indexing="ij")
    return f(np.stack([P, T, R], axis=-1))


FILES = ("rho002", "br002", "bt002", "bp002", "va002")


def load(data_dir, density_scale=2.0, r_max=None, va_density="scaled"):
    """Common grid (phi, theta, r) with n_e [cm^-3], B vector [G], v_A [km/s].

    The grid is the cell-centre mesh of va002 (only its axes are used), up to r_max [R_sun]
    (default: the whole run). density_scale multiplies the model density (Morosan et al. 2022
    scale the low-coronal density by 2 after Wang et al. 2017). v_A is computed from the scaled
    density when va_density is "scaled", and from the model density when it is "model".
    """
    missing = [f for f in FILES if not (Path(data_dir) / f"{f}.hdf").exists()]
    if missing:
        raise FileNotFoundError(f"{', '.join(missing)}.hdf missing in {data_dir}; run `corona3d fetch --model`")
    if va_density not in ("scaled", "model"):
        raise ValueError('va_density must be "scaled" or "model"')
    phi, theta, r, _ = read_mas("va002", data_dir)
    if r_max:
        r = r[r <= r_max]
    n_model = to_grid("rho002", phi, theta, r, data_dir) * N_UNIT
    ne = n_model * density_scale
    br = to_grid("br002", phi, theta, r, data_dir) * B_UNIT
    bt = to_grid("bt002", phi, theta, r, data_dir) * B_UNIT
    bp = to_grid("bp002", phi, theta, r, data_dir) * B_UNIT
    bmag = np.sqrt(br**2 + bt**2 + bp**2)
    va = bmag / np.sqrt(4 * np.pi * (ne if va_density == "scaled" else n_model) * M_P) / 1e5   # km/s
    return dict(phi=phi, theta=theta, r=r, ne=ne, br=br, bt=bt, bp=bp, bmag=bmag, va=va)


EXTRA = ("t002", "vr002", "vt002", "vp002")          # temperature and flow, used by the shock diagnostics


def load_plasma(data_dir, m):
    """Add the temperature T [K] and the flow vr, vt, vp [km/s] of the run to the grid m of load(),
    for each file of EXTRA that is present; returns the names that were found."""
    found = []
    for name, key, unit in (("t002", "T", T_UNIT), ("vr002", "vr", V_UNIT), ("vt002", "vt", V_UNIT),
                            ("vp002", "vp", V_UNIT)):
        if (Path(data_dir) / f"{name}.hdf").exists():
            m[key] = to_grid(name, m["phi"], m["theta"], m["r"], data_dir) * unit
            found.append(name)
    return found
