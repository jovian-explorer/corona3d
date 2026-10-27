"""A synthetic MAS run (uniform density, dipole field) written as HDF4 files, for offline tests."""
import numpy as np
import pytest


def write_mas(d, n_model=1.0, b0=1.0):
    """rho = n_model (code units), B = dipole of strength b0 (code units) on a small grid."""
    from pyhdf.SD import SD, SDC
    phi = np.linspace(0, 2 * np.pi, 33)
    th = np.linspace(0, np.pi, 61)
    r = np.linspace(1.0, 4.0, 61)
    P, T, R = np.meshgrid(phi, th, r, indexing="ij")
    fields = dict(rho002=np.full(P.shape, n_model), br002=b0 * 2 * np.cos(T) / R**3,
                  bt002=b0 * np.sin(T) / R**3, bp002=np.zeros(P.shape), va002=np.zeros(P.shape))
    d.mkdir(parents=True, exist_ok=True)
    for name, arr in fields.items():
        sd = SD(str(d / f"{name}.hdf"), SDC.WRITE | SDC.CREATE)
        for i, ax in enumerate((phi, th, r)):
            s = sd.create(f"fakeDim{i}", SDC.FLOAT32, (len(ax),))
            s[:] = ax.astype(np.float32)
            s.endaccess()
        s = sd.create("Data-Set-2", SDC.FLOAT32, arr.shape)
        s[:] = arr.astype(np.float32)
        s.endaccess()
        sd.end()
    return d


@pytest.fixture
def synthetic_home(tmp_path):
    """A working directory with config.json and the synthetic run under its default name."""
    from corona3d.cli import main
    from corona3d.paths import config, model_dir
    main(["--home", str(tmp_path), "init"])
    write_mas(model_dir(tmp_path, config(tmp_path)))
    return tmp_path
