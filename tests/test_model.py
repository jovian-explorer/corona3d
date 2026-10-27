"""The MAS units, the Alfven speed, the rotation of the model, and the streamer selection, on a
synthetic run (uniform density, dipole field)."""
import numpy as np
import pytest

from conftest import write_mas
from corona3d import mast_model
from corona3d.reconstruction import reconstruct

M_P = 1.67262192e-24


def test_units_and_alfven_speed(tmp_path):
    d = write_mas(tmp_path / "run", n_model=2.0, b0=1.0)
    m = mast_model.load(d, density_scale=2.0)
    i, j, k = 5, 30, 0                                   # equator (theta = 90 deg), r = 1
    assert m["ne"][i, j, k] == pytest.approx(2.0 * 1e8 * 2.0, rel=1e-6)      # rho 1e8 cm^-3, x2
    assert m["bmag"][i, j, k] == pytest.approx(2.2068914, rel=1e-5)          # B unit [G]
    va = 2.2068914 / np.sqrt(4 * np.pi * 4e8 * M_P) / 1e5
    assert m["va"][i, j, k] == pytest.approx(va, rel=1e-5)                   # km/s, scaled density
    m2 = mast_model.load(d, density_scale=2.0, va_density="model")
    assert m2["va"][i, j, k] == pytest.approx(va * np.sqrt(2), rel=1e-5)     # model density


def test_missing_run_is_reported_with_advice(tmp_path):
    with pytest.raises(FileNotFoundError, match="corona3d fetch --model"):
        mast_model.load(tmp_path)


def test_model_rotation_turns_the_surface_eastward_or_westward(synthetic_home):
    """A density bump at Carrington longitude 0 appears at longitude model_dlon after rotation."""
    from corona3d.paths import config, model_dir
    p = config(synthetic_home)
    d = model_dir(synthetic_home, p)
    from pyhdf.SD import SD, SDC
    sd = SD(str(d / "rho002.hdf"), SDC.WRITE)
    s = sd.select("Data-Set-2")
    a = s[:]
    ph = sd.select("fakeDim0")[:]
    a[np.abs(np.angle(np.exp(1j * ph))) < 0.3, :, :] *= 4.0     # denser around longitude 0
    s[:] = a
    s.endaccess(); sd.end()
    on = dict(open=False, closed=False, streamer=False, gcs=False)
    for dlon in (0.0, 40.0):
        iso = reconstruct(dict(p, model_dlon=dlon, ne_level=2.5e8), on, synthetic_home)["surface"]
        pts = iso.points
        far = pts[np.linalg.norm(pts, axis=1) > np.percentile(np.linalg.norm(pts, axis=1), 99)]
        lon = np.degrees(np.angle(np.mean(far[:, 0] + 1j * far[:, 1])))
        assert lon == pytest.approx(dlon, abs=6)


def test_streamer_lines_are_the_closed_lines_threading_the_shell(synthetic_home):
    from corona3d.paths import config
    p = config(synthetic_home)
    sc = reconstruct(p, dict(surface=False, open=False, closed=False, gcs=True), synthetic_home)
    assert sc["streamer"]                                  # a dipole closes every low line
    shell = sc["gcs"].reshape(-1, 3)
    for line in sc["streamer"][:20]:
        r = np.linalg.norm(line, axis=1)
        assert r[0] < 1.0 + 0.05 and r[-1] < 1.0 + 0.05     # both feet on the inner boundary
        assert np.min(np.linalg.norm(shell[:, None] - line[None], axis=2)) < 0.05   # through the shell
