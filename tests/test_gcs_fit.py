"""The height-time step of scripts/gcs_fit.py: the limb view, the selection of the front, and the
extrapolation to the time of the reconstruction, on synthetic numbers with a known answer."""
import importlib.util
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pytest

SCRIPT = Path(__file__).resolve().parent.parent / "scripts" / "gcs_fit.py"
pytestmark = pytest.mark.skipif(not SCRIPT.exists(), reason="source tree not present")


@pytest.fixture(scope="module")
def G():
    pytest.importorskip("skimage")
    spec = importlib.util.spec_from_file_location("gcs_fit", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_projection_factor_is_one_at_the_limb_and_zero_along_the_line_of_sight(G):
    assert G.projection_factor(-53.0, 0.0, [np.cos(np.radians(37)), np.sin(np.radians(37)), 0]) == \
        pytest.approx(1.0, abs=1e-9)                          # 90 deg from the observer
    assert G.projection_factor(0.0, 0.0, [1, 0, 0]) == pytest.approx(0.0, abs=1e-9)


def test_front_profile_and_extrapolation_recover_a_known_rise(G):
    """Front at 1.9 R_sun at 15:00 rising at 0.09 R_sun/min; noise frames before the eruption
    (small envelopes) and truncated fronts at the edge of a 4 R_sun field must be left out."""
    t0 = datetime(2013, 10, 25, 15, 0)
    rows = [(t0 - timedelta(minutes=10), 2.4, 0.3),          # pre-event noise at the seed
            (t0 - timedelta(minutes=2), 2.1, 0.6)]
    for m in (5, 8, 13, 18):
        rows.append((t0 + timedelta(minutes=m), 1.9 + 0.09 * m, 2.0 + m / 4))
    rows += [(t0 + timedelta(minutes=23), 3.95, 9.0),           # front at the edge: truncated
             (t0 + timedelta(minutes=28), 3.70, 9.5)]           # what remains after it left
    kept = G.front_profile(rows, r_out=4.0)
    assert [r for _, r, _ in kept] == pytest.approx([1.9 + 0.09 * m for m in (5, 8, 13, 18)])
    v, h0, err = G.extrapolate([t for t, _, _ in kept], [r for _, r, _ in kept], t0)
    assert v == pytest.approx(0.09, abs=1e-9) and h0 == pytest.approx(1.9, abs=1e-9)
    assert err == pytest.approx(0.0, abs=1e-6)
