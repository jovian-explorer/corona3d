"""Observer frames, deprojection, the radio heights, and the field-line tracer on a dipole."""
import numpy as np
import pytest

from corona3d import newkirk_r
from corona3d.geometry import FieldTracer, deproject, fibonacci_sphere, observer_basis, unit_vec


def test_newkirk_heights_of_the_paper():
    # 4x Newkirk at the harmonic: 1.576 R_sun at 173 MHz and 1.449 R_sun at 228 MHz
    assert newkirk_r(173) == pytest.approx(1.576, abs=1e-3)
    assert newkirk_r(228) == pytest.approx(1.449, abs=1e-3)
    assert newkirk_r(150) > newkirk_r(173) > newkirk_r(228) > newkirk_r(327)


@pytest.mark.parametrize("lon,lat", [(0, 0), (357.372, 4.993), (215.329, -0.984), (90, -7)])
def test_observer_basis_is_right_handed_with_x_west(lon, lat):
    ex, ey, ez = observer_basis(lon, lat)
    B = np.array([ex, ey, ez])
    assert np.allclose(B @ B.T, np.eye(3))
    assert np.allclose(np.cross(ex, ey), ez)
    assert ey[2] > 0                                   # solar north up
    # west limb: Stonyhurst longitude +90 deg = Carrington longitude of the observer + 90 deg
    assert ex @ unit_vec(lon + 90, 0) == pytest.approx(np.cos(np.radians(0)) * 1, abs=0.01) or \
        ex @ unit_vec(lon + 90, 0) > 0.99


def test_deproject_puts_sources_on_the_sphere_facing_the_observer():
    x, y = np.array([-1.1, 0.2]), np.array([-0.1, 0.5])
    p = deproject(x, y, 1.576, 357.372, 4.993)
    ex, ey, ez = observer_basis(357.372, 4.993)
    assert np.allclose(np.linalg.norm(p, axis=1), 1.576)
    assert np.allclose(p @ ex, x) and np.allclose(p @ ey, y)
    assert (p @ ez > 0).all()


def test_fibonacci_sphere():
    s = fibonacci_sphere(500, 1.03)
    assert s.shape == (500, 3) and np.allclose(np.linalg.norm(s, axis=1), 1.03)
    assert abs(s.mean(axis=0)).max() < 0.01


def _dipole():
    phi = np.linspace(0, 2 * np.pi, 33)
    th = np.linspace(0, np.pi, 181)
    r = np.linspace(1.0, 3.0, 201)
    P, T, R = np.meshgrid(phi, th, r, indexing="ij")
    return dict(phi=phi, theta=th, r=r, br=2 * np.cos(T) / R**3, bt=np.sin(T) / R**3,
                bp=np.zeros_like(R))


def test_tracer_follows_dipole_field_lines():
    # a dipole line obeys r = L sin^2(theta): from r0 = 1.05 at colatitude 45 deg, L = 2.1 (closed);
    # from colatitude 8 deg, L = 54 (leaves the 3 R_sun box: open)
    t = FieldTracer(_dipole())
    seeds = 1.05 * np.array([unit_vec(30, 45), unit_vec(30, 82)])
    (closed, o1), (open_, o2) = t.field_lines(seeds, ds=0.005, n_steps=1500)
    assert not o1 and o2
    assert np.linalg.norm(closed, axis=1).max() == pytest.approx(2.1, abs=0.02)
    lat = np.degrees(np.arcsin(closed[:, 2] / np.linalg.norm(closed, axis=1)))
    # both ends stop at r = 1, where sin^2(theta) = 1/2.1: latitude -46.4 deg in the south
    assert lat.min() == pytest.approx(-(90 - np.degrees(np.arcsin(np.sqrt(1 / 2.1)))), abs=0.5)
    r, z = np.linalg.norm(closed, axis=1), closed[:, 2]
    L = r / (1 - (z / r) ** 2)
    assert np.ptp(L) < 0.03
