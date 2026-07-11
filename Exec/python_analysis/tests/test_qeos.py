"""Debye ion-thermal limits + Mie-Grueneisen consistency."""

import numpy as np

from eos_tools.constants import M_D
from eos_tools.models.qeos import DebyeIonThermal, debye3

R = 1.380649e-16 / M_D
THETA0, RHO0, GG = 120.0, 0.2, 1.0


def make():
    return DebyeIonThermal(M_D, THETA0, RHO0, GG)


def test_debye3_limits():
    assert abs(debye3(1e-9) - 1.0) < 1e-8
    # large y: D3 -> (3/y^3) * pi^4/15
    y = 60.0
    assert abs(debye3(y) / (3.0 * np.pi ** 4 / 15.0 / y ** 3) - 1.0) < 1e-6


def test_cv_high_T_dulong_petit():
    m = make()
    cv = m.cv(RHO0, 50.0 * THETA0)
    assert abs(cv / (3.0 * R) - 1.0) < 1e-3


def test_cv_low_T_cubed():
    m = make()
    T = THETA0 / 50.0
    cv = m.cv(RHO0, T)
    cv_debye = 12.0 * np.pi ** 4 / 5.0 * R * (T / THETA0) ** 3
    assert abs(cv / cv_debye - 1.0) < 0.05


def test_mie_gruneisen_pressure():
    """theta ~ rho^gG implies P_ion = gG * rho * e_ion exactly."""
    m = make()
    rho = np.array([0.1, 0.3, 1.0])
    T = np.array([30.0, 150.0, 600.0])
    p = m.p(rho, T)
    e = m.e(rho, T)
    assert np.allclose(p, GG * rho * e, rtol=1e-5)


def test_entropy_positive_and_increasing_in_T():
    m = make()
    T = np.array([5.0, 20.0, 80.0, 300.0, 1200.0])
    s = m.s(0.25, T)
    assert np.all(s > 0.0)
    assert np.all(np.diff(s) > 0.0)
