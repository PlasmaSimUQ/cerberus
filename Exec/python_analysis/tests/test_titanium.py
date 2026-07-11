"""Titanium solid model (SS5b): fit sanity, Sommerfeld limits, consistency."""

import os

import numpy as np
import pytest

from eos_tools.constants import GPA_CGS, NA
from eos_tools.materials.titanium import (M_TI, SommerfeldElectrons,
                                          TitaniumColdModel,
                                          coldcurve_csv_path)


@pytest.fixture(scope="module")
def model():
    if not os.path.exists(coldcurve_csv_path()):
        pytest.skip("LLNL Ti cold-curve CSV not present")
    return TitaniumColdModel()


def test_fit_sanity(model):
    v0, B0, B0p = model.vinet
    assert model.vinet_rms < 0.05
    # pure-Vinet source: the fit's rho0 must be at the data's P=0 point
    assert abs(1.0 / v0 / 4.580 - 1.0) < 0.01
    # bulk modulus in the DAC-established window (no memory constant as an
    # input — this range is a sanity fence, an order-of-magnitude check)
    assert 50.0 < B0 / GPA_CGS < 200.0


def test_slater_theta_reasonable(model):
    th0 = model.theta0()
    assert 200.0 < th0 < 800.0  # literature ~420 K; Slater ignores shear


def test_sommerfeld_limits():
    e = SommerfeldElectrons(M_TI, 4.0)
    rho = 4.51
    T = np.array([100.0, 300.0, 1000.0])
    # cv_e = gamma T exactly (a = -gamma T^2/2); rtol reflects the base
    # class's central second difference — cancellation noise ~5e-7 relative
    # at H_REL = 1e-5, not a physics error
    cv = e.cv(rho, T)
    gam = e.gamma_spec(rho)
    assert np.allclose(cv, gam * T, rtol=1e-4)
    # entropy s_e = gamma T as well (Sommerfeld identity)
    assert np.allclose(e.s(rho, T), gam * T, rtol=1e-6)
    # magnitude: free-electron gamma within a factor ~3 of the ~3.5
    # mJ/mol/K^2 literature value (d-band enhancement expected)
    gam_mol = gam * M_TI * NA * 1e-7 * 1e3  # erg/g/K^2 -> mJ/mol/K^2
    assert 0.5 < gam_mol < 5.0, gam_mol


def test_model_consistency(model):
    rho = np.array([4.6, 6.0, 10.0])
    T = np.array([80.0, 400.0, 1500.0])
    # e = a + T s through the summed model
    assert np.allclose(model.e(rho, T),
                       model.a(rho, T) + T * model.s(rho, T), rtol=1e-9)
    # ion + electron entropy positive, cv positive
    assert np.all(model.s(rho, T) > 0.0)
    assert np.all(model.cv(rho, T) > 0.0)
    # compressed solid: pressure grows with rho at fixed T
    p = model.p(np.array([5.0, 7.0, 12.0]), np.full(3, 300.0))
    assert np.all(np.diff(p) > 0.0)
