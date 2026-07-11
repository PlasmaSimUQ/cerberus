"""Vinet cold-curve identities + fit recovery + LLNL data smoke-fit."""

import os

import numpy as np

from eos_tools.models.coldcurve import (VinetColdCurve, fit_vinet, vinet_e,
                                        vinet_p)

V0, B0, B0P = 5.7, 1.7e11, 6.9  # arbitrary but hydrogen-ish (cc/g, barye)


def test_vinet_identities():
    # P(v0) = 0
    assert abs(vinet_p(V0, V0, B0, B0P)) < 1e-6 * B0
    # bulk modulus at v0: B = -v dP/dv -> B0
    h = 1e-6 * V0
    B = -V0 * (vinet_p(V0 + h, V0, B0, B0P) - vinet_p(V0 - h, V0, B0, B0P)) / (2 * h)
    assert abs(B / B0 - 1.0) < 1e-6
    # e(v0) = e0
    assert abs(vinet_e(V0, V0, B0, B0P, e0=3.0) - 3.0) < 1e-20


def test_energy_is_integral_of_pressure():
    # P = -de/dv everywhere (the closed form must integrate the P form)
    v = np.linspace(0.2 * V0, 1.4 * V0, 300)
    h = 1e-6 * v
    dedv = (vinet_e(v + h, V0, B0, B0P) - vinet_e(v - h, V0, B0, B0P)) / (2 * h)
    P = vinet_p(v, V0, B0, B0P)
    scale = np.abs(P).max()
    assert np.allclose(-dedv, P, atol=1e-8 * scale)


def test_fit_recovers_planted_parameters():
    v = np.linspace(0.25 * V0, V0, 60)
    P = vinet_p(v, V0, B0, B0P)
    (v0f, B0f, B0pf), rms = fit_vinet(P, v)
    assert abs(v0f / V0 - 1) < 1e-3
    assert abs(B0f / B0 - 1) < 1e-3
    assert abs(B0pf / B0P - 1) < 1e-3
    assert rms < 1e-6


def test_model_consistency():
    m = VinetColdCurve(V0, B0, B0P)
    rho = np.array([0.2, 0.5, 1.0])
    T = np.array([100.0, 1000.0, 10000.0])
    # temperature-independent, zero entropy/cv
    assert np.all(m.s(rho, T) == 0.0)
    assert np.all(m.cv(rho, T) == 0.0)
    # p override agrees with differentiating a (base-class path)
    from eos_tools.models.base import HelmholtzModel
    p_base = HelmholtzModel.p(m, rho, T)
    assert np.allclose(p_base, m.p(rho, T), rtol=1e-6)


def test_llnl_coldcurve_smoke_fit():
    """Single Vinet over the merged LLNL curve (exp-6 + Vinet regions)."""
    csv = os.path.join(os.path.dirname(__file__), "..", "..", "testing",
                       "EOS-Table", "data", "raw", "llnl_coldcurve",
                       "h_coldcurve_0K.csv")
    csv = os.path.abspath(csv)
    if not os.path.exists(csv):
        import pytest
        pytest.skip("LLNL cold-curve CSV not present")
    from eos_tools.constants import GPA_CGS, M_D, NA
    P_GPa, V_mol = np.loadtxt(csv, delimiter=",", skiprows=4, unpack=True)
    P = P_GPa * GPA_CGS
    v = V_mol / (M_D * NA)  # cm3/mol-atom -> cc/g (deuterium mass basis)
    (v0f, B0f, B0pf), rms = fit_vinet(P, v)
    # merged exp-6+Vinet source: a single Vinet won't be exact; measured
    # first-pass rms is committed here as the regression number
    assert rms < 0.15, rms
    assert 0.8 / 0.20 < v0f < 1.0 / 0.15  # v0 ~ 1/rho0, rho0 ~ 0.17 g/cc
