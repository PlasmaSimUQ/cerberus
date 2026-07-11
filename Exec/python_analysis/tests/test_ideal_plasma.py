"""Ideal-plasma model: classical and degenerate electron limits, and the
grand-canonical Helmholtz consistency (base differencing vs analytic)."""

import numpy as np

from eos_tools.constants import H_PLANCK, KB, M_D, M_E
from eos_tools.models.base import HelmholtzModel
from eos_tools.models.ideal_plasma import (IdealPlasma, electron_gas,
                                           fermi_energy)


def test_electron_classical_limit():
    n_e, T = 1.0e19, 1.0e6  # kT/E_F ~ 4e4
    P, u, A_V, s_V = electron_gas(n_e, T)
    assert abs(P / (n_e * KB * T) - 1.0) < 1e-2
    # Sackur-Tetrode, g = 2
    lam3 = (H_PLANCK ** 2 / (2.0 * np.pi * M_E * KB * T)) ** 1.5
    s_st = n_e * KB * (np.log(2.0 / (n_e * lam3)) + 2.5)
    assert abs(s_V / s_st - 1.0) < 1e-2


def test_electron_degenerate_limit():
    n_e = 1.0e26  # ~340 g/cc deuterium equivalent
    EF = fermi_energy(n_e)
    T = 0.01 * EF / KB
    P, u, A_V, s_V = electron_gas(n_e, T)
    assert abs(P / (0.4 * n_e * EF) - 1.0) < 1e-3
    assert abs(u / (0.6 * n_e * EF) - 1.0) < 1e-3
    assert s_V * T / u < 0.05  # entropy nearly frozen out


def test_model_thermo_consistency():
    """Analytic p/e/s must agree with differencing the model's own a."""
    m = IdealPlasma(M_D, Z=1.0)
    rho = np.array([1e-3, 0.171, 10.0])
    T = np.array([5.0e7, 1.0e8, 4.0e8])
    p_num = HelmholtzModel.p(m, rho, T)
    s_num = HelmholtzModel.s(m, rho, T)
    assert np.allclose(p_num, m.p(rho, T), rtol=1e-5)
    assert np.allclose(s_num, m.s(rho, T), rtol=1e-5)
    # e = a + T s identity with all-analytic pieces
    assert np.allclose(m.e(rho, T), m.a(rho, T) + T * m.s(rho, T), rtol=1e-9)


def test_full_ionization_pressure_scale():
    """At 1e8 K, 0.171 g/cc: P ~ 2 n_i k T (ions + non-degenerate electrons)."""
    m = IdealPlasma(M_D, Z=1.0)
    rho, T = 0.171, 1.0e8
    n_i = rho / M_D
    p = float(m.p(rho, T))
    assert abs(p / (2.0 * n_i * KB * T) - 1.0) < 0.02
