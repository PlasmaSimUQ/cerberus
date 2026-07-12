"""Saha model: NIST parse, ionization limits, hydrogenic benchmark,
thermodynamic consistency of the (P, e) pair."""

import numpy as np
import pytest

from eos_tools.constants import EV_ERG, KB
from eos_tools.models.saha import SahaModel, load_nist_ie_csv, ti_ie_path


def test_nist_parse():
    import os
    if not os.path.exists(ti_ie_path()):
        pytest.skip("NIST Ti CSV not present")
    chi = load_nist_ie_csv(ti_ie_path())
    assert len(chi) == 22
    # first two stages: 6.828 and 13.576 eV; strictly increasing overall
    assert abs(chi[0] / EV_ERG - 6.82812) < 1e-3
    assert abs(chi[1] / EV_ERG - 13.5755) < 1e-3
    assert np.all(np.diff(chi) > 0)
    # K-shell stages in the keV range
    assert chi[-1] / EV_ERG > 6000.0


def hydrogen_like():
    m_h = 1.6735e-24
    return SahaModel(m_h, np.array([13.5984 * EV_ERG]))


def test_limits_hydrogenic():
    s = hydrogen_like()
    rho = 1e-7  # dilute
    assert s.zbar(rho, 300.0) < 1e-10
    assert abs(s.zbar(rho, 3.0e6) - 1.0) < 1e-6
    out_cold = s.eval(rho, 300.0)
    n_at = rho / s.m
    assert abs(out_cold["p"] / (n_at * KB * 300.0) - 1.0) < 1e-8
    out_hot = s.eval(rho, 3.0e6)
    assert abs(out_hot["p"] / (2.0 * n_at * KB * 3.0e6) - 1.0) < 1e-5
    # hot e contains the ionization energy investment
    e_kin = 1.5 * 2.0 * KB * 3.0e6 / s.m
    assert abs(out_hot["e"] - (e_kin + 13.5984 * EV_ERG / s.m)) \
        < 1e-6 * out_hot["e"]


def test_half_ionization_crossing():
    """Zbar = 0.5 where the Saha equation says n_e = S(T) with x0 = x1."""
    s = hydrogen_like()
    rho = 1e-8
    # scan T for the crossing, then verify the Saha relation closes there
    Ts = np.geomspace(5e3, 1e5, 400)
    zb = s.zbar(rho, Ts)
    Tc = np.interp(0.5, zb, Ts)
    n_at = rho / s.m
    logS = s._log_saha_rhs(Tc)[0]
    # at x0 = x1: n_e = 0.5 n_at... and Saha: n1 n_e / n0 = S -> n_e = S
    assert abs(np.exp(logS) / (0.5 * n_at) - 1.0) < 0.15


def test_maxwell_consistency_of_pe_pair():
    """The equilibrium (P, e) surfaces must satisfy the Gruneisen/Maxwell
    identity — the envelope-theorem check that ionization-fraction
    re-solves under differentiation are thermodynamically kosher."""
    from eos_tools.thermo import maxwell_residual
    from eos_tools.materials.titanium import M_TI
    import os
    if not os.path.exists(ti_ie_path()):
        pytest.skip("NIST Ti CSV not present")
    from eos_tools.models.saha import titanium_saha
    s = titanium_saha()
    # T spacing ~0.021 dex = the production 384^2 grid pitch: the metric's
    # PCHIP slopes need to resolve the steep ionization steps (measured:
    # median 2e-2 at 0.087 dex pitch -> 8e-4 here, pure h^2 truncation)
    lrho = np.linspace(-3.0, 0.0, 24)
    lT = np.linspace(4.5, 6.5, 96)
    R, T = np.meshgrid(10.0 ** lrho, 10.0 ** lT, indexing="ij")
    out = s.eval(R, T)
    Rm = maxwell_residual(lrho, lT, out["p"], out["e"])
    assert float(np.median(Rm)) < 2e-3, float(np.median(Rm))


def test_titanium_zbar_shape():
    import os
    if not os.path.exists(ti_ie_path()):
        pytest.skip("NIST Ti CSV not present")
    from eos_tools.models.saha import titanium_saha
    s = titanium_saha()
    rho = 1e-3
    zb = s.zbar(rho, np.array([1e4, 1e5, 1e6, 1e7, 3e8]))
    assert zb[0] < 2.0            # warm: first stages only
    assert 2.0 < zb[1] < 12.0     # 10 eV-ish: L-shell opening
    assert 8.0 < zb[2] < 21.0     # 100 eV-ish
    assert zb[4] > 21.9           # fully stripped
    assert np.all(np.diff(zb) > 0)
