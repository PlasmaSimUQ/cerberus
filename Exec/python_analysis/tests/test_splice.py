"""Splice building blocks: weight partition, fast plasma evaluator."""

import numpy as np
import pytest

from eos_tools.constants import M_D
from eos_tools.splice import nominal_weights, rho_capable


def grids():
    lrho = np.linspace(-4.0, 3.0, 40)
    lT = np.linspace(1.25, 9.0, 45)
    R = 10.0 ** lrho[:, None] * np.ones((1, len(lT)))
    TT = 10.0 ** lT[None, :] * np.ones((len(lrho), 1))
    return R, TT


def test_weights_partition_of_unity():
    R, TT = grids()
    W = nominal_weights(R, TT)
    assert W.shape[0] == 4
    assert np.all(W >= 0.0)
    assert np.allclose(W.sum(axis=0), 1.0, atol=1e-12)


def test_weights_limits():
    R, TT = grids()
    W = nominal_weights(R, TT)
    # cold dominates the cold/low-rho corner
    assert W[0][0, 0] > 0.99
    # ideal plasma owns the hot end at any density
    assert np.all(W[3][:, -1] > 0.99)
    # REOS.3 owns the mid-T band (e.g. 1e4 K) everywhere in rho
    j = np.argmin(np.abs(np.log10(TT[0]) - 4.0))
    assert np.all(W[1][:, j] > 0.99)
    # high-rho low-T belongs to REOS.3, not the cold model
    assert W[0][-1, 0] < 1e-6 and W[1][-1, 0] > 0.99


def test_rho_capable_constant_in_T():
    R, TT = grids()
    rc = rho_capable(R, TT)
    assert np.allclose(rc, rc[:, :1], atol=1e-14)  # no T dependence (v3)
    assert rc[0, 0] > 0.999 and rc[-1, 0] < 1e-3


def test_fast_plasma_matches_exact():
    from eos_tools.models.ideal_plasma import IdealPlasma, IdealPlasmaFast
    fast = IdealPlasmaFast(M_D, Z=1.0)
    exact = IdealPlasma(M_D, Z=1.0)
    rho = np.geomspace(1e-4, 1000.0, 6)
    T = np.full_like(rho, 5.0e7)
    for f, ex in ((fast.p, exact.p), (fast.e, exact.e), (fast.s, exact.s)):
        assert np.allclose(f(rho, T), ex(rho, T), rtol=1e-5)


def test_fast_plasma_range_guard():
    from eos_tools.models.ideal_plasma import IdealPlasmaFast
    fast = IdealPlasmaFast(M_D, Z=1.0)
    with pytest.raises(ValueError):
        fast.p(1e30, 10.0)  # absurd degeneracy -> outside tabulated range
