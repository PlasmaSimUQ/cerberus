"""Derivative blocks from condition() vs closed-form ideal-gas derivatives.

PCHIP node derivatives on smooth exponential data carry O(h^2) error
(~0.6% at this resolution) in the interior and one-sided error at the
axis edges, hence the split tolerances.
"""

import numpy as np

from eos_tools.condition import condition, inverse_maps, shift_energy
from eos_tools.constants import KB, M_D

GAMMA = 1.4
R = KB / M_D


def ideal_surfaces(nr=48, nt=48):
    lrho = np.linspace(-3.0, 0.0, nr)
    lT = np.linspace(3.0, 7.0, nt)
    rho = 10.0 ** lrho[:, None]
    T = 10.0 ** lT[None, :]
    p = rho * R * T
    e = R * T / (GAMMA - 1.0) * np.ones_like(rho)
    return lrho, lT, p, e


def test_condition_matches_analytic():
    lrho, lT, p, e = ideal_surfaces()
    dpdT, dpdrho, cv, dedrho, stats = condition(lrho, lT, p, e)

    rho = 10.0 ** lrho[:, None]
    T = 10.0 ** lT[None, :]
    interior = np.s_[1:-1, 1:-1]

    assert np.allclose(dpdT[interior], (rho * R * np.ones_like(T))[interior],
                       rtol=2e-2)
    assert np.allclose(dpdrho[interior], (R * T * np.ones_like(rho))[interior],
                       rtol=2e-2)
    assert np.allclose(cv[interior],
                       (R / (GAMMA - 1.0) * np.ones_like(p))[interior],
                       rtol=2e-2)
    # e is rho-independent -> exactly zero slope through PCHIP
    assert np.allclose(dedrho, 0.0, atol=1e-12 * np.abs(e).max())

    # edges: one-sided derivatives, looser
    assert np.allclose(dpdT, rho * R * np.ones_like(T), rtol=1e-1)
    assert np.allclose(cv, R / (GAMMA - 1.0) * np.ones_like(p), rtol=1e-1)

    assert stats["cv_floored"] == 0
    assert stats["monotonised"] == 0
    assert stats["maxwell"] == "none-needed"


def test_cv_floor_applies():
    lrho, lT, p, e = ideal_surfaces(16, 16)
    e_const = np.full_like(e, e.mean())  # cv = 0 everywhere
    _, _, cv, _, stats = condition(lrho, lT, p, e_const)
    assert stats["cv_floored"] == cv.size
    assert np.all(cv >= stats["cv_floor"])


def test_inverse_maps_seed_quality():
    lrho, lT, p, e = ideal_surfaces()
    le, lp, T_of_e, T_of_p = inverse_maps(lrho, lT, p, e)
    # seed check: T(rho, e(rho,T0)) ~ T0 for an interior probe
    i = 20
    j = 25
    T = 10.0 ** lT
    Ts = np.interp(np.log10(e[i, j]), le, np.arange(len(le)))
    seed = T_of_e[i, int(round(Ts))]
    assert abs(seed / T[j] - 1.0) < 0.15  # seed-only quality


def test_shift_energy_positive_on_hull():
    e = np.array([[-3.0, 1.0], [2.0, 5.0]])
    hull = np.ones_like(e)
    es, shift = shift_energy(e, hull)
    assert shift > 3.0
    assert es.min() > 0.0
    # already-positive e is untouched
    es2, shift2 = shift_energy(e + 10.0, hull)
    assert shift2 == 0.0
    assert np.array_equal(es2, e + 10.0)
