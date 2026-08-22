"""Quiet-foot crossover shape (doc/eos_quiet_foot_plan.md changes A/B).

Covers: floor-hugging bridge geometry, legacy byte-path when p_foot is
None, the anchor-cap fallback, ramp-only (allow_maxwell=False) routing,
and the strict-dpdT T-tilt across a conditioned surface.
"""

import numpy as np
import pytest

from eos_tools.maxwell import (KNEE_GAP, P_FOOT_TILT, condition_surface,
                               crossover_isotherm)


def _clipped_isotherm(n_flat=40, placeholder=1e-70):
    """311-style isotherm: placeholder flat through rho0, then a steep
    genuine rise (the Al/Ti low-T structure)."""
    rho = np.logspace(-5, np.log10(2.7), n_flat)
    rho = np.concatenate([rho, [2.727, 2.781, 2.835, 3.0, 3.5]])
    p = np.full(len(rho), placeholder)
    p[n_flat:] = [8.1e8, 1.7e10, 3.4e10, 8.0e10, 3.0e11]  # cgs
    e = np.linspace(1e9, 5e10, len(rho))
    return rho, p, e


def test_hug_places_foot_at_knee():
    rho, p, e = _clipped_isotherm()
    p_foot = 1e6  # 1 bar
    p2, e2, mask, info = crossover_isotherm(rho, p, e, p_foot=p_foot)
    knee = 39  # last flat cell (rho0)
    assert p2[knee] == pytest.approx(p_foot, rel=1e-9)
    # hug segment stays at/below the foot; the rise is the final cell
    assert (p2[1:knee + 1] <= p_foot * (1 + 1e-12)).all()
    assert p2[knee + 1] == pytest.approx(8.1e8)  # anchor untouched
    assert (np.diff(p2) > 0).all()               # strictly monotone
    assert (p2 > 0).all()


def test_legacy_shape_when_foot_absent():
    rho, p, e = _clipped_isotherm()
    p2, _, _, _ = crossover_isotherm(rho, p, e)  # p_foot=None
    knee = 39
    # legacy log-linear span parks anchor-scale pressure at the knee —
    # the defect this feature exists to remove; assert the contrast
    assert p2[knee] > 1e2 * 1e6                  # >> 1 bar (measured ~700 bar)
    assert (np.diff(p2) > 0).all()


def test_anchor_cap_falls_back_to_legacy():
    rho, p, e = _clipped_isotherm()
    p_anchor = 8.1e8
    p_foot = p_anchor / (KNEE_GAP * 0.99)        # too close to the anchor
    p2, _, _, _ = crossover_isotherm(rho, p, e, p_foot=p_foot)
    legacy, _, _, _ = crossover_isotherm(rho, p, e)
    np.testing.assert_allclose(p2, legacy)


def test_ramp_only_suppresses_maxwell():
    # a positive van der Waals loop that WOULD get an equal-area tie line
    rho = np.logspace(-3, 0.5, 60)
    v = 1.0 / rho[::-1]
    pv = 1e8 * (1.0 + 0.5 * np.sin(np.linspace(0, 3 * np.pi, 60))) \
        * np.linspace(1, 30, 60)
    p = pv[::-1].copy()
    e = np.linspace(1e9, 5e10, 60)
    p_mx, _, _, info_mx = crossover_isotherm(rho, p, e, allow_maxwell=True)
    p_rp, _, _, info_rp = crossover_isotherm(rho, p, e, allow_maxwell=False)
    assert info_rp["route"] == "ramp" and info_rp["pstar"] is None
    assert (np.diff(p_rp) > 0).all()
    if info_mx["route"] == "maxwell":            # the loop was tie-lined
        assert not np.allclose(p_mx, p_rp)


def test_tilt_makes_foot_strictly_increasing_in_T():
    rho, p1, e1 = _clipped_isotherm()
    nT = 6
    p = np.tile(p1[:, None], (1, nT))
    e = np.tile(e1[:, None], (1, nT))
    T = np.logspace(1.5, 3, nT)
    p_foot = 1e6
    p2, e2, mask, stats = condition_surface(rho, T, p, e, p_foot=p_foot)
    knee = 39
    foot = p2[knee, :]
    assert (np.diff(foot) > 0).all()             # strict dpdT at the foot
    expect = p_foot * (1.0 + P_FOOT_TILT * np.arange(nT) / (nT - 1))
    np.testing.assert_allclose(foot, expect, rtol=1e-9)
