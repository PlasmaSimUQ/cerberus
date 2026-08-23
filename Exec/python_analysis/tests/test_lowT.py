"""Sub-floor T extension (doc/eos_air_lowT_extension_plan.md, eos_tools.lowT).

Covers: native nodes untouched and native spacing kept; the three branch
regimes (vapor = ideal gas exactly, condensed = anchored slope, dome =
scaled row); dpdT continuity and positivity; p <= p0 (no native lift under
cummax-in-T); rho-monotonicity when the anchor row is monotone; constant
cv; no-op when the floor is already reached.
"""

import numpy as np
import pytest

from eos_tools.constants import AMU_G, KB
from eos_tools.lowT import extend_T_floor

M_AIR = 28.97


def _anchor_surface(nT=6, lT0=2.0, h=0.017081):
    """Synthetic molecular-gas conditioned surface: ideal vapor at low rho, a
    flat dome (p = p_sat) in the middle, a stiff condensed branch above.
    T axis starts at 100 K (lT0 = 2.0) on a uniform log spacing."""
    lrho = np.linspace(-7, 1.18, 60)
    rho = 10.0 ** lrho
    lT = lT0 + h * np.arange(nT)
    T = 10.0 ** lT
    R_s = KB / (M_AIR * AMU_G)
    p_sat = 8.0e6                                     # 8 bar dome
    p = np.empty((len(rho), nT))
    e = np.empty_like(p)
    for j, Tj in enumerate(T):
        ideal = rho * R_s * Tj
        cond = 3e9 * (rho / 1.0) ** 4                  # stiff above ~1 g/cc
        pj = np.where(ideal < p_sat, ideal, p_sat)    # vapor -> dome
        pj = np.maximum(pj, cond)                      # dome -> condensed
        p[:, j] = pj
        e[:, j] = 2.5 * R_s * Tj + 1e8 * np.log10(rho + 1e-9)  # any smooth e
    hull = np.ones_like(p)
    band = np.zeros(p.shape, dtype=bool)
    return lrho, lT, p, e, hull, band, R_s


def test_native_rows_untouched_and_spacing_kept():
    lrho, lT, p, e, hull, band, _ = _anchor_surface()
    lT2, p2, e2, hull2, band2, st = extend_T_floor(
        lrho, lT, p, e, hull, band, 1.25, M_AIR)
    k = st["rows"]
    assert k == 44                                   # (2.0 - 1.25)/0.017081 -> 44
    assert lT2[0] <= 1.25 and lT2[1] > 1.25 - 1e-12 + st["h"] - st["h"]
    np.testing.assert_allclose(np.diff(lT2), st["h"], rtol=0, atol=1e-12)
    assert np.array_equal(p2[:, k:], p) and np.array_equal(e2[:, k:], e)
    assert np.array_equal(hull2[:, k:], hull) and np.array_equal(band2[:, k:], band)
    assert np.all(hull2[:, :k] == 0.0) and not band2[:, :k].any()


def test_vapor_columns_are_the_ideal_gas_exactly():
    lrho, lT, p, e, hull, band, R_s = _anchor_surface()
    lT2, p2, e2, *_ , st = extend_T_floor(lrho, lT, p, e, hull, band, 1.25, M_AIR)
    k = st["rows"]
    rho = 10.0 ** lrho
    T = 10.0 ** lT2[:k]
    vap = rho * R_s * 100.0 < 8.0e6 * 0.5           # well inside the vapor regime
    ideal = rho[vap][:, None] * R_s * T[None, :]
    np.testing.assert_allclose(p2[vap, :k], ideal, rtol=1e-12)
    # e drops with the ideal cv: constant slope 5/2 R_s in T
    dedT = np.diff(e2[vap, :k + 1], axis=1) / np.diff(10.0 ** lT2[:k + 1])[None, :]
    np.testing.assert_allclose(dedT, 2.5 * R_s, rtol=1e-12)


def test_branches_dome_and_condensed():
    lrho, lT, p, e, hull, band, R_s = _anchor_surface()
    lT2, p2, *_ , st = extend_T_floor(lrho, lT, p, e, hull, band, 1.25, M_AIR)
    k = st["rows"]
    rho = 10.0 ** lrho
    T = 10.0 ** lT2[:k]
    T0 = 10.0 ** lT[0]
    p0 = p[:, 0]
    ideal0 = rho * R_s * T0
    dome = (p0 == 8.0e6)
    # the branch criterion is p0 vs the ideal-gas pressure at T0, not the
    # dome label: a just-exited liquid column with p_sat < p0 < rho R_s T0
    # is still sub-ideal and takes the scaled branch
    scaled = p0 < ideal0 * (1 - 1e-9)
    condensed = p0 > ideal0 * (1 + 1e-9)
    vapor = ~scaled & ~condensed
    assert dome.any() and condensed.any() and np.all(scaled[dome])
    # scaled row: p = p0 T/T0 (ideal gas at the dome-edge vapor density)
    np.testing.assert_allclose(p2[scaled, :k], p0[scaled][:, None] * (T / T0)[None, :],
                               rtol=1e-12)
    # condensed: anchored slope rho R_s, pressure stays cold-curve dominated
    exp = p0[condensed][:, None] - rho[condensed][:, None] * R_s * (T0 - T)[None, :]
    np.testing.assert_allclose(p2[condensed, :k], exp, rtol=1e-12)
    assert np.all(p2[condensed, :k] > 0.5 * p0[condensed][:, None])
    assert st["n_anchor"] == condensed.sum() * k
    assert st["n_scaled"] == scaled.sum() * k
    assert st["n_vapor"] == vapor.sum() * k and vapor.sum() > 0


def test_dpdT_positive_continuous_and_no_native_lift():
    lrho, lT, p, e, hull, band, R_s = _anchor_surface()
    lT2, p2, *_ , st = extend_T_floor(lrho, lT, p, e, hull, band, 1.25, M_AIR)
    k = st["rows"]
    T = 10.0 ** lT2
    rho = 10.0 ** lrho
    dp = np.diff(p2[:, :k + 1], axis=1)
    assert np.all(dp > 0)                                   # strictly increasing in T
    # slope is rho R_s (anchored) or p0/T0 (scaled); equal at the switch,
    # so the FD slope never exceeds max(rho R_s, p0/T0) and is continuous
    slope = dp / np.diff(T[:k + 1])[None, :]
    cap = np.maximum(rho * R_s, p[:, 0] / T[k])[:, None]
    assert np.all(slope <= cap * (1 + 1e-9))
    # a column-wise cummax in T (monotonise_T) cannot move the native row
    assert np.all(p2[:, k - 1] <= p2[:, k])
    assert np.all(np.maximum.accumulate(p2, axis=1)[:, k:] == p2[:, k:])


def test_rho_monotone_when_anchor_row_is():
    lrho, lT, p, e, hull, band, _ = _anchor_surface()
    assert np.all(np.diff(p[:, 0]) >= 0)
    lT2, p2, *_ , st = extend_T_floor(lrho, lT, p, e, hull, band, 1.25, M_AIR)
    k = st["rows"]
    assert np.all(np.diff(p2[:, :k], axis=0) >= -1e-12 * p2[:-1, :k])


def test_noop_when_floor_already_reached():
    lrho, lT, p, e, hull, band, _ = _anchor_surface()
    lT2, p2, e2, hull2, band2, st = extend_T_floor(
        lrho, lT, p, e, hull, band, 2.0, M_AIR)
    assert st["rows"] == 0 and lT2 is lT and p2 is p and e2 is e


def test_requires_ascending_axis():
    lrho, lT, p, e, hull, band, _ = _anchor_surface()
    with pytest.raises(ValueError):
        extend_T_floor(lrho, lT[::-1], p, e, hull, band, 1.25, M_AIR)
