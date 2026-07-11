"""Hull-aware regrid of an exactly-known ideal gas must reproduce it.

p = rho*R*T is log-linear in both log axes, so the PCHIP passes are exact
on it (PCHIP preserves linear data); e = R*T/(gamma-1) is exponential in
log T, so it carries a small interpolation error.
"""

import numpy as np

from eos_tools.constants import KB, M_D
from eos_tools.grids import regrid

GAMMA = 1.4
R = KB / M_D


def ideal_points(rhos, Ts):
    rho, T = np.meshgrid(rhos, Ts, indexing="ij")
    return {
        "rho": rho.ravel(),
        "T": T.ravel(),
        "p": (rho * R * T).ravel(),
        "e": (R * T / (GAMMA - 1.0) * np.ones_like(rho)).ravel(),
    }


def test_regrid_reproduces_ideal_gas():
    pts = ideal_points(np.logspace(-3, 0, 20), np.logspace(3, 7, 12))
    lrho, lT, p, e, hull, loo = regrid(pts, 48, 48)

    # rectangular source -> fully covered hull
    assert hull.mean() == 1.0

    rho = 10.0 ** lrho[:, None]
    T = 10.0 ** lT[None, :]
    assert np.allclose(p, rho * R * T, rtol=1e-12)
    # e is exponential in log T: PCHIP through 12 source points per decade-
    # third carries mid-cell error (measured max 4.1% at this raggedness)
    assert np.allclose(e, R * T / (GAMMA - 1.0) * np.ones_like(rho), rtol=5e-2)

    # leave-one-out stats exist and are small on smooth data
    assert loo.shape[1] == 2
    assert np.median(loo[:, 0]) < 1e-10  # p exact under log-log PCHIP


def test_regrid_marks_uncovered_cells():
    # ragged source: the high-T half only exists at high density
    lo = ideal_points(np.logspace(-3, 0, 16), np.logspace(3, 5, 8))
    hi = ideal_points(np.logspace(-1, 0, 8), np.logspace(5.2, 7, 8))
    pts = {k: np.concatenate([lo[k], hi[k]]) for k in ("rho", "T", "p", "e")}
    lrho, lT, p, e, hull, _ = regrid(pts, 40, 40)

    assert 0.0 < hull.mean() < 1.0
    # low-density high-T corner has no source data behind it
    assert hull[0, -1] == 0.0
    # but is filled with finite values
    assert np.all(np.isfinite(p)) and np.all(np.isfinite(e))
    # the covered low-T region is still exact
    in_lo = (lT <= 5.0)
    assert np.allclose(p[:, in_lo][hull[:, in_lo] > 0.5],
                       (10.0 ** lrho[:, None] * R * 10.0 ** lT[None, in_lo]
                        )[hull[:, in_lo] > 0.5], rtol=1e-10)
