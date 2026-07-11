"""Thermodynamic-consistency utilities on gridded (p, e) surfaces.

maxwell_residual: the Grueneisen/Maxwell identity on a (rho, T) grid,

    de/drho|T = (p - T dp/dT|rho) / rho^2

holds exactly for any consistent EOS. The residual map

    R = | dedrho - (p - T dpdT)/rho^2 | / ( max(p, T |dpdT|) / rho^2 )

is the Tier-1 acceptance metric (plan §1.2): measured per source first (the
yardstick — regrid + MC noise give a nonzero floor), then on the blended
surfaces. Slopes come from the same PCHIP construction condition() uses, so
the metric measures the data, not a differencing scheme mismatch.
"""

import numpy as np
from scipy.interpolate import PchipInterpolator


def grid_slopes(lrho, lT, f):
    """(df/drho|T, df/dT|rho) wrt LINEAR rho/T via PCHIP on the log axes."""
    rho = 10.0 ** lrho
    T = 10.0 ** lT
    ln10 = np.log(10.0)
    dfdT = np.empty_like(f)
    dfdrho = np.empty_like(f)
    for i in range(f.shape[0]):
        dfdT[i] = PchipInterpolator(lT, f[i]).derivative()(lT) / (T * ln10)
    for j in range(f.shape[1]):
        dfdrho[:, j] = PchipInterpolator(lrho, f[:, j]).derivative()(lrho) \
            / (rho * ln10)
    return dfdrho, dfdT


def maxwell_residual(lrho, lT, p, e):
    """Relative Maxwell/Grueneisen residual map, shape (n_rho, n_T)."""
    rho = 10.0 ** lrho[:, None]
    T = 10.0 ** lT[None, :]
    dedrho, _ = grid_slopes(lrho, lT, e)
    _, dpdT = grid_slopes(lrho, lT, p)
    lhs = dedrho
    rhs = (p - T * dpdT) / rho ** 2
    scale = np.maximum(np.abs(p), np.abs(T * dpdT)) / rho ** 2
    return np.abs(lhs - rhs) / np.maximum(scale, 1e-300)


def sound_speed_sq(rho, T, p, dpdrho, dpdT, cv):
    """cs^2 = dpdrho|T + (T/rho^2) dpdT^2 / cv  (the D5 combination)."""
    return dpdrho + (T / rho ** 2) * dpdT ** 2 / cv
