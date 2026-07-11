"""Rankine-Hugoniot principal-locus solver on a gridded (p, e) table.

Factored out of the QA plotting so SS3 acceptance gates and the splice QA
can share one canonical offline solver. (The test-suite check.py scripts
keep their OWN independent solvers on purpose — do not import this there.)
"""

import numpy as np
from scipy.interpolate import PchipInterpolator


def locus(lrho, lT, p, e, rho0, e0, p0):
    """Principal Hugoniot from the table and an initial state.

    Parameterise the locus by T (its natural parameter): at each grid
    temperature, solve for the density where the Rankine-Hugoniot
    energy relation closes. Iterating over grid *densities* instead
    undersamples badly — the whole in-table locus spans only ~0.1 dex
    in rho (compression saturates near 4.7). Points whose Hugoniot T
    lies below the table floor (the low-P foot) are absent by
    construction.

    e0 must be on the same energy zero-point as the table's e block
    (i.e. include the table's e_shift).
    Returns (compression rho/rho0, p) arrays sorted by pressure.
    """
    hr, hp = [], []
    for j in range(len(lT)):
        fe = PchipInterpolator(lrho, e[:, j])
        fp = PchipInterpolator(lrho, p[:, j])

        def resid(lr):
            return (fe(lr) - e0) - 0.5 * (fp(lr) + p0) * (
                1.0 / rho0 - 1.0 / 10.0 ** lr)

        lo_lim = np.log10(rho0 * 1.02)
        m = lrho > lo_lim
        if m.sum() < 2:
            continue
        lrs = lrho[m]
        rvals = resid(lrs)
        for k in np.where(rvals[:-1] * rvals[1:] <= 0)[0]:
            lo, hi = lrs[k], lrs[k + 1]
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                if resid(lo) * resid(mid) <= 0:
                    hi = mid
                else:
                    lo = mid
            root = 0.5 * (lo + hi)
            hr.append(10.0 ** root / rho0)
            hp.append(float(fp(root)))
    o = np.argsort(hp)
    return np.array(hr)[o], np.array(hp)[o]
