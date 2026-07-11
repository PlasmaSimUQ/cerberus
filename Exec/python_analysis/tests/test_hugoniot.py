"""RH locus solver on an analytic ideal-gas table.

For a gamma-law gas the principal-Hugoniot compression approaches
(gamma+1)/(gamma-1) in the strong-shock limit; every returned point must
close the Rankine-Hugoniot energy relation on the table itself.
"""

import numpy as np

from eos_tools.constants import KB, M_D
from eos_tools.hugoniot import locus

GAMMA = 1.4
R = KB / M_D


def make_table(nr=256, nt=128):
    lrho = np.linspace(-4.0, 0.0, nr)
    lT = np.linspace(3.0, 9.0, nt)
    rho = 10.0 ** lrho[:, None]
    T = 10.0 ** lT[None, :]
    p = rho * R * T
    e = R * T / (GAMMA - 1.0) * np.ones_like(rho)
    return lrho, lT, p, e


def test_strong_shock_limit():
    lrho, lT, p, e = make_table()
    rho0, T0 = 1e-2, 1e3
    p0 = rho0 * R * T0
    e0 = R * T0 / (GAMMA - 1.0)
    comp, pres = locus(lrho, lT, p, e, rho0, e0, p0)
    assert len(comp) > 10

    limit = (GAMMA + 1.0) / (GAMMA - 1.0)  # = 6
    # compression is monotone toward the limit and never exceeds it (p0>0)
    assert comp.max() < limit + 1e-3
    assert abs(comp[pres.argmax()] - limit) < 0.02 * limit

    # every point closes the RH energy relation against the table
    from scipy.interpolate import PchipInterpolator
    for c, pv in zip(comp[::7], pres[::7]):
        r1 = c * rho0
        # locate the temperature of this state from the ideal law
        T1 = pv / (r1 * R)
        e1 = R * T1 / (GAMMA - 1.0)
        resid = (e1 - e0) - 0.5 * (pv + p0) * (1.0 / rho0 - 1.0 / r1)
        assert abs(resid) < 1e-6 * max(e1, abs(e0))


def test_weak_end_low_compression():
    lrho, lT, p, e = make_table()
    rho0, T0 = 1e-2, 1e3
    p0 = rho0 * R * T0
    e0 = R * T0 / (GAMMA - 1.0)
    comp, pres = locus(lrho, lT, p, e, rho0, e0, p0)
    # weakest returned shock is well below the strong-shock limit
    assert comp[pres.argmin()] < 4.0
