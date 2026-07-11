"""T=0 cold curve: Vinet form + fit to tabulated V(P) data.

Vinet EOS in specific volume v = 1/rho, with x = (v/v0)^(1/3),
eta = 1.5 (B0p - 1):

    P(v)      = 3 B0 (1 - x) x^-2 exp(eta (1 - x))
    e_cold(v) = e0 + (9 B0 v0 / eta^2) [1 - (1 - eta (1 - x)) exp(eta (1 - x))]

e_cold is the exact integral of -P dv from v0 (unit test checks this
against numerical quadrature). The cold curve is temperature-independent,
so a = e_cold and s = cv = 0.

Units: B0 in erg/cc (barye), v0 in cc/g, e0 in erg/g.
"""

import numpy as np
from scipy.optimize import curve_fit

from .base import HelmholtzModel


def vinet_p(v, v0, B0, B0p):
    x = np.cbrt(np.asarray(v, float) / v0)
    eta = 1.5 * (B0p - 1.0)
    return 3.0 * B0 * (1.0 - x) / x ** 2 * np.exp(eta * (1.0 - x))


def vinet_e(v, v0, B0, B0p, e0=0.0):
    x = np.cbrt(np.asarray(v, float) / v0)
    eta = 1.5 * (B0p - 1.0)
    z = eta * (1.0 - x)
    return e0 + (9.0 * B0 * v0 / eta ** 2) * (1.0 - (1.0 - z) * np.exp(z))


def fit_vinet(P, v, p_weight_floor=None):
    """Least-squares Vinet fit to (P [erg/cc], v [cc/g]) points.

    Fits in log(P + shift) to balance the ~5-decade pressure span; the
    P = 0 point pins v0 through the residual at small P. Returns
    (v0, B0, B0p), rms relative P error over P > 0 points.
    """
    P = np.asarray(P, float)
    v = np.asarray(v, float)
    scale = p_weight_floor if p_weight_floor else 0.02 * P.max()

    def model(vv, v0, B0, B0p):
        # clip: curve_fit explores parameter sets where P + scale < 0
        return np.log(np.maximum(vinet_p(vv, v0, B0, B0p) + scale,
                                 1e-30 * scale))

    p0 = (v[np.argmin(np.abs(P))], 1e10, 7.0)
    popt, _ = curve_fit(model, v, np.log(P + scale), p0=p0, maxfev=20000)
    v0, B0, B0p = popt
    m = P > 0
    rel = (vinet_p(v[m], v0, B0, B0p) - P[m]) / P[m]
    return (v0, B0, B0p), float(np.sqrt(np.mean(rel ** 2)))


class VinetColdCurve(HelmholtzModel):
    """Cold (T=0) Helmholtz term; s = cv = 0 exactly."""

    name = "vinet_cold"

    def __init__(self, v0, B0, B0p, e0=0.0):
        self.v0, self.B0, self.B0p, self.e0 = v0, B0, B0p, e0

    def a(self, rho, T):
        v = 1.0 / np.asarray(rho, float)
        return vinet_e(v, self.v0, self.B0, self.B0p, self.e0) \
            * np.ones_like(np.asarray(T, float))

    # analytic overrides (exact; avoids differencing a T-constant)
    def s(self, rho, T):
        return np.zeros(np.broadcast(np.asarray(rho), np.asarray(T)).shape)

    def cv(self, rho, T):
        return np.zeros(np.broadcast(np.asarray(rho), np.asarray(T)).shape)

    def p(self, rho, T):
        v = 1.0 / np.asarray(rho, float)
        return vinet_p(v, self.v0, self.B0, self.B0p) \
            * np.ones_like(np.asarray(T, float))
