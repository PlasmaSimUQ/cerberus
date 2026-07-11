"""QEOS-style ion-thermal term: Debye solid with Grueneisen density scaling.

Specific Helmholtz per gram for 3 modes per nucleus (QEOS convention):

    a_ion(rho, T) = R [ 3 T ln(1 - exp(-th/T)) - T D3(th/T) ]  (+ 9/8 R th)

with R = kB/m_nucleus, Debye temperature th(rho) = th0 (rho/rho0)^gG, and
D3 the third Debye function D3(y) = (3/y^3) Int_0^y t^3/(e^t - 1) dt.

The zero-point term (9/8) R th is OFF by default: the LLNL cold curve this
pairs with was reduced from experiment and already carries the physical
zero-point contribution (see data/raw/llnl_coldcurve/).

Known limits (unit-tested): cv -> 3R for T >> th; cv ~ (12 pi^4/5) R (T/th)^3
for T << th. Electron-thermal (Thomas-Fermi) is deliberately absent in the
deuterium cold model — electronic excitation is negligible below the seam-1
band (~0.2 eV vs a ~10 eV gap); a TF term is a titanium-stage addition.
"""

import numpy as np
from scipy.integrate import quad

from ..constants import KB
from .base import HelmholtzModel


def debye3(y):
    """D3(y), vectorized; series-safe at small y, quad elsewhere."""
    y = np.asarray(y, float)

    def one(yy):
        if yy < 1e-6:
            return 1.0 - 3.0 * yy / 8.0 + yy ** 2 / 20.0
        val, _ = quad(lambda t: t ** 3 / np.expm1(t), 0.0, yy, limit=200)
        return 3.0 * val / yy ** 3

    return np.vectorize(one, otypes=[float])(y)


class DebyeIonThermal(HelmholtzModel):
    name = "debye_ion"

    def __init__(self, m_nucleus, theta0, rho0, gruneisen, include_zpe=False):
        self.R = KB / m_nucleus
        self.theta0 = theta0
        self.rho0 = rho0
        self.gG = gruneisen
        self.include_zpe = include_zpe

    def theta(self, rho):
        return self.theta0 * (np.asarray(rho, float) / self.rho0) ** self.gG

    def a(self, rho, T):
        T = np.asarray(T, float)
        th = self.theta(rho)
        y = th / T
        out = self.R * T * (3.0 * np.log(-np.expm1(-y)) - debye3(y))
        if self.include_zpe:
            out = out + 9.0 / 8.0 * self.R * th
        return out


class SlaterDebye(DebyeIonThermal):
    """Debye ion thermal with theta(rho) DERIVED from a cold curve:

        theta(rho) = (hbar/kB) (6 pi^2 n)^(1/3) sqrt(B_cold(rho)/rho)

    (Slater construction — no hand-entered Debye constants; B floored at
    1e-4 B0 in the tension region where the Vinet bulk modulus dips
    through zero). Shared by the deuterium and titanium materials."""

    def __init__(self, m_nucleus, cold):
        super().__init__(m_nucleus, theta0=1.0, rho0=1.0, gruneisen=0.0)
        self.cold = cold

    def theta(self, rho):
        from ..constants import HBAR
        from .coldcurve import vinet_p
        rho = np.asarray(rho, float)
        v = 1.0 / rho
        h = 1e-6 * v
        dPdv = (vinet_p(v + h, self.cold.v0, self.cold.B0, self.cold.B0p)
                - vinet_p(v - h, self.cold.v0, self.cold.B0, self.cold.B0p)) \
            / (2.0 * h)
        B = np.maximum(-v * dPdv, 1e-4 * self.cold.B0)
        n = rho / (KB / self.R)  # nuclei per cc (R = kB/m)
        return HBAR / KB * (6.0 * np.pi ** 2 * n) ** (1.0 / 3.0) \
            * np.sqrt(B / rho)
