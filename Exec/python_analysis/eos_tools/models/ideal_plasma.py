"""Fully ionized ideal plasma: classical ions + ideal Fermi-Dirac electrons.

The hot-extension source AND the Tier-2 entropy anchor (plan §1.2/§1.3):
above ~3e7 K deuterium is fully ionized and weakly coupled at every table
density, but the electrons are NOT classical at high rho (T_F ~ 2e7 K at
1000 g/cc) — hence Fermi-Dirac electrons, not Sackur-Tetrode.

Conventions (g_e = 2):

    n       = g (2 pi m kT / h^2)^(3/2) (2/sqrt(pi)) F_k=1/2(eta)
    u       = (2/sqrt(pi)) g (2 pi m kT/h^2)^(3/2) kT F_3/2(eta),  P = 2u/3
    A/V     = mu n - P            (grand-canonical identity; exact)

with F_k(eta) = Int_0^inf x^k / (1 + exp(x - eta)) dx and mu = eta kT.
eta is solved per (rho, T) by bracketed root finding on F_1/2.

Ions are classical ideal (exact Sackur-Tetrode Helmholtz). Coulomb
(Debye-Hueckel) corrections are deliberately absent in v1; the seam-3 QA
gate against iFPEOS/FPEOS measures whether they are needed.

Unit-tested limits: classical (P -> n kT, S -> Sackur-Tetrode) and
degenerate (P -> (2/5) n E_F) electron gas; entropy-integration recovery.
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from ..constants import H_PLANCK, HBAR, KB, M_E
from .base import HelmholtzModel

G_E = 2.0  # electron spin degeneracy


def fermi_integral(k, eta):
    """F_k(eta) = Int_0^inf x^k/(1+exp(x-eta)) dx, scalar eta."""
    pts = [eta] if eta > 0 else None
    hi = max(60.0, eta + 60.0)
    # tight tolerances: the base-class consistency check differentiates
    # a(rho,T) numerically, which amplifies quadrature noise by ~1/h
    val, _ = quad(lambda x: x ** k / (1.0 + np.exp(x - eta)), 0.0, hi,
                  points=pts, limit=400, epsabs=0.0, epsrel=1e-11)
    return val


def solve_eta(n_e, T, m=M_E, g=G_E):
    """Invert n(eta) for the reduced chemical potential eta = mu/kT."""
    pref = g * (2.0 * np.pi * m * KB * T / H_PLANCK ** 2) ** 1.5 \
        * (2.0 / np.sqrt(np.pi))
    target = n_e / pref
    # brackets: classical tail F ~ e^eta sqrt(pi)/2; degenerate F ~ (2/3) eta^1.5
    lo = min(-60.0, np.log(max(target, 1e-300)) - 5.0)
    hi = max(10.0, 1.6 * (1.5 * target) ** (2.0 / 3.0) + 20.0)
    return brentq(lambda et: fermi_integral(0.5, et) - target, lo, hi,
                  xtol=1e-12, rtol=1e-12)


def electron_gas(n_e, T):
    """(P, u, A_per_V, s_per_V) of the ideal FD electron gas (CGS)."""
    eta = solve_eta(n_e, T)
    pref = G_E * (2.0 * np.pi * M_E * KB * T / H_PLANCK ** 2) ** 1.5 \
        * (2.0 / np.sqrt(np.pi))
    u = pref * KB * T * fermi_integral(1.5, eta)
    P = 2.0 * u / 3.0
    A_V = eta * KB * T * n_e - P
    s_V = (u - A_V) / T
    return P, u, A_V, s_V


class IdealPlasma(HelmholtzModel):
    """Full ionization, charge Z per nucleus of mass m_nucleus."""

    name = "ideal_plasma"

    def __init__(self, m_nucleus, Z=1.0, g_ion=1.0):
        self.m = m_nucleus
        self.Z = Z
        self.g_ion = g_ion

    # --- pieces ---------------------------------------------------------
    def _a_ion(self, rho, T):
        n_i = rho / self.m
        lam3 = (H_PLANCK ** 2 / (2.0 * np.pi * self.m * KB * T)) ** 1.5
        return -(KB * T / self.m) * (np.log(self.g_ion / (n_i * lam3)) + 1.0)

    def _electron(self, rho, T, what):
        n_e = self.Z * np.asarray(rho, float) / self.m

        def one(ne, tt):
            P, u, A_V, s_V = electron_gas(ne, tt)
            return {"P": P, "u": u, "A_V": A_V, "s_V": s_V}[what]

        return np.vectorize(one, otypes=[float])(n_e, T)

    # --- Helmholtz + analytic overrides ---------------------------------
    def a(self, rho, T):
        rho = np.asarray(rho, float)
        return self._a_ion(rho, T) + self._electron(rho, T, "A_V") / rho

    def p(self, rho, T):
        rho = np.asarray(rho, float)
        T = np.asarray(T, float)
        return rho / self.m * KB * T + self._electron(rho, T, "P")

    def s(self, rho, T):
        rho = np.asarray(rho, float)
        T = np.asarray(T, float)
        n_i = rho / self.m
        lam3 = (H_PLANCK ** 2 / (2.0 * np.pi * self.m * KB * T)) ** 1.5
        s_ion = (KB / self.m) * (np.log(self.g_ion / (n_i * lam3)) + 2.5)
        return s_ion + self._electron(rho, T, "s_V") / rho

    def e(self, rho, T):
        rho = np.asarray(rho, float)
        T = np.asarray(T, float)
        return 1.5 * KB * T / self.m + self._electron(rho, T, "u") / rho

    def cv(self, rho, T):
        T = np.asarray(T, float)
        hT = 1.0e-4 * T
        return (self.e(rho, T + hT) - self.e(rho, T - hT)) / (2.0 * hT)


def fermi_energy(n_e):
    """E_F = (hbar^2 / 2 m_e) (3 pi^2 n_e)^(2/3)  [erg]."""
    return HBAR ** 2 / (2.0 * M_E) * (3.0 * np.pi ** 2 * n_e) ** (2.0 / 3.0)


class IdealPlasmaFast(IdealPlasma):
    """Grid-friendly variant: the per-point brentq/quad electron solve is
    replaced by 1-D interpolation in the degeneracy variable.

    eta depends on (n_e, T) only through t = n_e/pref(T), and F_3/2 only
    through eta — both genuinely one-dimensional. We tabulate eta(log t) and
    log F_3/2(eta) once on dense PCHIP nodes spanning far beyond the table's
    reach and interpolate. Accuracy vs the exact path is unit-tested
    (~1e-8 relative); build cost is ~1e3 quadratures, evaluation is fully
    vectorized.
    """

    N_NODES = 500
    LOG_T_RANGE = (-30.0, 14.0)  # log10 of t = n_e/pref; eta(1e14) ~ 4e9

    def __init__(self, m_nucleus, Z=1.0, g_ion=1.0):
        super().__init__(m_nucleus, Z=Z, g_ion=g_ion)
        from scipy.interpolate import PchipInterpolator
        lt = np.linspace(*self.LOG_T_RANGE, self.N_NODES)
        etas = np.empty_like(lt)
        for i, l in enumerate(lt):
            t = 10.0 ** l
            lo = min(-80.0, np.log(t) - 5.0)
            hi = max(10.0, 1.6 * (1.5 * t) ** (2.0 / 3.0) + 20.0)
            etas[i] = brentq(lambda et: fermi_integral(0.5, et) - t, lo, hi,
                             xtol=1e-12, rtol=1e-14)
        self._eta_of_lt = PchipInterpolator(lt, etas)
        lf32 = np.array([np.log(fermi_integral(1.5, et)) for et in etas])
        # F_3/2 is monotone in eta; parameterize by eta node values
        o = np.argsort(etas)
        self._lf32_of_eta = PchipInterpolator(etas[o], lf32[o])

    def _electron(self, rho, T, what):
        rho = np.asarray(rho, float)
        T = np.asarray(T, float)
        n_e = self.Z * rho / self.m
        pref = G_E * (2.0 * np.pi * M_E * KB * T / H_PLANCK ** 2) ** 1.5 \
            * (2.0 / np.sqrt(np.pi))
        lt = np.log10(n_e / pref)
        if np.any(lt < self.LOG_T_RANGE[0]) or np.any(lt > self.LOG_T_RANGE[1]):
            raise ValueError("degeneracy variable outside tabulated range")
        eta = self._eta_of_lt(lt)
        u = pref * KB * T * np.exp(self._lf32_of_eta(eta))
        P = 2.0 * u / 3.0
        if what == "P":
            return P
        if what == "u":
            return u
        A_V = eta * KB * T * n_e - P
        if what == "A_V":
            return A_V
        return (u - A_V) / T  # s_V
