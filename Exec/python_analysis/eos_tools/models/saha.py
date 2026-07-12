"""Saha average-ionization model (chemical picture, ideal mixture).

The titanium hot-side source (Ti SS3-prime, doc/eos_creation_plan.md):
ionization stages j = 0..Z in Saha equilibrium,

    n_{j+1} n_e / n_j = 2 (U_{j+1}/U_j) (2 pi m_e kB T / h^2)^{3/2}
                        exp(-chi_j / kB T)

with the NIST-harvested ionization energies chi_j
(data/raw/nist_ti/, manifest id 'nist-ti-ie'). Zbar is solved per (rho, T)
by bracketed root finding on Zbar = sum_j j x_j(n_e = Zbar n_atom);
populations are evaluated in log space.

Thermodynamics of the equilibrium ideal mixture:

    P = n_atom (1 + Zbar) kB T
    e = (3/2)(1 + Zbar) kB T / m  +  sum_j x_j E_j / m,   E_j = sum_{k<j} chi_k

This (P, e) pair is the consistent derivative set of the ideal-mixture
free energy at its Saha minimum (verified numerically by the Maxwell-
residual unit test — equilibrium re-solves under differentiation are
exactly what the envelope theorem licenses).

v1 simplifications, all recorded: partition functions U_j = 1 (no excited
states; shifts Zbar mildly at low T where the solid model owns the table);
Boltzmann electrons (no degeneracy — the FD ideal-plasma source takes over
at the hot seam, and high-rho/low-T degenerate territory belongs to the
solid model); no continuum lowering / pressure ionization (Saha is
unreliable above ~solid density in the WDM band — the Ti v0 table's
documented accuracy gap pending real WDM data).
"""

import csv
import os
import re

import numpy as np
from scipy.optimize import brentq

from ..constants import EV_ERG, H_PLANCK, KB, M_E
from .. import sources as _sources


def load_nist_ie_csv(path):
    """Parse the NIST ASD ionization-energy CSV (their ="..." quoting).

    Returns chi [erg], ordered by stage (Ti I first = neutral -> Ti II).
    """
    chis = []
    with open(path) as f:
        for row in csv.reader(f):
            if not row or row[0].startswith("Sp."):
                continue
            for cell in row:
                m = re.match(r'^="?([0-9]+(\.[0-9]+)?)"?$', cell.strip())
                if m:
                    chis.append(float(m.group(1)) * EV_ERG)
                    break
    if not chis:
        raise ValueError("no ionization energies parsed from %s" % path)
    return np.asarray(chis)


def ti_ie_path():
    return os.path.join(_sources.repo_root(), "Exec", "testing", "EOS-Table",
                        "data", "raw", "nist_ti",
                        "ti_ionization_energies.csv")


class SahaModel:
    """eval-style p/e/zbar for one element; chi in erg, m_atom in g."""

    name = "saha"

    def __init__(self, m_atom, chi):
        self.m = m_atom
        self.chi = np.asarray(chi, float)
        self.Z = len(self.chi)
        self.E_cum = np.concatenate([[0.0], np.cumsum(self.chi)])  # per stage

    # -- stage populations at given electron density ----------------------
    def _log_saha_rhs(self, T):
        """log S_j, j = 0..Z-1 (units: per cc)."""
        lam = 2.0 * (2.0 * np.pi * M_E * KB * T / H_PLANCK ** 2) ** 1.5
        return np.log(lam) - self.chi / (KB * T)

    def _zbar_pop(self, n_e, T, logS):
        """(Zbar, x_j) for a trial electron density."""
        # log relative populations: log x_{j+1} - log x_j = logS_j - ln n_e
        inc = logS - np.log(n_e)
        logx = np.concatenate([[0.0], np.cumsum(inc)])
        logx -= logx.max()
        x = np.exp(logx)
        x /= x.sum()
        j = np.arange(self.Z + 1)
        return float((j * x).sum()), x

    def zbar(self, rho, T):
        """Mean ionization, solved per point (vectorized wrapper)."""
        def one(r, t):
            n_at = r / self.m
            logS = self._log_saha_rhs(t)

            def g(zb):
                return zb - self._zbar_pop(max(zb, 1e-12) * n_at, t, logS)[0]

            # g(0+) <= 0, g(Z) >= 0; handle the un-ionized limit cheaply
            if g(1e-10) >= 0.0:
                return self._zbar_pop(1e-10 * n_at, t, logS)[0]
            return brentq(g, 1e-10, float(self.Z), xtol=1e-12, rtol=1e-11)

        return np.vectorize(one, otypes=[float])(rho, T)

    def eval(self, rho, T):
        """Dict with p [erg/cc], e [erg/g], zbar; ideal Saha mixture."""
        rho, T = np.broadcast_arrays(np.asarray(rho, float),
                                     np.asarray(T, float))

        def one(r, t):
            n_at = r / self.m
            logS = self._log_saha_rhs(t)

            def g(zb):
                return zb - self._zbar_pop(max(zb, 1e-12) * n_at, t, logS)[0]

            if g(1e-10) >= 0.0:
                zb, x = self._zbar_pop(1e-10 * n_at, t, logS)
            else:
                zb = brentq(g, 1e-10, float(self.Z), xtol=1e-12, rtol=1e-11)
                _, x = self._zbar_pop(zb * n_at, t, logS)
            p = n_at * (1.0 + zb) * KB * t
            e = 1.5 * (1.0 + zb) * KB * t / self.m \
                + float((x * self.E_cum).sum()) / self.m
            return p, e, zb

        p, e, zb = np.vectorize(one, otypes=[float] * 3)(rho, T)
        return {"p": p, "e": e, "zbar": zb, "rho": rho, "T": T}


def titanium_saha():
    from ..materials.titanium import M_TI
    chi = load_nist_ie_csv(ti_ie_path())
    assert len(chi) == 22, "expected all 22 Ti stages, got %d" % len(chi)
    return SahaModel(M_TI, chi)
