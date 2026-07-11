"""Titanium cold/low-T solid model (plan SS5b, the Ti SS2-prime).

Structure mirrors the deuterium QEOS piece on the proven machinery:

    a_Ti(rho, T) = e_cold(rho)                 [Vinet fit, LLNL Ti 0K curve]
                 + a_ion(rho, T)               [Slater-Debye, theta from the
                                                cold curve's bulk modulus]
                 + a_e(rho, T)                 [Sommerfeld electronic term]

The Sommerfeld term a_e = -gamma_e(rho) T^2 / 2 uses the free-electron
Fermi-Dirac gamma at n_e = z_c * n_atom (z_c = 4, titanium's valence
count; no fitted constants). KNOWN LIMITATION, recorded: transition-metal
d-band density of states enhances the true gamma_e(Ti) by ~2-3x over
free-electron — acceptable below ~2 kK where the electronic term is
<= 15% of the ion thermal energy, and the eventual hot-side source owns
everything above. Valid range of this model: solid Ti below melt
(1941 K at ambient, higher under pressure); no fluid/melt piece in v1
(melt smearing happens at the eventual splice seam, as for deuterium).

No memory-sourced physics constants: the cold curve is the harvested LLNL
column; theta(rho) is derived (Slater); z_c is a valence count. The
literature Debye temperature (~420 K) and bulk modulus (~110 GPa) are QA
*reports* for cross-checking, not inputs.
"""

import os

import numpy as np

from ..constants import GPA_CGS, KB, NA
from ..models.base import HelmholtzModel, SumModel
from ..models.coldcurve import VinetColdCurve, fit_vinet
from ..models.ideal_plasma import fermi_energy
from ..models.qeos import SlaterDebye
from .. import sources as _sources

M_TI = 47.867 / NA  # g per atom
Z_COND = 4.0        # conduction-electron count (3d^2 4s^2)

RHO0_AMBIENT = 4.506  # g/cc at 300 K (sanity report; 0 K curve gives 4.580)


def coldcurve_csv_path():
    return os.path.join(_sources.repo_root(), "Exec", "testing", "EOS-Table",
                        "data", "raw", "llnl_coldcurve",
                        "ti_coldcurve_0K.csv")


class SommerfeldElectrons(HelmholtzModel):
    """a_e = -gamma_e(rho) T^2/2, free-electron FD gamma at z_c e/atom.

    gamma_e per volume = pi^2/2 * n_e kB^2 / E_F(n_e); valid T << T_F
    (T_F ~ 10^5 K at solid Ti density, so the quadratic form holds far
    beyond the model's own solid-phase range).
    """

    name = "sommerfeld_e"

    def __init__(self, m_atom, z_c):
        self.m = m_atom
        self.z_c = z_c

    def gamma_spec(self, rho):
        rho = np.asarray(rho, float)
        n_e = self.z_c * rho / self.m
        return np.pi ** 2 / 2.0 * n_e * KB ** 2 / fermi_energy(n_e) / rho

    def a(self, rho, T):
        T = np.asarray(T, float)
        return -0.5 * self.gamma_spec(rho) * T ** 2


class TitaniumColdModel(SumModel):
    """Vinet cold curve + Slater-Debye ions + Sommerfeld electrons."""

    name = "cold_model_ti"

    def __init__(self, csv_path=None, z_c=Z_COND):
        csv_path = csv_path or coldcurve_csv_path()
        P_GPa, V_mol = np.loadtxt(csv_path, delimiter=",", skiprows=4,
                                  unpack=True)
        v = V_mol / (M_TI * NA)  # cc/g
        (v0, B0, B0p), rms = fit_vinet(P_GPa * GPA_CGS, v)
        self.vinet = (v0, B0, B0p)
        self.vinet_rms = rms
        cold = VinetColdCurve(v0, B0, B0p)
        self.ion = SlaterDebye(M_TI, cold)
        self.electrons = SommerfeldElectrons(M_TI, z_c)
        super().__init__([cold, self.ion, self.electrons],
                         name=self.name)
        self.cold = cold

    def theta0(self):
        """Slater-Debye temperature at the ambient 0 K density (report)."""
        return float(self.ion.theta(1.0 / self.vinet[0]))
