"""Ideal molecular D2 gas: translation + classical rigid rotor + harmonic
vibration. Covers the dilute-gas window between CoolProp's 600 K ceiling
and the seam-1 band (~2500 K), where D2 is still molecular.

    a = a_trans(ST, mass m_D2) + a_rot + a_vib
    a_rot = -R2 T ln(T / (sigma theta_rot)),  sigma = 2 (homonuclear)
    a_vib =  R2 T ln(1 - exp(-theta_vib / T))

R2 = kB / m_D2 (per gram of D2). Constants derive from the H2 spectroscopic
values by reduced-mass scaling (factor 1/2 for B_e, 1/sqrt(2) for omega_e):
theta_rot(D2) = 85.3/2 K ~ 43 K wrt H2's 85.3 K; theta_vib(D2) =
6332/sqrt(2) ~ 4477 K wrt H2's 6332 K (omega_e = 4401 cm^-1). Nuclear-spin
entropy is EXCLUDED (practical-absolute convention; constant offsets are
absorbed by the splice energy alignment anyway). Classical rotor requires
T >> theta_rot — satisfied wherever this piece carries weight (T > ~500 K).
Dissociation is deliberately absent: iFPEOS owns it above the seam-1 band.
"""

import numpy as np

from ..constants import H_PLANCK, KB, M_D
from .base import HelmholtzModel

M_D2 = 2.0 * M_D
THETA_ROT_D2 = 85.3 / 2.0     # K  (B_e(H2) = 60.85 cm^-1 -> 85.3 K, /2 for D2)
THETA_VIB_D2 = 6332.0 / np.sqrt(2.0)  # K (omega_e(H2) = 4401 cm^-1 -> 6332 K)
SIGMA = 2.0


class RotorGasD2(HelmholtzModel):
    name = "rotor_gas_d2"

    def __init__(self, include_vib=True):
        self.include_vib = include_vib
        self.R2 = KB / M_D2

    def a(self, rho, T):
        rho = np.asarray(rho, float)
        T = np.asarray(T, float)
        n = rho / M_D2  # molecules per cc
        lam3 = (H_PLANCK ** 2 / (2.0 * np.pi * M_D2 * KB * T)) ** 1.5
        a_tr = -self.R2 * T * (np.log(1.0 / (n * lam3)) + 1.0)
        a_rot = -self.R2 * T * np.log(T / (SIGMA * THETA_ROT_D2))
        out = a_tr + a_rot
        if self.include_vib:
            out = out + self.R2 * T * np.log(-np.expm1(-THETA_VIB_D2 / T))
        return out
