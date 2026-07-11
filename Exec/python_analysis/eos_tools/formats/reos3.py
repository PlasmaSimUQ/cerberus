"""H-REOS.3 ingest (Becker et al., ApJS 215, 21 (2014); VizieR J/ApJS/215/21).

table2.dat columns (whitespace separated, one row per point):
    rho [g/cc]   T [K]   P [GPa]   u [kJ/g]
organized as isotherms (rho inner). 4410 rows = 42 isotherms x 105 points.

Isotope scaling H -> D (exact for classical ions, which is what REOS.3's
DFT-MD uses — the lost NQE difference is a documented caveat): at equal
nuclear number density and T,
    rho_D = rho_H * (m_D/m_H),   P unchanged,   e_D = e_H * (m_H/m_D).

Reader invariants (transcription/download guards): expected row count,
strictly positive P, isotherm structure, P monotone nondecreasing in rho
along every isotherm.
"""

import numpy as np

from ..constants import GPA_CGS, M_D, M_H

N_ROWS = 4410
KJ_G = 1.0e10  # kJ/g -> erg/g


def read_reos3(path):
    """Parse table2.dat -> dict of 1-D arrays in CGS, HYDROGEN basis."""
    rho, T, P, u = np.loadtxt(path, unpack=True)
    if len(rho) != N_ROWS:
        raise ValueError("%s: expected %d rows, got %d"
                         % (path, N_ROWS, len(rho)))
    if not np.all(P > 0.0):
        raise ValueError("%s: non-positive pressures" % path)
    uT = np.unique(T)
    for t in uT:
        m = T == t
        o = np.argsort(rho[m])
        dP = np.diff(P[m][o])
        if np.any(dP < -1e-12 * P[m].max()):
            raise ValueError("%s: P not monotone in rho on the %g K "
                             "isotherm" % (path, t))
    return {
        "rho": rho, "T": T,
        "p": P * GPA_CGS,       # erg/cc
        "e": u * KJ_G,          # erg/g
        "n_isotherms": len(uT),
    }


def reos3_to_deuterium(raw):
    """Classical isotope mass scaling, hydrogen -> deuterium (CGS)."""
    s = M_D / M_H
    return {
        "rho": raw["rho"] * s,
        "T": raw["T"].copy(),
        "p": raw["p"].copy(),
        "e": raw["e"] / s,
    }
