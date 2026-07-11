"""FPEOS ingest (militzer.berkeley.edu first-principles EOS database)."""

import numpy as np

from ..constants import GPA_CGS, HA_ERG, M_D, M_H


def read_fpeos(path):
    """Parse an FPEOS *_EOS_*.txt element table.

    Line format:
    f= H N= 1 rho[g/cc]= v V[A^3]= v T[K]= v P[GPa]= v err E[Ha]= v err
    Returns dict of 1-D arrays (per-atom energies, hydrogen-mass densities)
    and the Hugoniot initial condition from the header.
    """
    rho, T, P, E = [], [], [], []
    e0_ha, v0_a3 = None, None
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                if "E0[Ha]=" in line:
                    e0_ha = float(line.split("E0[Ha]=")[1].split()[0])
                if "V0[A^3]=" in line:
                    v0_a3 = float(line.split("V0[A^3]=")[1].split()[0])
                continue
            t = line.split()
            if len(t) < 16 or t[0] != "f=":
                continue
            rho.append(float(t[5]))
            T.append(float(t[9]))
            P.append(float(t[11]))
            E.append(float(t[14]))
    return {
        "rho": np.array(rho), "T": np.array(T),
        "P_GPa": np.array(P), "E_Ha": np.array(E),
        "E0_Ha": e0_ha, "V0_A3": v0_a3, "n_atoms_per_formula": 1,
    }


def fpeos_to_deuterium(raw):
    """Isotope-scale the per-atom hydrogen table to deuterium (CGS).

    Equal nuclear number density and T: rho_D = rho_H * m_D/m_H,
    P unchanged, e[erg/g] = E[Ha/atom]*HA_ERG/m_D.
    """
    s = M_D / M_H
    pts = {
        "rho": raw["rho"] * s,
        "T": raw["T"].copy(),
        "p": raw["P_GPa"] * GPA_CGS,
        "e": raw["E_Ha"] * HA_ERG / M_D,
    }
    # Hugoniot initial condition, same scaling (V0 is per atom)
    pts["hug_rho0"] = M_D / (raw["V0_A3"] * 1e-24)
    pts["hug_e0"] = raw["E0_Ha"] * HA_ERG / M_D
    pts["hug_p0"] = 0.0
    return pts
