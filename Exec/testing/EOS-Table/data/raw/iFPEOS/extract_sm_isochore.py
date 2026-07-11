#!/usr/bin/env python3
"""Re-runnable extraction of the single published iFPEOS isochore.

Source: Supplemental Material of Mihaylov et al., PRB 104, 144104 (2021)
(SM/Supplmental_Material_iFPEOS.pdf, sha256 in ../sources.yaml) — Sec. 2,
Table S1, rho = 0.001 g/cm^3, the ONLY isochore APS published (see
README.md addendum). Used as a validation overlay for the spliced table,
not as a splice source.

Row format in the PDF: 'T  P +- sP  E +- sE' with T in K, P in Mbar,
E in eV/atom; sigma = 0.0 marks the authors' interpolated points.

Writes ifpeos_rho0.001_isochore.csv (dimensional CGS, deuterium-specific
energy) with per-row provenance retained via the sigma columns.
"""

import os
import re

import pypdf

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "SM", "Supplmental_Material_iFPEOS.pdf")
OUT = os.path.join(HERE, "ifpeos_rho0.001_isochore.csv")

RHO = 0.001  # g/cc (deuterium basis: iFPEOS is a deuterium table)
N_EXPECTED = 39  # T points per the main text and our page-4 extraction

MBAR_CGS = 1.0e12          # Mbar -> erg/cc
EV_ERG = 1.602176634e-12   # eV -> erg
M_D = 2.01355321271 * 1.66053906660e-24  # g

ROW = re.compile(
    r"^\s*(\d+)\s+([0-9.eE+-]+)\s*±\s*([0-9.eE+-]+)\s+(-?[0-9.eE+-]+)\s*±"
    r"\s*([0-9.eE+-]+)\s*$")


def main():
    r = pypdf.PdfReader(PDF)
    text = r.pages[3].extract_text(extraction_mode="layout")
    rows = []
    for line in text.splitlines():
        m = ROW.match(line.replace("−", "-"))
        if m:
            rows.append(tuple(float(g) for g in m.groups()))
    assert len(rows) == N_EXPECTED, "expected %d rows, parsed %d" \
        % (N_EXPECTED, len(rows))
    T = [q[0] for q in rows]
    assert T == sorted(T) and T[0] == 800.0 and T[-1] == 256000000.0
    # P monotone nondecreasing in T along the isochore
    P = [q[1] for q in rows]
    assert all(b >= a * (1 - 1e-12) for a, b in zip(P, P[1:]))

    with open(OUT, "w") as f:
        f.write("# iFPEOS rho = 0.001 g/cc isochore — the only published "
                "isochore (PRB 104, 144104 SM Table S1)\n")
        f.write("# extracted by extract_sm_isochore.py; sigma = 0 marks "
                "authors' interpolated points\n")
        f.write("# validation overlay only (see README.md addendum); "
                "e converted at m_D = %.6e g\n" % M_D)
        f.write("T_K,P_erg_cc,sigmaP_erg_cc,e_erg_g,sigmaE_erg_g\n")
        for tK, pMbar, sP, eeV, sE in rows:
            f.write("%.6g,%.10e,%.3e,%.10e,%.3e\n" % (
                tK, pMbar * MBAR_CGS, sP * MBAR_CGS,
                eeV * EV_ERG / M_D, sE * EV_ERG / M_D))
    print("wrote %s (%d rows)" % (OUT, len(rows)))


if __name__ == "__main__":
    main()
