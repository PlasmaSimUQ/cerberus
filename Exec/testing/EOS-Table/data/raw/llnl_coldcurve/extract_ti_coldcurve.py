#!/usr/bin/env python3
"""Re-runnable extraction of the titanium T=0 isotherm from LLNL-JRNL-686936.

Same source/report as extract_h_coldcurve.py (sha256 in ../sources.yaml).
Titanium is the SECOND element column of the Table 2 pages headed
'Element Sc Ti V Cr Mn Fe Co Ni Cu Zn' (PDF pages 65-66); the report's Ti
isotherm is reduced from DAC compression data (report section 3.22).

Token order per page: ... 'Element', <10 element symbols — one of which is
the symbol 'V' (vanadium), so anchoring on the 'P' marker token, NOT on
'V' labels — 'P', 10 x 'V' column labels, nP pressure integers, then 10
element columns x nP floats (Sc first, Ti second).

Usage: python3 extract_ti_coldcurve.py  (writes ti_coldcurve_0K.csv here)
"""

import os
import re

import pypdf

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "LLNL-JRNL-686936_zeroK_isotherms.pdf")
OUT = os.path.join(HERE, "ti_coldcurve_0K.csv")

PAGES = [(64, 50), (65, 51)]  # PDF pages 65 (P=0..49) and 66 (P=50..100)
COLUMN = 1                    # 0-based: Sc=0, Ti=1


def parse_page(reader, pnum, n_p):
    toks = reader.pages[pnum].extract_text().split()
    ip = toks.index("P", toks.index("Element"))
    rest = toks[ip + 1:]
    assert rest[:10] == ["V"] * 10, "expected 10 column labels after 'P'"
    rest = rest[10:]
    P = [int(x) for x in rest[:n_p]]
    assert P == list(range(P[0], P[0] + n_p)), "P column not contiguous"
    vals = [float(x) for x in rest[n_p:] if re.match(r"^\d+\.\d+$", x)]
    assert len(vals) == 10 * n_p, (pnum, len(vals))
    return P, vals[COLUMN * n_p:(COLUMN + 1) * n_p]


def main():
    r = pypdf.PdfReader(PDF)
    P, V = [], []
    for pnum, n_p in PAGES:
        p_, v_ = parse_page(r, pnum, n_p)
        P += p_
        V += v_
    assert P == list(range(101)), "expected P = 0..100 GPa"
    assert all(b < a for a, b in zip(V, V[1:])), "V(P) must be monotone"
    # sanity: ambient Ti atomic volume ~ 10.6 cm^3/mol (rho ~ 4.5 g/cc)
    assert 10.0 < V[0] < 11.5, V[0]
    with open(OUT, "w") as f:
        f.write("# Titanium T=0 K compression isotherm, LLNL-JRNL-686936 "
                "Table 2 (Ti column)\n")
        f.write("# extracted by extract_ti_coldcurve.py from PDF pages "
                "65-66; see sources.yaml entry 'llnl-coldcurve'\n")
        f.write("# V is cm^3 per mole of atoms; DAC-reduced isotherm, "
                "report section 3.22\n")
        f.write("P_GPa,V_cm3_per_mol_atom\n")
        for p, v in zip(P, V):
            f.write("%d,%.3f\n" % (p, v))
    print("wrote %s (%d points); V0=%.3f cm3/mol -> rho0=%.4f g/cc"
          % (OUT, len(P), V[0], 47.867 / V[0]))


if __name__ == "__main__":
    main()
