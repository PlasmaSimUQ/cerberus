#!/usr/bin/env python3
"""Re-runnable extraction of the hydrogen T=0 isotherm from LLNL-JRNL-686936.

Source: D. A. Young, H. Cynn, P. Soderlind, A. Lanza, "Zero-Kelvin
Compression Isotherms of the Elements 1 <= Z <= 92 to 100 GPa",
J. Phys. Chem. Ref. Data 45, 043101 (2016); LLNL-JRNL-686936 preprint,
https://www.osti.gov/pages/servlets/purl/1342023.

Table 2 tabulates atomic volumes V (cm^3 per mole of ATOMS) at T = 0 K in
1-GPa steps, 0..100 GPa. Hydrogen is the first element column on the two
pages whose header starts 'Element H He ...' (PDF pages 61-62). Their
hydrogen isotherm is built from Loubeyre et al. DAC data (Vinet fit) merged
with an exp-6 model at low pressure (report section 3.1) and includes the
H-mass zero-point contribution (see H_VS_D_CAVEAT in the emitted CSV
header).

Token order on those pages after the 10 'V' header tokens: the nP pressure
integers, then 10 element columns x nP floats each, hydrogen first
(verified: 500 = 10x50 and 510 = 10x51 volume tokens exactly).

Usage: python3 extract_h_coldcurve.py  (writes h_coldcurve_0K.csv here)
"""

import os
import re

import pypdf

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "LLNL-JRNL-686936_zeroK_isotherms.pdf")
OUT = os.path.join(HERE, "h_coldcurve_0K.csv")

# (0-based page index, number of pressure rows)
PAGES = [(60, 50), (61, 51)]  # PDF pages 61 (P=0..49) and 62 (P=50..100)


def parse_page(reader, pnum, n_p):
    toks = reader.pages[pnum].extract_text().split()
    idx = max(i for i, x in enumerate(toks) if x == "V")
    rest = toks[idx + 1:]
    P = [int(x) for x in rest[:n_p]]
    vals = [float(x) for x in rest[n_p:] if re.match(r"^\d+\.\d+$", x)]
    assert len(vals) == 10 * n_p, (pnum, len(vals))
    assert P == list(range(P[0], P[0] + n_p)), "P column not contiguous"
    return P, vals[:n_p]  # hydrogen is the first element column


def main():
    r = pypdf.PdfReader(PDF)
    P, V = [], []
    for pnum, n_p in PAGES:
        p, v = parse_page(r, pnum, n_p)
        P += p
        V += v
    assert P == list(range(101)), "expected P = 0..100 GPa"
    assert all(b < a for a, b in zip(V, V[1:])), "V(P) must be monotone"
    with open(OUT, "w") as f:
        f.write("# Hydrogen T=0 K compression isotherm, LLNL-JRNL-686936 "
                "Table 2 (H column)\n")
        f.write("# extracted by extract_h_coldcurve.py from PDF pages 61-62; "
                "see sources.yaml entry 'llnl-coldcurve'\n")
        f.write("# H_VS_D_CAVEAT: column is for H (H-mass zero-point "
                "included by the report's lattice reduction); used for D2 "
                "on the strength of Loubeyre 1996's finding that the "
                "H2/D2 EOS difference is small; revisit if the cold-model "
                "QA shows a mismatch against the fluid-D2 boundary.\n")
        f.write("P_GPa,V_cm3_per_mol_atom\n")
        for p, v in zip(P, V):
            f.write("%d,%.3f\n" % (p, v))
    print("wrote %s (%d points)" % (OUT, len(P)))


if __name__ == "__main__":
    main()
