# iFPEOS — characterisation memo (SS1, 2026-07-11)

**Source:** D. I. Mihaylov, V. V. Karasiev, S. X. Hu, J. R. Rygg,
V. N. Goncharov, G. W. Collins, *"Improved first-principles
equation-of-state table of deuterium for high-energy-density
applications"*, Phys. Rev. B **104**, 144104 (2021).
`iFPEOS_Deuterium_PhysRevB.104.144104.pdf` in this directory is the
accepted manuscript (sha256 in `../sources.yaml`), downloaded manually —
NSF-PAR, OSTI and APS all refuse scripted clients from this network.

**Status: the article PDF contains NO table data.** The full table is the
APS Supplemental Material (paper ref [36],
`https://link.aps.org/supplemental/10.1103/PhysRevB.104.144104`) — a
second **manual browser download** into `SM/` here, then checksum +
transcription per the PDF-mitigation recipe
(`doc/eos_creation_plan.md` Part 2). Everything below is
characterisation extracted from the article text, pending confirmation
against the SM data.

## Grid and methods (paper §III, Fig. 2)

- **53 ρ points, 0.001 ⩽ ρ ⩽ 1596.49 g/cm³; 39 T points, 800 K ⩽ T ⩽
  256 MK** (deuterium mass basis — no isotope scaling needed, unlike the
  FPEOS H table). Nominally rectangular (Fig. 2 shows every (ρ,T) point,
  labelled by method); raggedness/error bars to be confirmed from the SM.
- **Method per region:** KSMD (VASP, PAW, T-SCAN-L(+rVV10) XC) in the
  strongly-coupled/degenerate region; OFMD (PROFESS@QE, LKTFγTF with a
  Padé-fit γ(ρ)) at high T; **the corner 0.002 ⩽ ρ ⩽ 0.084 g/cm³,
  800 K ⩽ T ⩽ 182 kK is the authors' interpolation**, built from the
  ρ = 0.001, 0.1, 0.2, 0.3 g/cm³ isochores and the T = 182/250/400 kK
  rows — treat as lower-confidence in blend weights/QA (author-filled,
  analogous to our filled-hull cells).
- **Internal energy consistency is already handled by the authors**: OFMD
  energies are shifted to match KSMD at the switch temperature per
  isochore (weak T dependence reported), so the published table is a
  single consistent energy surface — our splice-plan §5.2 KS↔OF shifts do
  NOT need re-doing for iFPEOS.
- **NQEs via PIMD** applied from 800 K upward at select densities
  (~0.1–1.49 g/cm³) until they vanish; below/above that ρ range the
  correction is extrapolated/neglected by the authors.

## Implications for the splice plan (seam decisions)

- **Seam 1 (cold model ↔ iFPEOS, band 1000–2500 K)**: sound. The D₂
  reference density 0.171 g/cm³ sits in the *simulated* (KSMD+PIMD)
  region, not the interpolated corner (which ends at 0.084 g/cm³). The
  800 K floor leaves margin below the band bottom.
- **Seam 2 (iFPEOS ↔ FPEOS at ~200 kK)**: **provisionally obsolete.**
  iFPEOS spans the entire FPEOS domain (FPEOS: 0.002–1596 g/cc
  D-equivalent, 15.6 kK–64 MK) and reaches 4× hotter (256 MK). Unless the
  SM data shows worse raggedness/error bars than FPEOS in the plasma
  band, FPEOS reduces to a **validation overlay** (per the SS1
  contingency in `doc/eos_creation_plan.md` §1.3), and the deuterium
  table becomes a 3-source splice: cold model ↔ iFPEOS ↔ ideal-plasma.
- **Seam 3 (↔ ideal-plasma model)**: can move up from lT ≈ 7.6 toward
  lT ≈ 8.2–8.3 (iFPEOS tops out at 2.56×10⁸ K); the model still owns
  2.56×10⁸–10⁹ K and stays the Tier-2 entropy anchor.

Final seam decisions are taken after the SM table is transcribed and its
per-point coverage/error columns are inspected.

## Addendum (2026-07-12): the full table is unpublished

The APS Supplemental Material was confirmed to consist of the single
4-page PDF in `SM/` — Sec. 2 contains only the ρ = 0.001 g/cm³
demonstration isochore (39 T rows, T (K) / P (Mbar)±σ / E (eV/atom)±σ,
σ = 0 marking interpolated points; it also corrects the
interpolation-region isochores to 0.001/0.1/0.3/0.38 g/cm³). The full
53-isochore table was never published.

Consequences (see `doc/eos_creation_plan.md` §1.3 v2): iFPEOS is demoted
to **validation data** (this isochore + the paper's Hugoniot/dissociation
figures); **H-REOS.3** (`../reos3/`, Becker et al. ApJS 215, 21) takes the
seam-1/WDM role, and the seam-2 analysis above is superseded — seam 2 is
back, now REOS.3 ↔ FPEOS at lT ≈ 5.5. Upgrade path if desired: request
the full table from the LLE authors (Mihaylov/Karasiev/Hu).
