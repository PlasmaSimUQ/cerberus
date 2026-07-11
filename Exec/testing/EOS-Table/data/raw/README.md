# Raw EOS source data

## fpeos_10-26-25.tar.gz

- **Source:** https://militzer.berkeley.edu/FPEOS/ (link `fpeos_10-26-25.tar.gz`,
  released 2025-10-26), retrieved 2026-07-07.
- **sha256:** `b42c44298e8630afe73e30db4755bd9a4b8a2f41ad147a82130accd3476a1ab1`; the archive is kept byte-identical to the
  download — all edits happen in `Exec/python_analysis/eos_table_prep.py`.
- **Reference:** B. Militzer, F. González-Cataldo, S. Zhang, K. P. Driver,
  F. Soubiran, "First-Principles Equation of State Database for Warm Dense
  Matter Computation", Phys. Rev. E 103 (2021) 013203.
- **License:** no explicit license file; the authors request citation of the
  six articles listed in each table's header (for H these are the
  Militzer/Ceperley/Hu deuterium PIMC papers, incl. Hu et al. PRB 84 224109
  (2011) — the original FPEOS deuterium table).
- Contents: per-material EOS text tables for 11 elements + 10 compounds,
  plus the authors' C++ interpolation/Hugoniot code (`FPEOS/*.h`, `*.C`)
  and Python plot scripts. Extract with `tar xzf` (extracted copy is not
  committed; only this archive is).

## FPEOS/H_EOS_09-18-20.txt — the table we condition (deuterium via mass scaling)

Characterisation (2026-07-07, feeds the W2 spec and regrid design):

- **401 (rho, T) points**: 33 densities from 0.000983 to 798.913 g/cc
  (~5.5 per decade), 17 distinct temperatures from 15,625 K to 6.4e7 K
  (1.35 eV to 5.5 keV), roughly factor-2 (log-uniform) spacing with
  irregular extras (95,250 K; 181,825 K; 500,000 K).
- **Ragged, not rectangular**: isochores carry 5-16 temperature points
  (9 isochores have all 16-17; low- and high-rho isochores as few as 5).
  The regrid MUST be hull-aware; the unsimulated corners are outside the
  source hull, not zero.
- **Line format** (whitespace tokens):
  `f= H N= 1 rho[g/cc]= <v> V[A^3]= <v> T[K]= <v> P[GPa]= <v> <err> E[Ha]= <v> <err>`
  i.e. pressure in GPa and internal energy in Hartree **per atom**, both
  with one-sigma error bars. No derivative columns — cv, dP/dT, dP/drho are
  produced offline by the conditioning tool (plan D4/D5).
- **Grid coarseness**: ~12 points per isochore over ~3.6 decades of T is
  coarse (the review's warning holds). The conditioning tool must measure
  regrid error (leave-one-out) and report it in the QA output; the
  "interpolation-level tolerance" used by Stage-3 gates inherits this
  number.
- **Deuterium mapping**: the H table is derived from the same PIMC/DFT-MD
  simulation campaign as the 2011 FPEOS *deuterium* table (its citation
  list says as much), expressed per nucleus. Isotope conversion used by the
  tool: rho_D = rho_H * (m_D/m_H) at equal nuclear number density and T;
  P(rho_D, T) = P_H(rho_H, T); specific energy e_D [erg/g] =
  E[Ha/atom] * E_Ha / m_D. With m_D = 2.013553 u, m_H = 1.007825 u:
  rho_D/rho_H = 1.99792. D-equivalent range: 0.00196-1596 g/cc — matching
  Hu et al. (2011) exactly.

## Hugoniot reference data

The archive's own tools compute Hugoniots from the tables (initial
conditions are recorded in each table header: for H, E0 = -0.5838 Ha,
V0 = 19.576 A^3/atom, i.e. rho0 = 0.0855 g/cc H ~ 0.171 g/cc D — liquid-D
initial density). The Stage-1 QA Hugoniot overlay therefore compares our
conditioned-table Hugoniot against the published FPEOS Hugoniot locus
(Hu et al. 2011, Fig. 13-15 points) and can cross-check against the
authors' own `fpeos` binary if built.

## Machine-checkable manifest

`sources.yaml` in this directory is the machine-checkable record of every
raw building block (citation, URL/DOI, retrieval date, sha256, license,
acquisition mode). Verify or fetch with:

```sh
python3 Exec/python_analysis/eos_table_prep.py sources --verify
python3 Exec/python_analysis/eos_table_prep.py sources --fetch
```

## iFPEOS (Mihaylov et al., PRB 104, 144104) — status 2026-07-11

The seam-1 WDM building block (see `doc/eos_creation_plan.md`).
Article PDF acquired manually (scripted fetch impossible: par.nsf.gov
unreachable, OSTI purl 500s, APS 403s to non-browser clients);
characterisation memo from the article text: `iFPEOS/README.md`
(53 ρ × 39 T points, 800 K–256 MK; seam-2 provisionally obsolete —
iFPEOS covers the whole FPEOS domain).
**Outstanding manual step:** the full table is the APS Supplemental
Material — download in a browser from
`https://link.aps.org/supplemental/10.1103/PhysRevB.104.144104` into
`iFPEOS/SM/`, add checksums to `sources.yaml`, then transcribe per the
plan's PDF-mitigation recipe.
