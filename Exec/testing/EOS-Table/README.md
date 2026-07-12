# EOS-Table — Stage 1: data groundwork for the tabulated EOS

Detailed plan for Stage 1 of `doc/eos_implementation_plan.md` (work items
W1–W3). This directory holds the Stage-1 artifacts: the baseline test-suite
record, raw source EOS data, the conditioned canonical table(s), and the QA
plots that gate the stage exit.

**Harness note:** this directory deliberately contains **no `check.py` yet** —
`run_tests.py` executes every subdirectory that has one, and Stage 1 is
Python-only (no MFP executable involved). The `check.py` + `run` pair arrives
at the end of Stage 2 with the one-zone self-test, which is hosted **in this
same directory** and consumes the table produced here — see `STAGE2.md`.

## Directory contract

```
EOS-Table/
  README.md        this plan
  baseline/        W1: recorded pass/fail state of the existing test suite
  data/raw/        W3: source tables as downloaded (FPEOS deuterium first)
  data/            W3: conditioned canonical output, e.g. D_fpeos.eostab
  qa/              W3: QA plots emitted by the conditioning tool
```

The conditioning tool itself lives in `Exec/python_analysis/eos_table_prep.py`
(shared tooling, not test data). Raw + conditioned tables are committed here
as long as they stay small (FPEOS-class tables are ~100s of kB); revisit if a
table exceeds a few MB.

---

## W1 — Baseline test-suite record  (effort: S)

Purpose: a regression anchor *before* any EOS-related source change lands.

1. Build the standard test executable (`Exec/local/`, default flags).
2. `cd Exec/testing && python3 run_tests.py` (writes `run_log.txt`).
3. Copy the result into `baseline/run_log_<date>_<git-sha>.txt` together with
   a one-line note of the build config (DIM, EB, PARTICLES, compiler, and
   whether `MFP_PRIM_FLOOR` is defined — floor behaviour is compile-gated, so
   every later comparison run must match this configuration).
4. Record which cases pass/fail *today* — pre-existing failures (if any) are
   part of the record, not something Stage 1 fixes.
5. **Determinism check** (review amendment): run the suite a second time from
   a clean tree and diff the two logs, plus a plotfile checksum on one
   single-rank case. This decides which later regression gates may demand
   bit-identity and which only verdict-identity (multi-rank MPI reduction
   order and per-case rebuilds make blanket bit-identity an unverified
   assumption). Record the verdict in `baseline/`.

Exit: committed baseline logs (both runs) + the determinism verdict. Every
later stage re-runs the suite and diffs against this record.

## W2 — Canonical `.eostab` format specification  (effort: S)

One self-describing ASCII format; the C++ reader (Stage 2) parses only this.
**Spec v1 — FROZEN 2026-07-07**, informed by the FPEOS source
characterisation (`data/raw/README.md`): the source is *ragged* (isochores
carry 5–16 T points), so the hull-mask block is **required**, and the
conditioned rectangular grid is clamped to the source hull with
nearest-hull-value fill outside it (marked 0 in the mask; runtime treats
those cells as clamp+flag territory from Stage 4 on).

```
EOSTAB 1                                # magic + format version
# --- provenance (free-form key: value lines, '#' comments allowed) ---
material:    D
source:      FPEOS (Hu et al., PRB 84 224109 (2011)); file <name>, retrieved <date>
generator:   eos_table_prep.py <git-sha>, run <date>
composition: A=2.014 Z=1
units:       cgs                        # values stored DIMENSIONAL;
                                        # Cerberus nondimensionalises at load (plan D7)
conditioning: cv_floor=<value> monotonised=dPdrho,cv maxwell=none
# --- grid (uniform in log10) ---
grid: n_rho=<Nr> n_T=<Nt>
lrho: <log10 rho_min> <log10 rho_max>
lT:   <log10 T_min>   <log10 T_max>
# optional inverse-map axes (required iff the T_of_e / T_of_p blocks below
# are present; each map lives on (lrho x its own log-uniform axis)):
le:   <log10 e_min> <log10 e_max>  n_e=<Ne>
lp:   <log10 p_min> <log10 p_max>  n_p=<Np>
# --- data blocks, each 'block: <name>' followed by Nr*Nt reals ---
# storage order: i_rho outer, j_T inner  (idx = i*n_T + j, matching the
# draft header MFP_eos_tab2T.H)
block: p          # pressure
block: e          # specific internal energy
block: dpdT       # (dP/dT)|rho      -- smoothed, for Riemann-contract outputs
block: dpdrho     # (dP/drho)|T      -- smoothed/monotonised
block: cv         # (de/dT)|rho      -- smoothed, cv>0 floored
block: dedrho     # (de/drho)|T      -- smoothed
# optional inverse-map blocks (Athena++-informed, plan D5/D8: runtime
# inversions become lookup + 1-2 Newton polish steps instead of cold Newton;
# computed offline by scanning the forward surface):
block: T_of_e     # T(rho, e)  -- seeds the re inversion (cons2prim)
block: T_of_p     # T(rho, p)  -- seeds the rp inversion + face evals (W10)
block: hull       # REQUIRED: 1.0 = inside source-data hull, 0.0 = filled
# other optional blocks: zbar, entropy
```

Decisions encoded: ASCII for inspectability/diff-ability (tables are small);
axes given by log-uniform range + count (reader reconstructs, no stored axis
arrays to disagree with); derivatives are *data*, pre-conditioned offline
(plan D4/D5); dimensional CGS storage with nondimensionalisation at load
(plan D7); inverse maps are *data* too, computed offline (plan D5 —
Athena++'s pre-inverted-table pattern, hybridised: our forward (ρ,T) surface
stays canonical and the maps only seed the Newton polish). Athena++'s
dimensionless-ratio storage (log p/e etc., exact for ideal gas) is a noted
possible refinement, not in v1.

Exit: spec section above frozen and mirrored in a docstring in the writer.

## W3 — Offline conditioning tool + FPEOS deuterium ingest  (effort: M)

Tool: `Exec/python_analysis/eos_table_prep.py` (numpy/scipy/matplotlib; add to
`Exec/python_analysis/environment.yml`). Pipeline:

1. **Ingest** — readers for (a) FPEOS published ASCII layout, (b) a generic
   `raw` reader for whitespace columns `rho T P e`. Each returns points +
   metadata.
2. **Regrid** — resample onto the uniform (log10 rho, log10 T) target grid.
   PCHIP (monotone piecewise-cubic Hermite) along isotherms then isochores;
   refuse to extrapolate beyond the source hull (clamp grid to hull instead).
   If the source grid proves coarse (first-principles PIMC/DFT tables can be
   ~10² points, not 10⁴), regrid error may dominate everything the
   conditioning controls — measure it (leave-one-out residuals on source
   points) and record the number in the QA output.
3. **Condition** —
   - cv = de/dT via PCHIP slopes along isochores; floor cv >= cv_floor > 0,
     record where flooring triggered;
   - (dP/drho)|T via PCHIP slopes along isotherms; monotonise (never negative)
     and record where;
   - Maxwell construction: **verify first whether FPEOS needs it** (PIMC/DFT
     tables are usually loop-free); implement only if the QA isotherm plots
     show van-der-Waals loops. Record `maxwell:` status in the header either way.
4. **Write** — emit `.eostab` per W2 spec into `data/`, including the
   inverse-map blocks `T_of_e`/`T_of_p` (computed by inverting the
   conditioned forward surface column-by-column in rho).
5. **QA plots** (into `qa/`, all gated on visually sane output):
   - isotherms P(rho) and e(rho) at ~10 temperatures spanning the hull
     (log-log), raw points overplotted on the conditioned grid;
   - cv(rho,T) heatmap with floored cells marked;
   - hull/coverage map (source points vs target grid);
   - **principal Hugoniot computed from the conditioned table** (solve the
     Rankine–Hugoniot energy relation e−e0 = ½(P+P0)(1/rho0 − 1/rho) per
     initial state) overlaid on the published FPEOS Hugoniot points — this is
     the physics gate, and the same curve Stage 4 validates in-code;
   - numerical round trip: invert e(rho,T) for T over a grid sweep in Python
     (prototype of the Stage-2 C++ inversion) and plot residuals.

Data acquisition (TBD — first concrete task of W3):
- FPEOS deuterium: published as supplementary material of Hu, Militzer et al.
  (PRB 84, 224109 (2011); iFPEOS update 2021). Locate the downloadable table,
  note its license/citation requirements in `data/raw/README` alongside the
  file, and keep the file byte-identical to the download (all edits happen in
  the tool).
- Fallback if download access stalls: generate a synthetic ideal-gas +
  analytic-correction table with the `raw` reader to unblock Stage 2 C++ work;
  swap in FPEOS when available (format identical, so nothing downstream moves).

Exit (= Stage 1 exit): conditioned `D_fpeos.eostab` committed in `data/`,
QA plots committed in `qa/`, Hugoniot overlay agrees with published points
(eyeball first; commit a numeric threshold into the QA script the first time
it passes, so later regenerations are gated mechanically — review amendment),
generating command recorded in this README.

## Stage-1 results (2026-07-07)

Stage 1 executed; all three work items closed. Generating commands (from
this directory):

```sh
python3 ../../python_analysis/eos_table_prep.py synthetic \
    --out data/ideal_synthetic.eostab --qa qa
python3 ../../python_analysis/eos_table_prep.py fpeos \
    --src data/raw/FPEOS/H_EOS_09-18-20.txt \
    --out data/D_fpeos.eostab --qa qa \
    --hug-ref data/raw/hugoniot_MC2000_PRL85_1890.txt
```

(`data/raw/FPEOS/` is the untracked extraction of the committed archive —
`tar xzf fpeos_10-26-25.tar.gz` in `data/raw/` first.)

- **W1**: 3-case suite (Couette, Double-Rarefaction, Viscous-Vortex), all
  PASS twice. Determinism (corrected 2026-07-08): whole-file checksums
  differ run-to-run, but `fcompare` shows **all physical fields
  bit-identical** — only the `cost` load-balance diagnostic and
  rank-to-file packing differ. Gates: verdict-identity in bulk +
  `fcompare` field-identity (excluding `cost`) as the spot-check. See
  `baseline/BASELINE.md`.
- **W2**: spec v1 frozen above; mirrored in the `eos_table_prep.py`
  docstring.
- **W3**: `data/D_fpeos.eostab` (96×96, deuterium via isotope scaling of
  the FPEOS-2021 H table): hull coverage 76.9% (source is ragged); cv
  flooring hit 51/7083 hull cells (0.7%; all 2133 filled cells floored by
  construction); `maxwell: none-needed` (zero dP/drho<0 pre-monotonise);
  leave-one-out regrid error P median 0.23% / max 10.5%, e (span-normed)
  median 0.08% / max 4.3%.
- **Hugoniot gate**: table-computed locus vs published PIMC points
  (Militzer–Ceperley PRL 85 1890, same rho0=0.171 g/cc): +0.08…+0.26%
  compression error for T ≥ 125 kK (pure-PIMC in both datasets);
  +1.0…+4.4% at 31–62 kK, the region FPEOS-2011/2021 deliberately revised
  with DFT-MD (peak compression ≈4.5 vs MC2000's ≈4.29) — a dataset
  difference, not a conditioning error. Locus points below ~94 GPa fall
  under the table's T floor (1.35 eV) by construction.
- **Committed numeric threshold** (replacing "eyeball", per review
  amendment): regenerations must reproduce the published compressions to
  **≤0.5% for the T ≥ 125 kK points** (and the ≤5% low-T deviation is
  expected — investigate if it *shrinks*, since that would mean the
  conditioning is reverting toward MC2000).

## Stage-1 execution order

1. **W3.data first** (review amendment — this was the plan's least-verified
   external dependency): locate/download FPEOS deuterium and characterise it
   — units (eV vs K, GPa vs erg/cc), grid density, whether derivatives ship,
   license. Its shape feeds both the W2 spec and the regrid design, and it
   has the longest external latency. Kick off the W1 suite run in parallel
   (it is slow).
2. W2 spec freeze (fast; now informed by the real source-data shape). The
   inverse-map blocks are part of the freeze — deciding them later would
   re-open the format.
3. W1 second (determinism) run + baseline record.
4. W3.tool: ingest → regrid → condition → write → QA (bulk of the work;
   synthetic-table fallback keeps this moving if 1 stalls).
5. Review QA plots together before declaring Stage 1 done.

## SS1 record (2026-07-11) — offline-tool refactor + harvest infrastructure

Plan: `doc/eos_creation_plan.md`. `eos_table_prep.py` became a CLI shim
over the new `Exec/python_analysis/eos_tools/` package (code moved verbatim;
the RH locus solver factored into `eos_tools/hugoniot.py`); raw-source
manifest `data/raw/sources.yaml` + `sources --verify/--fetch` subcommand;
pytest unit tests under `Exec/python_analysis/tests/` gated by the new
`Exec/testing/EOS-PyTools` adapter case.

Measured gates (all green):
- regenerated `ideal_synthetic.eostab` and `D_fpeos.eostab` byte-identical
  to the committed tables modulo the `generator:` provenance line;
- pytest: 15 passed (regrid e-interpolation tolerance committed at the
  measured 4.1% max for the 12-points-per-isochore synthetic source);
- `sources --verify`: 0 failures (fpeos, hugoniot-mc2000);
- affected suite verdict-identical: `EOS-Table` PASS (self-test, D_fpeos
  fd-vs-blocks reported-not-gated as before), `EOS-Sod-Ideal` PASS
  (twin/roundoff/analytic/conservation/wall gates; geos roundoff 5e-15,
  wall-geos 1.48x), `EOS-Hugoniot` PASS (21 gates incl. clamps-fired
  15598 > 0). No C++ or committed-table changes, so the remaining cases
  run identical binaries and inputs.

Open item: iFPEOS (PRB 104 144104) acquisition is a **manual** step from
this network — see `data/raw/README.md` and the `ifpeos` manifest entry.

## SS3 record (2026-07-12) — spliced wide-range deuterium table

Artifact: `data/D_spliced.eostab.gz` (sha256
`72a6cd22b4e78407b8c490a265d14ce6fd3490baf67cf9c12e03a36f4c79f7f5`,
6.5 MB; decompressed 22.6 MB, gitignored — `gunzip -k` in run scripts).
384x384, lrho [-4,3], lT [1.25,9]. Pipeline: `eos_tools/splice.py`
(`splice --material D`); sources cold-composite | H-REOS.3(D) | FPEOS(D) |
ideal-plasma; seams lT = 2.55/5.5/7.6 (T-only tanh, hw 0.15/0.25/0.20);
rho hand-off to REOS.3 at 0.24 g/cc; energy chain-aligned (constant cE per
pair); tension clip at 1e3 barye (0 cells fired); monotone-in-T value
surfaces enforced (monoT_p=5887 cells max 3.3 rel in hull-0 fills,
monoT_e=3596 max 0.9).

Measured gates, all green (SPLICE-QA):
- alignment constancy |dev|/kT: 0.0090 / 0.0154 / 0.0107 (gate < 0.10);
- convexity: cs^2 > 0 at all 126,725 in-hull cells (85.9% hull);
- Maxwell residual: blended median 3.1e-3 vs REOS.3 yardstick 3.8e-3
  (gate < 2x source baseline);
- C1 seams (in-band p99 vs 3x out-of-band p99): 3.1e-3 / 2.3e-3 / 3.1e-4
  vs 1.7e-2;
- cold-start principal Hugoniot from table (0.171 g/cc, 20 K): peak
  compression 4.398 in [4.2, 4.9]; max adjacent stride 0.136 for
  P >= 0.3 GPa (gate 0.15); MC2000 PIMC anchor median 1.19% over 8
  in-range points (gate 5%);
- resolution convergence 384 vs 768: Hugoniot compression shift median
  0.0028% (gate 0.2%), max 0.37%;
- iFPEOS rho=0.001 isochore overlay (validation): median 1.92% over 35
  sigma-filtered points;
- C++ one-zone self-test (DEBUG exe) on the emitted file: reader,
  roundtrip-e (5.4e-11), roundtrip-p (9.9e-11), identities (4.6e-16),
  hull, corner (7.3e-11) all PASS; fd-vs-blocks 49.1 reported-not-gated
  (tier-2 policy; localized to fill-plateau cells where the secant is
  ~0/0).

Design decisions forced by measurement (full log in
doc/eos_creation_plan.md §4):
- seam 1 moved 3.2 -> 2.55 (355 K): the cold model's internal
  CoolProp->rotor ramp (505-600 K) put a 24% p(T) dip at liquid densities;
  REOS.3's 60 K floor makes the rotor piece unnecessary as a table source;
- rho hand-off 1.0 -> 0.24 g/cc: QEOS-solid-above-melt vs REOS.3 mismatch
  at (0.4 g/cc, 1-2.5 kK) put a 0.43 compression discontinuity on the
  Hugoniot at ~10 GPa; a T-sliding hand-off over-corrected (diagonal kink)
  and was reverted;
- fallback ladder (never data-source fills): the FPEOS nearest-fill below
  its 2e-3 g/cc floor planted non-monotone p/e bumps across dilute rows —
  caught by the C++ round-trip self-test, not the offline QA;
- monotonise-in-T conditioning pass added (the v1 single-branch inversion
  contract); satL fill extended below CoolProp's 18.72 K floor.

Known limitations (v1): melt smeared (no latent heat); REOS.3 classical
ions (no NQE) and PBE-class XC; dome + sub-60 K/high-rho corners are
hull-0 fills; experimental Hugoniot compilations (Nellis/Knudson/Hicks/
Boriskov) still pending manual harvest — validation currently rests on
the MC2000 PIMC anchor + iFPEOS isochore.

## SS4 record (2026-07-12) — cold-start shock validation (EOS-ColdShock)

New case `Exec/testing/EOS-ColdShock/` (EOS-Hugoniot structure): a hot
driver slab shocks cryogenic liquid deuterium on `data/D_spliced.eostab`
(gunzipped by `run`). Initial state (0.171 g/cc, 24 K, ~90 bar): 24 K not
~20 K because the (0.171, T) cell is fully in-hull only for T >= ~23.5 K
(dome edge satL(20 K) = 0.1718 > 0.171, plus the hull-flag shoulder of the
splice rho hand-off at 0.185-0.23 g/cc below 60 K) — the sanctioned
un-ionized numerical-relief start (shock-notes section 5).

Measured gates, all 24 green (check.py):
- three drive strengths from ONE cold start: gas-gun 6.67 GPa
  (compression 2.562, RH-locus err 0.20%), multi-Mbar 176.2 GPa (4.156,
  1.46%), plasma 6.26 TPa (4.064, 1.32%) — solid -> dissociating fluid ->
  plasma emerges from the conservative update alone;
- pre-shock closure: drho <= 2e-13, dT <= 9e-12, de = 3.48e-7 (identical
  in all three runs — nondim unit-chain round-off, committed at 1e-6),
  |p_sim - p_tab|/p_plateau <= 3e-10 (the reference cell's own p is
  resolution-sensitive; plateau-normalized per the plan's recorded trap);
- **noclamp gates: 0 hull clamps in all three physical runs** — the whole
  24 K -> plasma path lives inside the table hull (the SS4-specific gate);
- anchor-PIMC (from the table itself, no offline artifact): median 1.19%
  over 8 in-range MC2000 points;
- conservation drift <= 2.2e-15; boundary-flux-predicted mass loss on the
  abusive run to 1.3e-6;
- abusive run = VACUUM-FORMING expansion (|u| = 2 code ~ 40x the liquid
  cs): the wide-range table's 17.8 K floor is hydrodynamically
  unreachable (a 5x-cs dome expansion completed with ZERO clamps —
  measured), so density-below-table-edge is the only abusive regime
  left. It completes in ~2 s wall with 769,777 tallied clamps and finite
  fields — the regime the FPEOS-table case deferred to W16 (absolute
  floors collapsed dt) is closed by the hull-physical floors (W8.1) on
  this table.

Suite re-run after SS4: EOS-PyTools, EOS-Table, EOS-Sod-Ideal,
EOS-Hugoniot all PASS (recorded in doc/eos_creation_plan.md §4).

## SS5b record, part 1 (2026-07-12) — titanium solid model (Ti SS2')

Ti pipeline started on the proven machinery (doc/eos_creation_plan.md §2.2).
`eos_tools/materials/titanium.py`: LLNL Ti 0 K isotherm (Table 2 Ti column,
`data/raw/llnl_coldcurve/ti_coldcurve_0K.csv`, re-runnable extraction) +
Vinet fit + Slater-Debye ions (shared `models/qeos.SlaterDebye`) + a
Sommerfeld electronic term (free-electron FD gamma at z_c = 4; d-band
enhancement ~3x documented as the known v1 limitation).

Measured QA (`coldmodel --material Ti`), all gates green:
- Vinet fit rms 0.05% over 0-100 GPa; **B0 = 109.4 GPa, B0' = 3.65 —
  inside the DAC literature window without being an input**;
  rho0(fit) = 4.5802 vs the extracted 4.580 g/cc datum;
- Slater theta0 = 562 K vs literature Debye 420 K (the expected
  shear-blind Slater overestimate; report-only);
- gamma_e(free-electron, z_c=4) = 1.03 mJ/mol/K^2 vs literature ~3.5
  (the d-band factor; report-only);
- cs(4.51 g/cc, cold) = 4.79 km/s vs Ti bulk sound speed ~4.9-5.2
  (report-only);
- cv > 0, cs^2 > 0 at all p > 0 cells of the (3-25 g/cc) x (17.8 K-4 kK)
  rectangle; Maxwell residual median 1.5e-4.

Harvested for the Ti hot side: NIST ASD Ti I-XXII ionization energies
(`data/raw/nist_ti/`, scripted) — input for the planned Saha
average-ionization model. The ML-MD melt constraints (arXiv 2603.04680)
ship no data files; figure digitization or an author request is the
remaining manual step. Ti SS3' (Saha hot side + splice + Ti_spliced
.eostab) is the next stage.

## SS5b record, part 2 (2026-07-12) — Saha model + Ti_spliced.eostab (Ti SS3')

`eos_tools/models/saha.py`: chemical-picture average-ionization model on
the NIST-harvested Ti I-XXII energies (Zbar by bracketed root-finding,
log-space populations; U_j = 1, Boltzmann electrons, no continuum
lowering — all recorded). Unit tests: hydrogenic half-ionization
benchmark, dilute/full-ionization limits, Ti Zbar staging, and a
Maxwell-residual consistency check of the equilibrium (P, e) pair
(median 8e-4 at production grid pitch — envelope-theorem consistency
verified numerically).

`eos_tools/splice_ti.py` -> `data/Ti_spliced.eostab.gz` (sha256
d02ac84b10c691c3444130d9c6872af2b200d7543116e97b2223b6bef882bbe5),
384x384, 1e-3..50 g/cc, 17.8 K..1e9 K: solid | Saha | FD ideal plasma
(Z=22), seams lT = 4.50/8.40, rho hand-off 3.2 g/cc, expanded-condensed
wedge (0.02-3.2 g/cc below 4 kK) and tension clip (3311 cells) flagged
hull-0.

Measured gates (structural gates hard, WDM-band metrics reports — the
**Ti v0 accuracy statement**: the seam-1 band is uncontrolled pending
OFMD/ML-MD data, exactly the splice-plan §4 "genuinely new work" gap):
- convexity: cs^2 > 0 at all 125,115 in-hull cells (84.8% hull); finite;
- Hugoniot from ambient solid (4.51 g/cc, 293 K): max compression stride
  0.139 for P >= 1 GPa (grid-striding-verified: 0.524 at 96^2); peak
  compression 5.589 (report — no Ti reference data harvested yet);
- seam-3 Saha vs FD ideal plasma: P agreement median 0.000%, max 0.033%
  (a genuine cross-validation of the two independent hot models);
- alignment: ip->saha |dev|/kT = 0.001; saha->solid |dev|/kT = 4.8
  (report — the cohesion/WDM gap, by construction);
- C++ one-zone self-test: reader/roundtrip-e (9.6e-11)/roundtrip-p
  (8.9e-11)/identities (4.2e-16)/hull/corner (9.3e-11) all PASS;
  fd-vs-blocks reported-not-gated (metric degenerates on the larger Ti
  fill plateaus).

Path to Ti v1 (recorded): OFMD or author-supplied WDM data for the seam-1
band ('ti-mlmd-melt' manifest entry, manual); Ti experimental Hugoniot
compilation for the anchor gate; optional degeneracy/continuum-lowering
upgrades to the Saha side.
