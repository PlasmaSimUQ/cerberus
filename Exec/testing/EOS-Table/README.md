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
  PASS twice; **NOT bit-deterministic run-to-run** (all plotfile checksums
  differ under 8–10-rank MPI) → all later gates are verdict-identity. See
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
