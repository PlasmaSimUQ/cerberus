# Offline EOS-creation toolchain + building-block harvesting — Plan & progress

**Living document** — the plan for producing wide-regime `.eostab` tables
(cold solid/liquid → WDM → classical plasma) for single-fluid hydro, and the
harvesting of their building-block data. Updated as stages land; the
progress log is §4. Approved 2026-07-11.

Companion documents:
- `doc/eos_splice_plan.md` — data-source & splice strategy (Route B implemented here; Route A/SESAME **out of scope** by decision)
- `doc/shock-initialization-notes.md` — the acceptance physics (cold-start Hugoniot, C¹ seams, global convexity)
- `doc/eos_implementation_plan.md` — the Cerberus-side engineering plan (Stages 1–5 done)
- `Exec/testing/EOS-Table/README.md` — the frozen `.eostab` v1 spec + per-stage gate records

> **[2026-07-27]** Requirements + acceptance gates for the next
> (SESAME-era) table generation are consolidated in
> `doc/eos_table_conditioning_requirements.md` — headed by the surface-
> monotonicity requirement (H1/G1) the SS3-era gates did not cover.

## Decisions (fixed)

1. **Deuterium first**, titanium second on the proven pipeline.
2. **Cold reference: cryogenic liquid D₂ at ~20 K, ρ₀ = 0.171 g/cc** — matches published Hugoniot
   reference states and the existing `EOS-Hugoniot` `ref_density = 171` (kg/m³).
3. **No SESAME dependency.** The cold/low-T building block is an in-house semi-analytic Python
   module (CoolProp/Richardson-2014 fluid D₂ + Vinet cold curve + Debye/TF thermal).
4. **Tension region (p ≤ 0):** the offline tool clips to a small positive floor (`p_floor = 1e3`
   barye) with hull = 0 there; the C++ in-hull `p_min > 0` constructor check stays (D10 deferred).

Scope: single-material, one-temperature tables for single-fluid hydro.

---

## 1 — Python toolchain (`Exec/python_analysis/eos_tools/`)

`eos_table_prep.py` is a path-stable CLI shim over the package (existing
test `run` invocations unchanged):

```
eos_tools/
  cli.py            # synthetic, fpeos, qa, sources (built) + splice, coldmodel (planned)
  constants.py
  formats/          # eostab.py (frozen-spec writer/reader), fpeos.py (built);
                    # ifpeos.py, expt.py (planned)
  grids.py          # hull-aware PCHIP regrid
  condition.py      # condition, inverse_maps, shift_energy
  thermo.py         # (planned) entropy_from_isochores, maxwell_residual
  splice.py         # (planned) Seam, blend_weight, energy_align, splice()
  hugoniot.py       # RH locus solver
  qa.py             # plots + acceptance checks
  sources.py        # raw-data manifest verify/fetch
  models/           # (planned, SS2) Helmholtz-level analytic sources:
                    #   coolprop_d2, coldcurve (Vinet/BM), qeos (Debye+TF), ideal_plasma
  materials/        # (planned) deuterium.py, titanium.py — MaterialSpec dataclasses
tests/              # pytest; run by the Exec/testing/EOS-PyTools adapter case
```

Do **not** consolidate the test `check.py` embedded readers/Hugoniot solvers into the package —
their independence from the generator is a deliberate gate property.

### Pipeline (`splice --material D`)

1. **Ingest** each source (analytic model or data reader) onto the common grid with its own hull.
2. **Energy-origin alignment** (splice-plan §5.1): constant offset per adjacent pair measured in
   the overlap band; gate: offset constancy std ≤ a few % of local kT (also catches transcription
   errors structurally).
3. **Blend** with `w(lT) = 0.5(1 + tanh(1.47(lT − c)/Δ))` per seam, **T-only**. Key identity: with
   T-only weights, blending specific Helmholtz energy `a` gives `P_blend = w·P_hot + (1−w)·P_cold`
   *exactly*; the only Tier-1↔Tier-2 difference is one E-correction in the band, `−T·w′(lT)·Δa`.
   - **Tier 1 (first):** blend P, E directly; `condition()` the blended surfaces; *measure* the
     Maxwell/Grüneisen residual `R = |dedrho − (p − T·dpdT)/ρ²|/(max(p, T|dpdT|)/ρ²)` per source
     (yardstick) and blended, plus `max|T w′ Δa|/e` in-band. Acceptance (commit measured):
     band median R ≤ 2× per-source baseline; `|T w′ Δa|/e ≤ 1%`.
   - **Tier 2 (`--mode helmholtz`):** entropy per source by isochore integration of conditioned cv
     down from a **single analytic anchor** `s_a(ρ) = s_ideal_plasma(ρ, T_a)`, T_a ~ 3–6×10⁷ K
     (FD electrons; one smooth analytic anchor avoids the per-isochore-constant/Maxwell trap);
     `a = e − Ts`, blend, rederive `E = a − T∂a/∂T` (P unchanged). Gate: band R drops to the
     per-source floor. The cold composite needs no Tier 2 (Helmholtz-level by construction).
4. **Tension clip:** `p < p_floor` → `p_floor`, hull = 0; two-phase dome hull = 0 wholesale;
   monotonise dpdρ. Result satisfies the C++ in-hull `p_min > 0` check.
5. **Condition + emit:** existing machinery; provenance records sources, seams, offsets, clips.
6. **QA:** convexity map, C¹ seam metrics, Maxwell residual, cold-start Hugoniot vs experiment,
   reshock locus, resolution convergence.

### Deuterium sources and seams

*(v2, 2026-07-12 — the full iFPEOS table proved unpublished: the APS SM holds only the
ρ = 0.001 g/cc demonstration isochore. **H-REOS.3** takes the WDM role; iFPEOS demotes to
validation data. Optional upgrade path: request the full table from the LLE authors.)*

1. **Cold composite model** (internal Helmholtz-level join, not a table seam — **built, SS2**):
   CoolProp/Richardson D₂ ⊕ ideal D₂ rotor gas ⊕ QEOS solid (LLNL Vinet cold curve +
   Slater-Debye); exact-derivative blends.
2. **Seam 1** cold model ↔ **H-REOS.3** (mass-scaled H→D): **lT = 2.55 ± 0.15 (band
   ~250–500 K; revised from 3.2 during SS3)** — the C++ self-test caught a 24% p(T) dip from
   the cold model's internal CoolProp→rotor ramp (505–600 K) at liquid densities; with
   REOS.3's 60 K floor the table hands over *below* that ramp and the rotor gas carries no
   splice weight. REOS.3 owns dissociation. Caveats vs the unavailable iFPEOS: classical ions
   (no NQE), PBE-class XC — mitigated by the Hugoniot gates and the iFPEOS isochore
   spot-check (1.9% median). The ρ hand-off to REOS.3 sits at 0.24 g/cc (also
   measurement-driven; see §4).
3. **Seam 2** H-REOS.3 ↔ FPEOS: lT ≈ 5.5 ± 0.25 (~300 kK; REOS.3 tops out at 10⁷ K, FPEOS
   PIMC is the better source at high T). The Hugoniot compression maximum (~30–60 kK) is
   owned by REOS.3 alone.
4. **Seam 3** FPEOS ↔ ideal-plasma model: lT ≈ 7.6 ± 0.2 (FPEOS ceiling 6.4×10⁷ K; the model
   owns the range to 10⁹ K and doubles as the Tier-2 anchor). Gate: P, E agree < 1% in-band —
   pre-validated at SS2 (median 0.06% vs FPEOS at T ≥ 3×10⁷ K).

### Grid / file size

**384×384**, lρ ∈ [−4, 3], lT ∈ [1.25, 9] (20 K = 10^1.30 interior). Bilinear mid-cell error
≈ (n·ln10·h)²/8 → 0.5–1.1% P at the stiff 20 K cold curve (n ≈ 5–7); n_T set by ≥ 10 cells per
seam half-band. ~23 MB ASCII → **commit `D_spliced.eostab.gz`** (~7 MB) + sha256; test `run`
scripts `gunzip -k` (decompressed form gitignored). No git-LFS; no v2 binary format (would reopen
the frozen spec). Committed artifact keeps all test gates CoolProp-free.

Recorded trap: at 20 K the isotherm crosses p = 0 within ~1 cell below ρ₀, so the reference cell's
p is resolution-sensitive at the tens-of-bar level → the cold-shock test gates the pre-shock state
on (ρ, T, e) closure and |p_sim − p_tab| normalized by the **shocked plateau** pressure, not the
existing `PRESHOCK_TOL = 1e-5` relative-p.

---

## 2 — Building-block data harvesting

### Deuterium (all public)

| Piece | Source | Acquisition | Status |
|---|---|---|---|
| WDM/dissociation 60 K–10⁷ K | **H-REOS.3**: Becker et al., ApJS 215, 21 (2014), VizieR J/ApJS/215/21 table2 | scripted | **done 2026-07-12** (106 ρ × 42 T) |
| WDM validation | **iFPEOS**: Mihaylov, Karasiev, Hu, Rygg, Goncharov, Collins, PRB 104, 144104 (2021) | manual | article PDF + partial SM acquired; **full table unpublished** — demoted to validation (ρ=0.001 isochore); optional author request |
| Hot plasma 15.6 kK–64 MK | FPEOS H table (`data/raw/FPEOS/`, tarball archived) | scripted | done |
| Cross-check | arXiv 1306.1902 wide-range DFT H EOS (validation overlay) | scripted | pending |
| Fluid D₂ 18.7–600 K | Richardson-2014 Helmholtz EOS via **CoolProp** (generation-time dep only) | library | pending (SS2) |
| Compressed cold curve | Loubeyre et al. DAC isotherms (Vinet fit input) | manual/digitized | pending (SS2) |
| Hot extension + anchor | in-house `ideal_plasma.py` (FD electrons, Boltzmann ions) | code | pending (SS2) |
| Hugoniot validation | Nellis gas-gun, Knudson Z, Hicks Omega, Boriskov + existing MC2000 | manual/digitized | MC2000 done; rest pending |

### Titanium (SS5, schedule-decoupled)

Vinet fit to Ti DAC + Debye (θ₀ ≈ 420 K) + TF electron via generic `qeos.py`; ML-MD Ti melt data
(arXiv 2603.04680) as constraints; hot end via `ideal_plasma.py` + literature QMD points;
validated vs experimental Ti Hugoniot. No SESAME.

### Manifest & provenance

`Exec/testing/EOS-Table/data/raw/sources.yaml` (`id, citation, doi/url, retrieved, sha256,
license, acquisition: scripted|manual|digitized, notes`). CLI: `sources --verify` (checksums;
always runnable) and `sources --fetch` (scriptable subset). PDF-transcription risk (iFPEOS SM)
mitigated mechanically: pristine PDF + transcribed `.txt` + the re-runnable extraction script all
committed; reader invariants (points-per-isochore, monotonicity, spot values vs paper body); QA
overlay vs digitized figure; and the SS3 alignment-constancy gate catches structural errors for
free. Digitized datasets carry a per-point provenance column.

---

## 3 — Stages SS1–SS5 (gates committed as measured, house style)

**SS1 — Refactor + harvest infra** *(no physics change)*: package + shim; pytest +
`Exec/testing/EOS-PyTools` adapter; `sources.yaml` + fetch/verify; iFPEOS acquisition +
characterization memo; experimental Hugoniot CSVs.
*Gates:* regenerated `D_fpeos.eostab`/`ideal_synthetic.eostab` diff-clean modulo provenance;
existing suite verdict-identical; pytest green; checksums verify.

**SS2 — Cold/low-T semi-analytic module**: `models/*` + composite; `coldmodel --material D --qa`.
*Gates:* ρ(20 K, 1 bar) = 0.171 ± 0.5%; liquid cs within 10% of ~1.1 km/s; 300 K isotherm within
5% of DAC; cv > 0, cs² > 0 on the model rectangle; composite Maxwell residual ≤ 1e-3;
ideal-plasma vs iFPEOS/FPEOS P < 2% at T ≥ 3×10⁷ K.

**SS3 — Splice + `D_spliced.eostab` (Tier 1)**: alignment; clip; 384² emission + gz + QA.
*Gates:* cs² > 0 at every in-hull cell (hard); C¹ seams — max in-band inter-cell jump of
dlnp/dlnT (also dlne/dlnT, dlnp/dlnρ) ≤ 3× the 90th percentile outside; Maxwell residual per
Tier-1 acceptance; cold-start Hugoniot from (0.171, 20 K, 1 bar): gas-gun u_s(u_p) within 2%,
40–200 GPa compressions within experimental 1σ for ≥ 80% of points / 5% all, peak compression
∈ [4.2, 4.9], no d(compression)/d(logP) spike at seam pressures; reshock-locus QA plot;
768² convergence < 0.2% Hugoniot shift.

**SS4 — `Exec/testing/EOS-ColdShock/`** (cloned from `EOS-Hugoniot`): three drive strengths from
(0.171 g/cc, 20 K) + one abusive rarefaction.
*Gates:* pre-shock (ρ,T,e) closure + plateau-normalized p; on-Hugoniot plateau vs the check's own
locus; **hull-clamp counter = 0 for physical runs** and > 0 for the abusive run; conservation;
completion; full suite (incl. untouched `EOS-Hugoniot`) green.

**SS5 — Application wiring + titanium**: (a) the application case → tabulated gas on the spliced
table; smoke gates (completion, clamp budget, compression/timing baseline artifact).
(b) Ti via `materials/titanium.py` — its own SS2′/SS3′ pass, same gate types.
*(Tier-2 Helmholtz mode lands whenever the measured Tier-1 residual motivates it.)*

## Risks

- ~~iFPEOS coverage drives seam placement~~ — retired: coverage characterised (53 ρ × 39 T,
  800 K–256 MK); residual risk is SM raggedness/error bars, checked at transcription.
- CoolProp is generation-time only; the committed `.gz` keeps test gates CoolProp-free.
- QEOS quality at ρ > 2 g/cc, T < 800 K (re-compression region) is only DAC-gated — flagged in
  provenance; revisit if implosion trajectories dwell there.
- Ti hot end is genuinely new work — isolated in SS5, never blocks deuterium.
- `e_shift` cancels in Hugoniot e₂ − e₁ only if every quantity comes from the same table —
  preserved by construction (single spliced table).

---

## 4 — Progress log

### 2026-07-11 — SS1 complete

- `eos_tools/` package landed; `eos_table_prep.py` reduced to the CLI shim (code moved verbatim,
  RH locus solver factored into `hugoniot.py`).
- pytest suite (15 tests) + `Exec/testing/EOS-PyTools` adapter case (pass and fail paths
  verified); `pytest`/`pyyaml`/`coolprop` added to `environment.yml` (pip section).
- `sources.yaml` manifest + `sources --verify/--fetch`; verify green (fpeos, hugoniot-mc2000,
  ifpeos PDF).
- **Measured gates, all green** (full record: `Exec/testing/EOS-Table/README.md` §SS1):
  tables regenerate byte-identical modulo the `generator:` line; EOS-Table PASS; EOS-Sod-Ideal
  PASS (geos round-off 5e-15, wall-geos 1.48×); EOS-Hugoniot PASS (21 gates, clamps-fired
  15598 > 0). No C++/committed-data changes → remaining suite provably unaffected.
- **iFPEOS acquired + characterised** (`data/raw/iFPEOS/README.md`): article PDF (manual
  download — NSF-PAR unreachable, OSTI 500, APS 403 to scripted clients); 53 ρ × 39 T points,
  800 K–256 MK; authors' interpolation corner at ρ ≤ 0.084 g/cc below 182 kK; OFMD↔KSMD energy
  shifts already applied by the authors. **Provisional plan revision: seam 2 dropped (3-source
  splice), seam 3 moves to lT ≈ 8.2** — final call at SM transcription.
- **Open item (manual):** APS Supplemental Material (the actual table) →
  `data/raw/iFPEOS/SM/`, then transcription + invariant checks.

### 2026-07-11 — SS2 complete (cold/low-T semi-analytic module)

- `eos_tools/models/`: `base` (Helmholtz-differencing consistency), `coldcurve` (Vinet +
  `fit_vinet`), `qeos` (Debye ion thermal), `ideal_plasma` (classical ions + Fermi-Dirac
  electrons; doubles as hot extension and Tier-2 anchor), `rotor_gas` (dilute molecular D₂:
  translation + classical rotor σ=2 + harmonic vib, H₂ spectroscopic constants mass-scaled),
  `coolprop_d2` (Richardson-2014 wrapper, two-phase/out-of-range masking), `composite`
  (exact-derivative blending of aligned Helmholtz sources — no numerical differencing of
  blends). `materials/deuterium.py`: `DeuteriumColdModel` = CoolProp ⊕ rotor gas (T-ramp at
  505–600 K) ⊕ QEOS solid (ρ-ramp at 0.355 g/cc, the cold curve's 1.2 GPa point).
- **No memory-sourced constants**: the cold curve is fit to the harvested LLNL-JRNL-686936
  hydrogen 0 K isotherm (`data/raw/llnl_coldcurve/`, re-runnable extraction script +
  invariants); the Debye temperature is *derived* from it (Slater construction,
  θ(ρ₀) ≈ 119 K, consistent with solid-D₂ literature). H-vs-D zero-point caveat recorded in
  the CSV header.
- `coldmodel --material D --qa` — **measured SS2 gates, all PASS** (192×192):
  ρ(20 K, 1 bar) = 0.17177 g/cc (+0.45%, tol 0.5%); cs(20 K) = 1.06 km/s (tol 10% of 1.1);
  Vinet rms 0.5% over LLNL 0–100 GPa; cs² > 0 at all 34 802 in-hull cells; cv > 0 in-hull;
  Maxwell residual median 6.7×10⁻⁵ (gate 10⁻³; p95 1.5×10⁻² sits in the blend bands and is
  grid-truncation of the metric); ideal-plasma vs FPEOS P at T ≥ 3×10⁷ K: **median 0.06%,
  max 3.2%** (gates 2%/5%) — seam 3 and the Tier-2 anchor are pre-validated.
- Recorded nuances: fitted Vinet (v₀ = 4.54 cc/g, B₀ = 0.665 GPa, B₀′ = 5.96) are *effective
  compressed-branch* parameters (the LLNL curve is an exp-6+Vinet merge; the fluid piece owns
  the low-P foot). The solid↔fluid alignment scatter (std_e ≈ 0.45 kT in the 25–50 K
  compressed-liquid window) is the v1 melt smearing — soft-gated, revisit only if the SS3
  Hugoniot foot needs it. CoolProp↔rotor-gas alignment is constant to 0.003 kT (validates the
  rotor-gas model).
- pytest grew to 37 tests (Vinet identities/energy-integral, Debye limits + Mie-Grüneisen,
  FD classical/degenerate limits + grand-canonical consistency, CoolProp closure = two SS2
  gates, exact-blend identities, alignment recovery). `EOS-PyTools` adapter PASS
  (ANSI-stripping fix landed after a colourised-log false negative).
- `thermo.py` (maxwell_residual, sound_speed_sq) landed early — it is the SS3 Tier-1 metric.

### iFPEOS SM status (2026-07-11)

`SM/Supplmental_Material_iFPEOS.pdf` (manual download, in the manifest) is **partial**: 4
pages, Sec. 2 contains only the ρ = 0.001 g/cc isochore (39 T rows — confirming the 39-point
T grid, 800 K → 256 MK, and the exact column format: T (K), P (Mbar) ± σ, E (eV/atom) ± σ;
σ = 0.0 marks the authors' interpolated points). It also corrects the interpolation-region
isochore list to 0.001/0.1/**0.3/0.38** g/cc. **The remaining 52 isochores are in additional
file(s) on the APS SM page — still a manual download.** Transcription reader + seam decisions
stay blocked on that.

### 2026-07-12 — iFPEOS table unpublished; H-REOS.3 promoted (SS3 unblocked)

- User confirmed the 4-page SM is the *only* supplemental material APS provides — the full
  53-isochore iFPEOS table was never published. iFPEOS demoted to validation data (its
  ρ = 0.001 isochore + paper figures); optional upgrade path: e-mail the LLE authors.
- **H-REOS.3 harvested** (scripted, VizieR J/ApJS/215/21 `table2.dat` + ReadMe, checksummed in
  the manifest): 106 ρ × 42 T ≈ rectangular, 3×10⁻⁸–1800 g/cc, 60 K–10⁷ K, columns
  ρ/T/P[GPa]/u[kJ/g], no P ≤ 0 rows, seam-1 band covered at every density. Sources-and-seams
  section revised to v2 (4-source splice: cold ↔ REOS.3 ↔ FPEOS ↔ ideal-plasma; seam 2 back
  at lT ≈ 5.5, seam 3 back at lT ≈ 7.6).

### 2026-07-12 — SS3 complete: `D_spliced.eostab.gz` emitted, all gates green

Full gate record: `Exec/testing/EOS-Table/README.md` §SS3 (artifact sha256 there; raw
`.eostab` gitignored, `.gz` 6.5 MB committed-ready). Highlights: convexity clean over
126,725 in-hull cells; Maxwell median 3.1e-3 ≈ REOS.3's own 3.8e-3 baseline; C1 seams 1–2
orders below threshold; cold-start Hugoniot peak compression 4.398, MC2000 anchor median
1.19%, 384→768 convergence 0.0028% median; iFPEOS isochore overlay 1.92% median; C++
one-zone self-test round-trips at 1e-10 (fd-vs-blocks 49.1 reported-not-gated, tier-2
policy).

Three defects were found by measurement and fixed (details in the README record): seam 1
moved to lT = 2.55 (rotor-ramp dip), ρ hand-off to 0.24 g/cc (Hugoniot discontinuity at
10 GPa from QEOS-solid-above-melt exposure; a T-sliding hand-off over-corrected and was
reverted), and a physics-priority fallback ladder replacing data-source fills (the FPEOS
nearest-fill below its density floor planted non-monotone bumps that only the **C++
round-trip self-test** caught — the offline QA missed it; the one-zone self-test is now a
standing part of the table acceptance loop). A monotone-in-T conditioning pass enforces the
v1 single-branch inversion contract (counts in the table's provenance line).

New tooling: `formats/reos3.py` (invariant-checked reader), `grids.regrid_onto` (FPEOS
regen verified diff-clean), `IdealPlasmaFast` (tabulated FD evaluator, 7e-7 vs exact),
saturated-liquid dome fill (+ sub-18.72 K extension), `splice.py`/`qa_splice.py`,
`splice` CLI subcommand; 45 pytest tests; iFPEOS isochore transcribed
(`extract_sm_isochore.py`).

### 2026-07-12 — SS4 complete: EOS-ColdShock all 24 gates green

Full record: `Exec/testing/EOS-Table/README.md` §SS4. The target physics works end-to-end:
**one cryogenic-liquid start (0.171 g/cc, 24 K) shocked into the gas-gun (6.67 GPa,
compression 2.562), multi-Mbar (176 GPa, 4.156) and classical-plasma (6.26 TPa, 4.064)
regimes, each plateau on the table's own RH locus to ≤ 1.5%, with ZERO hull clamps** —
ionization emergent from the conservative update, exactly the shock-notes requirement.
Initial T is 24 K (not ~20 K): the (0.171, T) cell is fully in-hull only above ~23.5 K
(dome edge + ρ-handoff hull shoulder) — the sanctioned un-ionized relief start.

Notable findings: (i) the table's 17.8 K floor is *hydrodynamically unreachable* — the
FPEOS-style "cool below the floor" abusive recipe produces zero clamps by construction, so
the abusive gate became a **vacuum-forming expansion** (40× cs), which completes with
769,777 tallied clamps and finite fields — closing, for the tabulated gas, the dt-collapse
regime the FPEOS case deferred to W16; (ii) the pre-shock e-closure carries a constant
3.5×10⁻⁷ nondim round-off (committed at 1e-6); (iii) plateau-normalized pre-shock p gate
measured at ≤ 3×10⁻¹⁰ (the §1 trap was real: the cold cell's own p is tens-of-bar
sensitive).

Post-SS4 suite: EOS-PyTools / EOS-Table / EOS-Sod-Ideal / EOS-Hugoniot re-run — all PASS.

### 2026-07-12 — SS5b part 1 complete: titanium solid model (Ti SS2′)

Full record: `Exec/testing/EOS-Table/README.md` §SS5b. `materials/titanium.py` = LLNL Ti
0 K isotherm (extracted, invariant-checked) + Vinet fit + shared `SlaterDebye` +
Sommerfeld electronic term (free-electron, z_c = 4). All QA gates green; the fit's
**B₀ = 109.4 GPa lands inside the DAC literature window with no fitted-to-literature
inputs**, Slater θ₀ = 562 K (vs 420 K literature — known shear-blind bias, reported), cold
sound speed 4.79 km/s (vs ~4.9–5.2). 49 pytest tests. Harvested: NIST Ti I–XXII ionization
energies (scripted) for the Saha hot side.

### 2026-07-12 — SS5b part 2 complete: Saha model + Ti_spliced.eostab (Ti SS3′)

Full record: `Exec/testing/EOS-Table/README.md` §SS5b part 2. `models/saha.py` (NIST Ti
I–XXII, envelope-theorem-consistent (P, e), unit-tested against hydrogenic and limit cases)
+ `splice_ti.py` → `Ti_spliced.eostab.gz` (384², 10⁻³–50 g/cc, 17.8 K–10⁹ K; solid | Saha |
FD ideal plasma, seams lT = 4.5/8.4). All structural gates green (convexity clean over
125k in-hull cells, Hugoniot stride 0.139, C++ round-trips 10⁻¹⁰); seam-3 Saha↔FD
cross-validation at 0.03% max. **Ti v0 accuracy statement**: the seam-1/WDM band is
uncontrolled pending real data — its 4.8 kT alignment scatter is reported, not gated.
54 pytest tests.

### Next
- **Ti v1 data**: ML-MD melt constraints (arXiv 2603.04680 — digitize or author request);
  Ti experimental Hugoniot compilation; then re-seam and gate the WDM band.
- **SS5a**: point the application case at `D_spliced.eostab` (Lua gas block +
  deuteron `ref_*`), smoke gates per §3.
- Manual items open: experimental D₂ Hugoniot compilations (paywalled); optional iFPEOS
  author request.
