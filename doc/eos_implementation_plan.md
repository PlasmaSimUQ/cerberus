# Tabulated EOS for Cerberus — Implementation Plan

Status: planning (2026-07-05; revised 2026-07-07 after a critical review of
the stage plans and a study of Athena++'s general-EOS implementation —
`src/eos/general/` in PrincetonUniversity/athena — which informs D2, D5, D8,
W10 and the new §6 MHD plan; 2026-07-12 folded in W19–W28 and Stages 8–9
from `doc/eos_wide_range_and_mixtures.md`). Companion documents:
- `doc/eos_table_reader.md` — standalone documentation of the `EosTable`
  reader and how the code uses it (kept current as stages land)
- `doc/eos_wide_range_and_mixtures.md` — design note for the wide-range
  single-material stitching (Stage 8, W19–W23) and the composition-weighted
  mixture closure (Stage 9, W24–W28); the physics rationale for both lives
  there, this plan carries the work-item/stage bookkeeping
- `doc/eos_mixture_dalton_plan.md` — detailed Stage-9 bring-up plan for the
  Dalton (partial-pressure) 1T mixture closure (W24–W26, W28-first-cut)
- `Cerberus_TabularEOS_Implementation_Plan.txt` — physics/2T design study (FLASH/Microphysics provenance)
- `MFP_eos_tab2T.H` — annotated draft header for the eventual two-temperature design
- `EOS_Summary.txt`, `EOS_OpenSource_Port.txt` — background and data-supply notes

This document is the *engineering* plan: what gets built in Cerberus, in what
order, against the code as it exists today. Scope for v1 is **CPU-only**,
single-material, one-temperature, hydro states only.

---

## 1. Code audit: where Cerberus assumes the ideal gas law

Key structural fact: Cerberus already has the EOS abstraction. `HydroGas`
(`Source/states/Eulerian/hydro/gas/MFP_hydro_gas.H`) is a virtual interface
(`cons2prim`, `prim2cons`, `define_rho_p_T`, `get_temperature_*`, `get_cp_*`,
`get_gamma_*`, `get_speed_*`) behind a `ClassFactory`, with two backends:
`ThermallyPerfectGas` (tag `thermally_perfect`) and `EilmerGasModel` (tag
`eilmer`). A tabulated EOS is a third backend, not a new framework.

Prim vector layout (`MFP_hydro_defs.H`):
`Density, Xvel, Yvel, Zvel, Prs, Temp, Gamma, SpHeat` + tracers.
- `Gamma` is written by every gas model's `cons2prim` and **read only by the
  three Riemann solvers** — the single flux/thermo coupling point.
- `SpHeat` is written but read by no live code — a **dead slot** available for
  repurposing (it will carry sound speed in `eos` flux mode, §2 D2).

### Hard gamma-law algebra downstream of the gas interface (switch sites)

| Site | Hard-coded expressions |
|---|---|
| `riemann/MFP_hllc.cpp:44,50,59,60,71-72,78,87` | face energy `p/(γ−1)+ke`; `a=√(γp/ρ)`; PVRS p* estimate; two-shock coefficient `(γ+1)/2γ` |
| `riemann/MFP_hydro_hlle.cpp:43-44,57-58,65-66` | energy + sound-speed reforms; Einfeldt `u±a` speeds |
| `riemann/MFP_ausmdv.cpp:61-69,77` | energy + sound speed; enthalpy `(E+p)/ρ` |
| `actions/MFP_CTU_Braginskii.cpp:2573-2574,2606,2743,2755` | recomputes `p=(γ−1)(E−ke)` and ideal-gas `T=pm/ρ` from a cached per-cell γ instead of calling the gas interface |

### Soft (structural) ideal-gas assumptions

- Positivity machinery tests **p > 0** everywhere: `get_positive_prim`
  (`MFP_hydro.H:56`), `apply_prim_floor`, `prim_valid`, the reconstruction
  fallback (`MFP_eulerian.cpp:648-666`), the pre-Riemann clamp
  (`MFP_eulerian.cpp:1100`). For an ideal gas p>0 ⇔ e_int>0; for a table with
  a cold curve they are **not** equivalent (cold compressed solid can carry
  valid energy at near-zero or negative pressure). The floor *variable*, not
  the mechanism, is the ideal-gas assumption.
- `cons_valid` checks total `Eden > 0`, never internal energy — weak for any
  EOS, wrong shape for a table hull.
- Reconstruction (`MFP_eulerian.cpp:614`) treats `Gamma/Temp/SpHeat` as
  independent linear face fields; harmless plumbing, but it is what feeds the
  solvers' γ. (Athena++ hard-aborts *characteristic-projection* reconstruction
  with a general EOS — `reconstruction.cpp:104` — because the eigenvector
  algebra assumes γ-law. Cerberus reconstructs primitives per-component, so
  that failure class does not arise here; noting it in case characteristic
  limiting is ever added.)

### Confirmed clean (through the interface, or EOS-agnostic)

CTU flux path (`MFP_CTU_hydro.cpp`), CFL timestep
(`get_allowed_time_step` → `get_speed_from_cons`), viscous dt + transport
coefficients (interface `get_gamma/get_cp/get_temperature`), shock detector
(pure pressure ratio), hydro BCs (no characteristic/Riemann-invariant γ),
plasma5 / Lorentz / acceleration / reactions / gas-kinetics, refinement,
plotfile output, ICs (`define_rho_p_T` + `prim2cons`).

### MHD: out of scope, guarded

The MHD state has **no gas interface**: it carries a constant `Real gamma`
member (`MFP_mhd.H:33`) and hard-codes gamma-law algebra throughout
(`MFP_mhd.cpp:327,360,403,411,434`; `mhd/riemann/MFP_mhd_hlle.cpp:57,77,91,
106,131,142`). The tabulated EOS will **not** support MHD initially; the
eventual lifting path (Athena++-style, table sound speed inside the fast
magnetosonic speed) is planned in §6.
Two guard actions are part of this plan (work item W14):
1. Comments in `MFP_mhd.H` (at the `gamma` member and class header)
   documenting that MHD assumes a constant-γ ideal gas and is incompatible
   with the tabulated EOS.
2. A config-time check: if any MHD state coexists with a hydro state whose
   gas model is `tabulated`, `amrex::Abort` with a clear message. Failing at
   `read_config` time beats silently running inconsistent physics.

---

## 2. Design choices

**D1 — Placement.** `TabulatedEOS : public HydroGas`, factory tag
`"tabulated"`, one `Register()` line, configured per-state:
`gas = { type='tabulated', table='D.eostab', small_temp=..., small_dens=...,
ttol=... }`.

**D2 — Flux thermodynamics: a three-mode switch.** Three named modes,
realised as bookkeeping over the (gas type, Lua `flux` key) pair — not as
branches inside solvers (revised; see mechanism below):

- **`gamma`** — the existing code path, byte-for-byte untouched (the
  existing solver files are never edited). Active whenever the tabulated EOS
  is not in use; γ-law and Eilmer states behave exactly as today, so the
  existing test suite is a regression guarantee.
- **`effective_gamma`** — the tabulated gas fills the `Gamma` slot with the
  energy-consistent γₑ = 1 + p/(ρ e); solvers unchanged. **Default for
  `type='tabulated'` during first implementation.** (Note: the acoustic
  Γ₁ = ρc²/p differs from γₑ off the ideal corner; face energy is kept exact,
  the solver's internal wave-speed estimate is approximate — acceptable for
  HLL-family estimates. CFL always uses the true table sound speed via
  `get_speed_from_cons`.)
- **`eos`** — face sound speed and internal energy come from direct
  gas-model evaluations inside a **dedicated Riemann solver** (below) instead
  of γ-algebra. **Becomes the default for tabulated states after Stage 5 is
  validated**; `gamma` / `effective_gamma` remain selectable overrides for
  debugging and A/B comparison.

Mechanism for `eos` mode — **a separate solver file, not branches in the
existing solvers** (revised after the Athena++ review; supersedes the earlier
SpHeat/e-slot idea): a new factory-registered solver
`riemann/MFP_hllc_general_eos.{H,cpp}` (tag `HLLC_general_eos`), selected
through the existing per-state Lua `flux` key (`MFP_hydro.cpp:240`,
`HydroRiemannSolverBuilder` dispatch at `MFP_hydro_riemann.H:27`). It follows
the structure of Athena++'s `hllc.cpp` `GENERAL_EOS` branches:

- face internal energy `e(ρ, p)` and sound speed `a²(ρ, p)` from direct
  gas-model evaluations at the reconstructed face states (one inverse-map
  table lookup each, see D5 — this is why they are affordable per face);
- the two-shock q-factor uses a *local* effective gamma `γ = a²ρ/p`
  evaluated at the PVRS mid-state pressure (Athena++ `hllc.cpp:89-103`),
  not the face `Gamma` slot;
- requires two new `HydroGas` virtuals with ideal-gas base defaults,
  `get_internal_energy_from_rho_p` and `get_sound_speed_from_rho_p`, so
  *every* gas model (including TPG) can run the new solver — which gives the
  key validation: `HLLC_general_eos` on a γ-law gas must agree with `HLLC`
  to round-off, isolating solver bugs from EOS bugs;
- the solver holds a non-owning pointer to its state's `HydroGas`, wired in
  `HydroState::set_flux()` (`MFP_hydro.cpp:226`).

Existing solvers (`MFP_hllc.cpp`, `MFP_hydro_hlle.cpp`, `MFP_ausmdv.cpp`) are
**never edited** — the `gamma` path stays byte-identical by construction.
General-EOS variants of HLLE/AUSMDV are optional later work on the same
pattern.

**D3 — Scope (v1).** Single material, one temperature, `total_equi`-style
table (FPEOS deuterium first) on hydro states only. Deferred: 2T ion/electron
split, `rho_lookup` cross-state coupling, mixture/MTMMMT driver, MHD.

**D4 — Table format and curation.** One canonical pre-conditioned `.eostab`
format, produced by an **offline Python tool**: ingest FPEOS/raw dumps
(SESAME later), perform monotonisation (PCHIP), cv > 0 flooring, Maxwell
construction where needed, and emit QA plots. All curation is offline,
inspectable, and versioned; the C++ reader parses only the canonical format.
No EOSPAC dependency in v1.

**D5 — Interpolation and inversion.** Bilinear in (log₁₀ρ, log₁₀T) with hard
hull clamping. The Newton slope for inversions comes from the **analytic
derivative of the interpolant itself** (prevents stalls when smoothed
derivative tables disagree with the value surface); the smoothed derivative
tables are used only for outputs (cs, dpde, dpdr_e — the Riemann contract).
Inversions (re, rp, tp) are hull-bracketed Newton with bisection fallback;
non-convergence → clamp to floor + flag. e(T) monotonicity is NOT assumed.
Bicubic Hermite is a later upgrade on the same layout.

**Inverse maps (Athena++-informed).** Athena++ eliminates runtime inversion
entirely by tabulating pre-inverted surfaces offline (`p/e(e/ρ,ρ)`,
`e/p(p/ρ,ρ)`, `a²ρ/p(p/ρ,ρ)`). We adopt this in hybrid form: the offline
tool also emits optional blocks `T(ρ,e)` and `T(ρ,p)` on their own
log-uniform e/p axes. Runtime inversions then start from one bilinear
inverse-map lookup and take ~1–2 Newton polish steps against the forward
surface to reach `ttol` — the guarded Newton+bisection driver survives as
the fallback and as the consistency guarantee (Athena++ skips the polish and
accepts a finite p→e→p round-trip error; we do not, because the Godunov
update conserves e). The (ρ,T) surface stays canonical — T is needed for
transport, collisions, and the Stage-6 2T split. Off-table queries **clamp
and flag** (Athena++ silently linear-extrapolates off its grid; unsafe below
a cold curve). A guarded linear extrapolation off the hull edge — continuing
`p`, `e` along the stored edge derivatives for a bounded distance — is
recorded as the alternative if the clamp plateau itself causes artefacts;
see STAGE4.md W8.1 for the trade-off. Clamp-and-flag stays the v1 default. Athena++'s dimensionless-ratio storage (log of `p/e` etc.,
exactly constant for ideal gas) is noted as a possible conditioning
refinement for the derived blocks, not adopted in v1.

**D6 — CPU-only (v1).** Table storage is plain host memory
(`std::vector<Real>` in an owner class; raw-pointer POD "view" struct for
evaluation). No GPU arenas, no managed memory, no device decorations.
Forward-awareness is preserved by *structure only*: the init/evaluate split
and the POD-view calling convention are what a future GPU pass needs; the
allocation swap is then a one-file change.

**D7 — Units.** Nondimensionalise the entire table once at load using the
state's Lua reference quantities (ρ→ρ/ρ₀, T→T/T₀, P→P/(ρ₀u₀²), e→e/u₀²,
cv→cv·T₀/u₀²). Evaluation code never converts units.

**D8 — Newton initial guess.** Primary: the D5 inverse-map lookup (already
within interpolation error of the answer, so polish is 1–2 iterations).
Fallbacks where the maps are absent: analytic ideal-gas estimate from the
incoming state, then hull midpoint; the prim `Temp` slot as a warm start only
where a call site provably has valid previous primitives. No persistent
temperature state variable in v1 (avoids checkpoint/plotfile/regrid churn for
a second-order benefit).

**D9 — Parallel loading.** Every rank reads the table file (consistent with
Cerberus's every-rank-runs-Lua config model). IOProcessor + Bcast, or the
node-shared MPI-3 window pattern from the STL/SDF work, are later drop-ins.

**D10 — Floors.** The floor/validity *variable* becomes gas-model-defined:
γ-law states keep p>0 exactly as now (base-class default, zero behavior
change); the tabulated model validates against its hull (ρ, T bounds) and
internal energy, since p>0 is the wrong test below the cold curve. Maps
`small_temp`/`small_dens` onto the existing `effective_zero` conventions.
*Stage-3 refinement (see STAGE4.md):* `prim_valid`/`cons_valid`
(`MFP_hydro_gas.cpp:152,166`) are **non-virtual and `Abort`** — `prim_valid`
tests `p>0`, `cons_valid` only `Eden>0` (never internal energy), and
`apply_prim_floor` (`:88`) floors to an absolute `1e-14` (the open dt-collapse
pathology). For the v1 FPEOS table the whole hull is `p>0`, so v1 does **not**
virtualize `prim_valid`/`cons_valid`; instead it (i) makes `apply_prim_floor`
virtual and floors the tabulated gas to the *physical* hull edge, and
(ii) makes the `cons2prim` hull-`e` clamp the real guard (it catches the
below-hull `e_int` that `cons_valid`'s `Eden>0` cannot). Virtualizing the
validity tests is the cold-curve lifting task, deferred until a cold-curve
table is actually added.

---

## 3. Work items

| # | Item | Effort | Depends |
|---|---|---|---|
| W1 | Baseline `Exec/testing/run_tests.py` pass-state record (regression anchor) + **run-to-run determinism check** (decides which later gates may demand bit-identity vs verdict-identity) | S | — |
| W2 | Canonical `.eostab` format spec (grids, values, smoothed derivatives, **optional inverse-map blocks T(ρ,e)/T(ρ,p)**, unit metadata, provenance header) | S | — |
| W3 | Offline Python conditioning tool + FPEOS-D ingest + QA plots (`Exec/python_analysis/`) | M | W2 |
| W4 | `EosTable` owner/reader + POD view: bilinear interp, guarded re/rp/tp inversions, host-only storage, nondimensionalisation | M | W2 |
| W5 | One-zone unit tests: rt→re→rt, rt→rp→rt round trips over the hull; cs/Γ₁ vs centred finite differences (debug-gated Lua `self_test`, mirroring the SDF pattern) | S–M | W4 |
| W6 | `TabulatedEOS : HydroGas`: all virtual overrides, effective-γ fill of `Gamma` slot, Lua schema, factory registration | M | W4 |
| W7 | Flux-mode resolution plumbing only: config-time defaults by gas type + Lua override + validation (`flux='HLLC_general_eos'` aborts until W10). **No edits to existing solvers** — the `eos` path is W10's separate solver file | S | — (parallel) |
| W8 | Floors generalization (D10): (a) `HydroGas::apply_prim_floor` made `virtual`, tabulated override floors to hull `ρ_min`/`p(ρ,T_min)` not `1e-14` (resolves the open absolute-floor dt-collapse pathology, scoped to the table); (b) `cons2prim` clamps `e_int` to the hull `e`-column before inversion — the load-bearing guard `cons_valid`'s `Eden>0` test cannot provide; (c) the `clamped` flag made load-bearing (per-box counter + verbosity report); (d) `prim_valid`/`cons_valid`/`get_positive_prim` **verified never-fire** on FPEOS (whole hull `p>0`) rather than virtualized — cold-curve virtualization documented-but-deferred | M | W6 |
| W9 | Test cases in `Exec/testing/`: Sod into the table's ideal-gas corner vs γ-law analytic; single-material strong-shock Hugoniot vs published FPEOS locus | M | W6, W3 |
| W10 | **DONE (Stage 5)** `eos` flux mode = new solver `riemann/MFP_hllc_general_eos.{H,cpp}` (Athena++ `hllc.cpp` GENERAL_EOS pattern, see D2): face (e, a) from a combined gas-model (ρ,p) evaluation, q-factor from local γ=a²ρ/p at the **face** (not PVRS p*, changed after the G8 wall measurement); `HydroGas` virtuals `get_internal_energy_from_prim`/`get_sound_speed_from_prim_rp`/`get_face_eval_from_prim` with Gamma-slot ideal defaults; round-off vs `HLLC` on TPG = 5e-15; A/B vs `effective_gamma` on W9 (Hugoniot gates unchanged); default `flux` flipped for tabulated states | L | W7, W9 |
| W11 | **DONE (Stage 5)** Braginskii consistency: `T_e`/`T_i` reform in `rhs` (`:2573-75`, `:2605-07`) + positivity bailouts in `check_invalid` (`:2791-94`, `:2802-05`). Corrected from plan: these are `static` functions → dispatch via static `s_*_gas_tab` pointers set in `calc_time_derivative`, not `ion_state`; only **T** needed reforming (`p_e`/`p_i` are intermediates), via `species_temperature_tabulated` → the table `T(ρ,e)` map. γ-law path bit-identical (G7 fcompare) | M | W7 |
| W12 | 2T design memo (paper before code): `rho_lookup` data flow — which state's data, at which RK stage — and the electron wave-speed chain rule (∂P_e/∂ρ_e vs table axis ∂P_e/∂ρ_material) | M | parallel |
| W13 | 2T components (SESAME 303/304 split, per-state component binding, 303+304=301 init assertion, exchange-term cv). The mixture/MTMMMT driver is no longer here: W13 **specialises the Stage-9 1T machinery** (W24–W27) to electron/ion at separate temperatures, reusing its per-component binding, mixing-rule interface, and mixture-sound-speed chain rule | L | W10, W12, W24–W27 |
| W14 | MHD guard: incompatibility comments in `MFP_mhd.H` (class header + `gamma` member); config-time `amrex::Abort` if an MHD state coexists with a `tabulated` hydro gas | S | W6 |
| W15 | MHD with general EOS (future, demand-driven; Athena++ `general_mhd.cpp` pattern — see §6): MHD gas-closure abstraction, table-backed fast magnetosonic speed, `MFP_mhd_hlle_general_eos` solver, guard replaced | L | W10 |
| W16 | Separate `rho_floor`/`p_floor` (later goal; see STAGE4.md deferred section): per-variable floors replacing the single `effective_zero` for all gas models + MHD, closing the absolute-floor dt-collapse pathology globally (Stage-4 W8.1 closes it for the tabulated gas only); includes the vacuum-cell velocity decision | M | W8 |
| W17 | Checkpoint/restart test with the tabulated gas (later goal; see STAGE4.md deferred section): stop/restart/`fcompare`-vs-uninterrupted bit-identity gate; no new StateData expected, table reloads per rank at config time | S | W9 |
| W18 | **Later** — store `.eostab` payload in **SI** instead of CGS (units = the dimensional system a quantity is expressed in; SI = kg·m·s, CGS = g·cm·s, so e.g. density switches g/cm³→kg/m³, pressure erg/cm³→Pa, a factor of 10). Currently the table columns are written and read in CGS while Cerberus reference quantities are SI, so the reader ctor does an SI→CGS round-trip; storing SI removes that mismatch and one conversion class. Scope: (a) `.eostab` format spec `units` metadata field (W2) declares `SI`, with a version bump + a back-compat reader branch so existing CGS tables still load by their header tag; (b) the Python conditioning tool (W3, `eos_table_prep.py`) emits SI columns and the SI provenance header; (c) the `EosTable`/`TabulatedEOS` reader (W4/W6) drops the SI→CGS step and nondimensionalises directly against the SI reference quantities; (d) regenerate `EOS-Table/data/*.eostab`; (e) gate: one-zone self-test (W5) round-trips + `EOS-Sod-Ideal`/`EOS-Hugoniot` verdicts unchanged after regeneration (numbers are dimensionless post-nondimensionalisation, so gates should be bit-identical modulo table round-off) | M | W2, W3, W4, W6 |

### Stage 8 — wide-range single-material stitching (W19–W23)

Rationale and physics in `doc/eos_wide_range_and_mixtures.md` §2–4. Offline-heavy;
no C++ solver/reader changes with W23 deferred. Mnemonics WA/WB/WC/WE/WD.

| # | Item | Effort | Depends |
|---|---|---|---|
| W19 (WA) | `eos_table_prep.py` **`stitch` subcommand**: ingest N source sub-EOS files (each with its own ρ,T support + region tag), build the union log-uniform (ρ,T) grid, regrid each source onto it. `source` metadata becomes a **list**; add a per-cell **provenance-map block** (dominant source per cell) for QA. `.eostab` format version bump + back-compat reader branch so single-source CGS/SI tables still load by header tag | M | W2, W3 |
| W20 (WB) | **Overlap blending of p and e**: per-seam weight function (smoothstep/tanh in log-T, or along a supplied boundary curve) across the solid/WDM/plasma seams | M | W19 |
| W21 (WC) | **Continuity + monotonicity conditioning across crossovers**: preserve ionisation/melt softening and the α→ω / ω→β energy jumps as steep **monotone** crossovers; keep p(ρ)_T monotone; cv > 0 in single-phase regions; latent-heat energy retained in e. No Maxwell construction (UC1-d) | M | W20 |
| W22 (WE) | **Validation gates**: (a) principal shock Hugoniot through the transition region vs published Ti data; (b) seam continuity of p, e, cₛ within tolerance; (c) a solid→plasma shock tube runs stably; (d) reader inversion + CFL survive the softest (ionisation/melt) point | M | W19, W21, W3 |
| **W23 (WD)** | **DEFERRED — exact coexistence plateau + energy-map-primary inversion.** Only needed for split-shock (two-wave) structure near a transition threshold; out of scope for the overdriven, strength-free target regime (`doc/eos_wide_range_and_mixtures.md` §3). Adds a Maxwell/tie-line construction + an energy-primary inversion branch in degenerate mixed-phase cells + a coexistence flag | (deferred) | W21 |

### Stage 9 — composition-weighted mixtures, 1T (W24–W28)

Rationale in `doc/eos_wide_range_and_mixtures.md` §5; detailed Dalton bring-up
plan in `doc/eos_mixture_dalton_plan.md`. Hot-loop C++; builds on Stage 5
(`HLLC_general_eos`, done) and the already-advected α-fractions (done). This is
the **1T parent of W13**. Mnemonics WF–WJ.

**Status 2026-07-12: W24, W25, W26 and the Dalton subset of W28 are DONE**
(`MFP_mixture_gas.{H,cpp}`, tag `tabulated_mixture`; validation in
`Exec/testing/EOS-Mixture/` — pure-cell round-off gates pass **bitwise**,
per-component conservation exact, wall ratio 2.39 ≤ 2.5 for N=2). One
implementation addition beyond the written plan: below-hull partial densities
use the ideal-gas low-density limit (p, dpdT scaled by ρₖ/rho_hull_min at the
hull-edge row) instead of a plain clamp — a plain clamp gives a dilute
component the FULL p(rho_hull_min,T) partial pressure and a visible
drop_tol kink; with the scaling the measured kink is 3e-12. W27 (Amagat)
remains open.

| # | Item | Effort | Depends |
|---|---|---|---|
| W24 (WF) **[DONE]** | **Mixing-rule interface + N-table binding**: a gas (`type='tabulated_mixture'`, or `tabulated` with a `components={{table=,name=},…}` list) owns N `EosTable`s bound to component names matching the αₖ tracer slots; define the mixing-rule virtual interface (forward p,e / `cons2prim` inverse / sound speed). Generalise the `set_flux` default-flux selection from a literal `TabulatedEOS::tag` compare to a `needs_general_eos_solver()` virtual | M | W6, W10 |
| W25 (WG) **[DONE]** | **Dalton (partial-pressure) closure** — the cheap first cut: ρₖ=αₖρ, p=Σpₖ(ρₖ,T), e=Σαₖeₖ(ρₖ,T); `cons2prim` = 1-D Newton in T from total e_int (component densities known from αₖρ). Pure-cell short-circuit; per-component hull clamp+flag (+ dilute ideal-limit scaling, see status note) | M | W24 |
| W26 (WH) **[DONE]** | **Mixture frozen sound speed + face eval**: a² = (∂p/∂ρ)\|_{s,α} via chain rule through the mixing rule (frozen composition; the wave-speed chain-rule subtlety flagged in the W12 2T memo); implement `get_face_eval_from_prim`/`get_speed_from_*` and the Gamma/SpHeat slot fill for the mixture so `HLLC_general_eos` drives it | M | W24, W25 |
| W27 (WI) | **Amagat (additive-volume, P–T-equilibrium) closure** — the accurate mode, selectable against W24's interface: 1/ρ=Σαₖ/ρₖ(p,T), e=Σαₖeₖ(p,T); ~~`cons2prim` = 2-D Newton in (p,T)~~ nested guarded 1-D solves (outer T, inner p from the strictly-monotone volume constraint — see the plan doc §2 for the revision rationale); sound-speed chain rule for the volume-constraint form; per-component non-overlapping-hull fallback. **Detailed implementation + validation plan: `doc/eos_amagat_plan.md` (2026-07-26; AM1–AM3 done — drivers complete; AM10 solver-optimisation item added for production-scale cost: the measured mixed-cell factor is 8.4× Dalton and the drop_tol-band diffusion tails set the effective mixed fraction)** | L | W24, W26 |
| W28 (WJ) **[Dalton subset DONE]** | **Validation**: (a) a mixture of two **identical** component tables reproduces the single-table result to round-off (mixing-rule analog of the Stage-5 G1 gate) — done, on IDENTICAL ideal-gas tables where Dalton is exact; (b) a binary-mixture shock tube — done (`EOS-Mixture/`); (c) Dalton-vs-Amagat A/B on a dilute case where they should nearly agree — needs W27; (d) mass-fraction conservation + positivity across a steepening composition gradient — done | M | W25, W26, W27 |

### Note (2026-07-26) — W29 **[DONE 2026-07-26]**: contact-resolution switches + attribution for the general-EOS HLLC fallback

**[Landed 2026-07-26: switches `fallback_guard_{eval,sound,wave}` + counters/report always compiled; neutrality bitwise; controls verified — record in `doc/eos_amagat_plan.md` §8 result note.]**

The in-solver HLLE fallback (`doc/general_eos_hllc_fallback_plan.md`) makes
`HLLC_general_eos` crash-proof, but every substituted face trades the
contact-resolving three-wave flux for the two-wave HLL average — and
contact smearing is exactly what the EOS work must not suffer *silently*:
the mat_a/mat_b application case's ℓ≥2 drain is an HLLE-at-contact pathology
(`doc/t4_fix_plan_floorB_cavitation.md` Part I-B), and the Amagat
validation gates measure interface width directly. Wanted, when scheduled:

1. **Per-guard Lua switches** on the fallback (enable/disable guards
   1/2/4 individually, plus a global off; default = current behaviour), so
   a run can (a) demonstrate where HLLC alone is genuinely singular and
   (b) A/B contact quality with and against specific guards.
2. **When-and-where accounting in production builds** — the per-guard,
   per-block counts with the block-identity line (designed and working
   under `MFP_SOLVER_DIAG`, fallback plan §5) promoted to a
   verbosity-gated standard report, so a contact-resolution violation is
   attributable in *every* run, not only in diagnostic builds.
3. The same accounting contract binds **any future flux limiting (option
   A)**: a limiter that scales a face flux must report when and where it
   acted exactly as the fallback does, or contact-quality regressions
   become undiagnosable.

Cross-refs: fallback counters `general_eos_hllc_fallback_plan.md` §5;
guard-2 diversion counter D-R7 in `t4_fix_plan_floorB_cavitation.md`
§II.2.0; `doc/eos_amagat_plan.md` §5 note (the B/C-tier measurements this
protects).

---

## 4. Staged implementation

**Stage 1 — Data groundwork** *(small; W1–W3; no Cerberus source changes)*
Offline tool, canonical format, conditioned FPEOS deuterium table.
Exit: a `.eostab` whose isotherms, Hugoniot, and cv-positivity have been
plotted, eyeballed, and committed alongside the generating script.
Detailed Stage-1 plan and artifact home: `Exec/testing/EOS-Table/README.md`
(baseline record in `baseline/`, raw + conditioned tables in `data/`,
QA plots in `qa/`; no `check.py` until Stage 2's one-zone test).

**Stage 2 — Core machinery + config plumbing** *(medium; W4, W5, W7)*
Detailed Stage-2 plan: `Exec/testing/EOS-Table/STAGE2.md`.
Table code exercised by one-zone tests only; W7 config plumbing lands with
**no existing-solver edits**. Exit: round trips converge to `ttol` over the
whole hull measured as the **residual in e/p** (the conserved quantities; the
T error is monitored but unbounded where cv→0 makes T ill-posed given e),
derivative identities match finite differences, **and** the full existing
test suite matches the W1 baseline verdict-for-verdict, plus an `fcompare`
field-identity spot-check (W1 finding: all physical fields are
bit-reproducible run-to-run; only the `cost` load-balance diagnostic
differs) — proving the additions inert before new physics is reachable.

**Stage 3 — First plasma: `effective_gamma` mode** *(medium; W6, W9-Sod, W14)*
Detailed Stage-3 plan: `Exec/testing/EOS-Table/STAGE3.md`.
Tabulated gas wired end-to-end via the existing interface. MHD guard lands
here (W14) since this is the first point a user could configure the
incompatible combination. Exit: Sod-into-ideal-corner matches the γ-law
analytic; suite still green.

**Stage 4 — Robustness + validation gate** *(medium; W8, W9-Hugoniot)*
**DONE 2026-07-08** — detailed plan + results: `Exec/testing/EOS-Table/STAGE4.md`.
Key deviation from plan: the *density* axis, not the energy axis, was the
load-bearing clamp (`locate()` clamped off-hull ρ silently; the fix is
`clamp_rho()` at every gas-model entry before anything is derived).
Hugoniot gate: measured compressions on the table-RH locus to 0.7–1.3%
across 247 GPa–25 TPa; abusive rarefaction completes with 10,966 tallied
clamps and zero validity aborts.
Swaps the synthetic table for the **real conditioned FPEOS deuterium table**.
Hull-aware floors, survival of a deliberately abusive strong-shock run
(bisection fallback exercised, flagged, non-fatal), Hugoniot overlay against
the published FPEOS locus. Reframed after the Stage-3 verification reads
(see STAGE4.md): the FPEOS hull is entirely `p>0` (no cold curve at T≥1.35 eV),
so the load-bearing robustness is **hull-edge / negative-internal-energy
clamping**, not literal `p<0`. The three floor/validity mechanisms
(`apply_prim_floor` → hull-physical floor; the `cons2prim` hull-`e` clamp made
load-bearing with a counter; `prim_valid`/`cons_valid`/`get_positive_prim`
verified never-fire rather than virtualized) are detailed there; the
cold-curve virtualization is documented-but-deferred.
**Milestone gate: everything after this stage may be revised based on what
Stages 1–4 teach.**

**Stage 5 — `eos` flux mode becomes default** *(large; W10, W11)* —
**DONE 2026-07-09** (STAGE5.md has the full results notes). Detailed plan +
design decisions D-a…D-g and gates G1–G8 in `Exec/testing/EOS-Table/STAGE5.md`.
The `HLLC_general_eos` solver (Athena++ pattern, D2) is now the default for
tabulated states, validated in three steps: (1) round-off agreement with
`HLLC` on a γ-law gas — L∞/range 5e-15, isolating solver bugs from EOS bugs;
(2) A/B against `effective_gamma` on the Stage-4 cases — the 21 Hugoniot
gates re-passed unchanged (compressions 0.67/1.07/1.27%, jump is
conservation-set so both modes agree); (3) default `flux` flipped for
tabulated states (γ-law states run the untouched existing solvers, byte-
identical by fcompare). Braginskii temperature (W11) reformed behind static
dispatch pointers — γ-law transport path bit-identical (G7). Key deviations
from the written plan, both in W11: the reform sites are `static` functions
(dispatch via static pointers, not `ion_state`), and only *temperature*
needed reforming (`p_e`/`p_i` are intermediates). q-factor uses the face Γ₁
(not PVRS-p*) after the G8 wall-time measurement; wall budget 2.0× (measured
1.57×). The combined `get_face_eval_from_prim` virtual (one inversion for
both face e and a) was added for that budget.

**Stage 6 — Two-temperature** *(large; W12 memo written early, W13 gated on it)*
Per `Cerberus_TabularEOS_Implementation_Plan.txt`, sequenced 2T-first (true
multifluid is the target use case); the mixture driver stays parked until a
concrete need. Note: write the W12 memo during Stages 1–2 regardless — whether
an evaluation ever takes two densities affects the `EosTable` view API.

**Stage 7 — MHD with general EOS** *(large; W15; optional, demand-driven)*
See §6. Not scheduled until a concrete MHD+tabulated use case exists; the
W14 guard stays in force until this stage replaces it.

**Stage 8 — wide-range single-material stitching** *(medium; W19–W22; W23 deferred)*
Detailed design + physics: `doc/eos_wide_range_and_mixtures.md` §4. Offline-heavy
— stitch several source sub-EOS (solid → WDM → plasma) for one material into a
single wide-range `.eostab`; the reader, gas model, and solver are unchanged
(one log-uniform (ρ,T) grid out). Sub-stage split: **8a** = W19 + W20 +
single-phase gates (ships a working wide-range table); **8b** = W21 + the
transition/cusp gates (gets the α→ω→β / ionisation crossovers right).
Exit: principal Ti shock Hugoniot through the transition region matches
published data; a solid→plasma shock tube runs stably; suite still green.

**Stage 9 — composition-weighted mixtures, 1T** *(large; W24–W28)*
Detailed design: `doc/eos_wide_range_and_mixtures.md` §5; Dalton bring-up plan:
`doc/eos_mixture_dalton_plan.md`. Hot-loop C++ — a mixture gas binds one
`EosTable` per component to the existing α-fraction tracer slots and closes the
state with a selectable mixing rule (Dalton first, Amagat accurate) driven
through the Stage-5 `HLLC_general_eos` solver. No new transport (composition is
already advected). The **1T parent of W13**. Recommended after Stage 8 so
mixtures can bind realistic wide-range component tables. Exit: identical-tables
round-off gate (W28a) passes; a binary shock tube runs; mass fractions conserve
and stay positive across a steepening composition gradient.

---

## 5. Reference: source-code touch map

| Area | Files | Stage |
|---|---|---|
| New gas model | `Source/states/Eulerian/hydro/gas/MFP_tabulated_gas.{H,cpp}`, `MFP_eos_table.{H,cpp}` | 2–3 |
| Offline tool | `Exec/python_analysis/eos_table_prep.py` (+ format doc) | 1 |
| Flux mode | config plumbing only (`MFP_hydro.cpp` defaults/validation); **new file** `riemann/MFP_hllc_general_eos.{H,cpp}`; existing solvers untouched | 2, 5 |
| Gas-model face evals | `MFP_hydro_gas.{H,cpp}` — new virtuals `get_internal_energy_from_rho_p` / `get_sound_speed_from_rho_p` with ideal-gas defaults | 5 |
| Floors | `MFP_hydro_gas.{H,cpp}`, `MFP_hydro.H`, `MFP_eulerian.cpp` (fallback interplay) | 4 |
| Braginskii | `Source/actions/MFP_CTU_Braginskii.cpp` | 5 |
| MHD guard | `Source/states/Eulerian/mhd/MFP_mhd.H` (comments), config-time check (`MFP_config.cpp` or `MFP_mhd.cpp` init) | 3 |
| MHD general EOS (future) | `mhd/MFP_mhd.{H,cpp}` closure abstraction, **new file** `mhd/riemann/MFP_mhd_hlle_general_eos.{H,cpp}` | 7 |
| Wide-range stitching (offline) | `Exec/python_analysis/eos_tools/` — new `stitch` subcommand + blending/conditioning; `.eostab` format version bump + back-compat reader branch in `MFP_eos_table.cpp` | 8 |
| Mixture gas model | **new file** `Source/states/Eulerian/hydro/gas/MFP_mixture_gas.{H,cpp}` (owns N `EosTable`s + mixing-rule); `MFP_hydro_gas.{H,cpp}` (`needs_general_eos_solver()` virtual); `MFP_hydro.cpp` (`set_flux` default-flux generalisation) | 9 |
| Tests | `Exec/testing/EOS-Table/` (Stage-1 data + Stage-2 one-zone; see its `README.md`/`STAGE2.md`), `EOS-Sod-Ideal/`, `EOS-Hugoniot/`, `EOS-Mixture/` (Stage 9) | 2–4, 9 |

---

## 6. Future: MHD with a general EOS (Athena++ pattern, W15 / Stage 7)

Athena++ demonstrates that general-EOS MHD is a bounded extension once the
hydro machinery exists: its single `EquationOfState` object serves hydro and
MHD alike, `general_mhd.cpp` implements cons2prim via the `p(ρ,e)` lookup,
and the **fast magnetosonic speed needs only the table sound speed** —
`a² = AsqFromRhoP(ρ,p)` fed into the standard positive-definite fast-speed
formula (`general_mhd.cpp:184-191`); the HLLD solvers gain `GENERAL_EOS`
branches analogous to HLLC's. Nothing about the MHD wave structure requires
γ beyond a².

Cerberus's MHD state is structurally further from this than its hydro state
was (no gas abstraction at all — a bare `Real gamma` member, `MFP_mhd.H:33`),
so the port is a real project, sequenced only when a concrete MHD+tabulated
use case exists:

1. **Memo first** (W12 discipline): one page on what the MHD state actually
   needs from an EOS — for single-fluid MHD it is only `p(ρ,e)`, `e(ρ,p)`,
   `a²(ρ,p)`, i.e. exactly the three surfaces the inverse-map blocks already
   provide; no per-species mixing, no temperature in the flux path.
2. **Closure abstraction with an inert default**: introduce a minimal MHD gas
   closure (three virtuals, mirroring the two added to `HydroGas` in W10 plus
   `p_from_rho_e`); the constant-γ implementation reproduces today's
   behaviour byte-identically — same inert-switch discipline as D2, gated by
   the same baseline suite.
3. **Replace the audited γ sites** (`MFP_mhd.cpp:327,360,403,411,434`) with
   closure calls.
4. **New solver file** `mhd/riemann/MFP_mhd_hlle_general_eos.{H,cpp}`
   mirroring `HLLC_general_eos`: face energies from `e(ρ,p)`, fast speed from
   table `a²` via the positive-definite formula. Existing
   `MFP_mhd_hlle.cpp` untouched.
5. **Validation twin-run**: Brio–Wu shock tube, constant-γ MHD vs tabulated
   closure on a synthetic ideal table — the same isolation logic as the
   Stage-3 Sod twin.
6. **Guard evolution**: the W14 abort relaxes to "abort unless the MHD state
   itself uses a general-EOS closure consistent with the hydro states"; the
   W14 abort message should already point at this section as the lifting
   path.
