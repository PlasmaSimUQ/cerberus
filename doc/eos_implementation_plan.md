# Tabulated EOS for Cerberus — Implementation Plan

Status: planning (2026-07-05; revised 2026-07-07 after a critical review of
the stage plans and a study of Athena++'s general-EOS implementation —
`src/eos/general/` in PrincetonUniversity/athena — which informs D2, D5, D8,
W10 and the new §6 MHD plan). Companion documents:
- `doc/eos_table_reader.md` — standalone documentation of the `EosTable`
  reader and how the code uses it (kept current as stages land)
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
a cold curve). Athena++'s dimensionless-ratio storage (log of `p/e` etc.,
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
| W8 | Floors generalization (D10): gas-model-defined validity; hull clamps; reconstruction-fallback interplay; non-convergence flag path | M | W6 |
| W9 | Test cases in `Exec/testing/`: Sod into the table's ideal-gas corner vs γ-law analytic; single-material strong-shock Hugoniot vs published FPEOS locus | M | W6, W3 |
| W10 | `eos` flux mode = new solver `riemann/MFP_hllc_general_eos.{H,cpp}` (Athena++ `hllc.cpp` GENERAL_EOS pattern, see D2): face (e, a²) from gas-model (ρ,p) evaluations, q-factor from local γ=a²ρ/p at PVRS p*; two new `HydroGas` virtuals with ideal-gas defaults; round-off agreement vs `HLLC` on TPG; A/B vs `effective_gamma` on W9 cases; flip default `flux` for tabulated states | L | W7, W9 |
| W11 | Braginskii consistency: replace the four `p=(γ−1)(E−ke)` / `T=pm/ρ` reform sites (`MFP_CTU_Braginskii.cpp:2573,2574,2606,2743,2755`) with gas-interface calls, behind the same switch | M | W7 |
| W12 | 2T design memo (paper before code): `rho_lookup` data flow — which state's data, at which RK stage — and the electron wave-speed chain rule (∂P_e/∂ρ_e vs table axis ∂P_e/∂ρ_material) | M | parallel |
| W13 | 2T components (SESAME 303/304 split, per-state component binding, 303+304=301 init assertion, exchange-term cv) and, if needed, the mixture/MTMMMT driver | L | W10, W12 |
| W14 | MHD guard: incompatibility comments in `MFP_mhd.H` (class header + `gamma` member); config-time `amrex::Abort` if an MHD state coexists with a `tabulated` hydro gas | S | W6 |
| W15 | MHD with general EOS (future, demand-driven; Athena++ `general_mhd.cpp` pattern — see §6): MHD gas-closure abstraction, table-backed fast magnetosonic speed, `MFP_mhd_hlle_general_eos` solver, guard replaced | L | W10 |

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
Hull-aware floors, survival of a deliberately abusive strong-shock run
(bisection fallback exercised, flagged, non-fatal), Hugoniot overlay against
the published FPEOS locus.
**Milestone gate: everything after this stage may be revised based on what
Stages 1–4 teach.**

**Stage 5 — `eos` flux mode becomes default** *(large; W10, W11)*
The `HLLC_general_eos` solver (Athena++ pattern, D2), validated in three
steps: (1) round-off agreement with `HLLC` on a γ-law gas (isolates solver
bugs from EOS bugs), (2) A/B against `effective_gamma` on the Stage-4 cases,
(3) default `flux` flipped for tabulated states (γ-law states continue to run
the untouched existing solvers). Braginskii reform sites (W11) converted
behind the same mode resolution.

**Stage 6 — Two-temperature** *(large; W12 memo written early, W13 gated on it)*
Per `Cerberus_TabularEOS_Implementation_Plan.txt`, sequenced 2T-first (true
multifluid is the target use case); the mixture driver stays parked until a
concrete need. Note: write the W12 memo during Stages 1–2 regardless — whether
an evaluation ever takes two densities affects the `EosTable` view API.

**Stage 7 — MHD with general EOS** *(large; W15; optional, demand-driven)*
See §6. Not scheduled until a concrete MHD+tabulated use case exists; the
W14 guard stays in force until this stage replaces it.

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
| Tests | `Exec/testing/EOS-Table/` (Stage-1 data + Stage-2 one-zone; see its `README.md`/`STAGE2.md`), `EOS-Sod-Ideal/`, `EOS-Hugoniot/` | 2–4 |

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
