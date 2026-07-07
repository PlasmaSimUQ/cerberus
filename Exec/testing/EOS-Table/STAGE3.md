# EOS-Table — Stage 3: first plasma (`effective_gamma` mode)

Detailed plan for Stage 3 of `doc/eos_implementation_plan.md` (work items
W6, W9-Sod, W14). Stage 2 delivered the table engine (`EosTable`) proven in
isolation; Stage 3 wires it into the hydro machinery as a gas model and gets
the first end-to-end shock solution running against it.

Exit criteria (from the main plan):
1. A Sod shock tube run with the tabulated gas on a synthetic ideal-gas
   table matches the analytic gamma-law solution — AND matches a twin run
   using `thermally_perfect` to interpolation-level tolerance;
2. the full existing test suite still matches the W1 `baseline/` record
   (verdict-for-verdict; see the W1 determinism check);
3. the MHD guard (W14) aborts a config that combines an MHD state with a
   tabulated-gas hydro state;
4. the twin-run wall-time ratio (tabulated vs TPG) is measured and within
   budget (see the performance amendment under W6).

---

## W6 — `TabulatedEOS : public HydroGas` (effort: M)

New files: `Source/states/Eulerian/hydro/gas/MFP_tabulated_gas.{H,cpp}`.
Registration is one line, following `MFP_thermally_perfect_gas.cpp:5-7`:
`GetHydroGasFactory().Register("tabulated", HydroGasBuilder<TabulatedEOS>)`.
The builder (`MFP_hydro_gas.H:116`) hands the ctor `(global_idx, gas_table)`
and sets `effective_zero` from the state def — no factory changes needed.

### Lua schema

```lua
gas = {
    type   = 'tabulated',
    table  = 'D_fpeos.eostab',       -- path relative to run dir
    mass   = 2.014,                  -- scalar or per-tracer array (as TPG)
    charge = 0.0,
    names  = {'D'},                  -- optional
    small_temp = 1e-6,               -- code units; -> EosTable floors
    small_dens = 1e-10,
    ttol       = 1e-8,
    max_newton = 100,
}
```

**Composition semantics (v1):** one *thermodynamic* material per state.
`mass`/`charge` follow the TPG array convention so the non-virtual base
helpers (`get_mass_from_cons` etc., used by plasma5/Lorentz/collisions) work
unchanged with tracers present — but tracers are **passive** for
thermodynamics: the single table closes the state regardless of alpha.
Document this in the header; multi-material closure is the Stage-6 mixture
driver.

### Constructor / loading

Read keys, `EosTable::load(path)` + `nondimensionalise()` (Stage-2 W4).
**Ordering check:** the reference quantities (`MFP::rho_ref` etc.,
`MFP.H:224`) must be finalized (`MFP::update_ref`) before any state's
`set_gas` runs in `read_config` (`MFP_config.cpp:112-140`) — verify and, if
not guaranteed, assert in the ctor rather than nondimensionalising garbage.

### Override map

| Virtual (base `MFP_hydro_gas.H`) | Tabulated implementation |
|---|---|
| `get_tag()` | `"tabulated"` (also the hook for the W14 guard) |
| `cons2prim(U,Q)` | e_int = (Eden − ke)/rho; **re inversion** (T_of_e inverse-map seed + Newton polish, plan D5/D8) → T; fill `Prs`=p(rho,T), `Temp`=T, `Gamma`=γₑ, `SpHeat`=cp. Return-bool semantics copied from TPG (caller uses it to detect trouble); non-convergence → floors + `false`. |
| `prim2cons(Q,U)` | **rp inversion** (rho,p) → T,e; `Eden = rho·e + ke`. |
| `define_rho_p_T(Q)` | Close the third of (rho, p, T) from the given two: rt→p, rp→T, **tp→rho** (first live use of the density-axis inversion; ICs like "material at p₀, T₀"). Mirror TPG's convention for which fields count as already-set. |
| `get_temperature_from_cons(U)` | re inversion → T. |
| `get_gamma_from_cons/prim` | **Energy-consistent effective gamma** γₑ = 1 + p/(rho·e). This is what lands in the `Gamma` prim slot, hence what the (unchanged) Riemann solvers consume — face energy exact, solver wave-speed estimate approximate (plan D2). |
| `get_cp_from_cons/prim` | General-EOS identity cp = cv + (T/rho²)·(dp/dT)²/(dp/drho) — reduces to the ideal value on the synthetic table; consumed by the viscous kappa (`MFP_hydro_viscous.cpp:151`). |
| `get_speed_from_cons/prim` | \|u_d\| + cs with the **true table sound speed** → CFL timestep is exact even in `effective_gamma` mode. |
| `write_info(json)` | Table provenance block (source, generator sha, hull) into the run's info output. |

**Newton guess in cons2prim:** primary seed is the `T_of_e` inverse-map
lookup (plan D8 — typically 1–2 polish iterations). Where the map is absent:
analytic guess from the table's ideal-corner gamma (stored as a
header-derived member at load), falling back to the hull midpoint. Incoming
`Q` is never trusted (may be stale or uninitialized); the
warm-start-from-`Temp`-slot refinement only where a caller demonstrably
supplies a valid prior `Q`.

**Performance budget (review amendment):** the granular `HydroGas` getters
(`get_gamma_from_cons`, `get_temperature_from_cons`, `get_cp_from_cons`, …)
each perform an independent (rho,e)→T solve — cheap with the inverse-map
seed, but a viscous/Braginskii run still repeats the *same* solve several
times per cell per step, and for TPG these getters cost nanoseconds so the
interface was never designed to cache. Stage-3 exit therefore includes a
measured gate: **tabulated vs TPG wall time on the twin run ≤ 3×** (replace
with the measured number once known). If it misses, the fix is inversion
de-duplication — getters that accept an already-known T, or a last-solve
memo keyed on (rho,e) — designed then, not built speculatively now.

**Verification reads before coding (review amendment):** three short reads
whose conventions the overrides must copy, none verified yet:
1. TPG's `define_rho_p_T` — *how* it decides which two of (rho,p,T) count as
   "given" (nonzero? sentinel?); a mismatch silently breaks every IC;
2. what each `cons2prim` caller actually does with a `false` return — the
   recent floor/fallback commits touched this area, so the failure contract
   must be read, not assumed;
3. which build configurations define `MFP_PRIM_FLOOR` — the twin run and the
   W1 baseline must be floor-consistent or they diff for non-EOS reasons.

**Deliberately NOT in Stage 3:** `prim_valid`/`cons_valid` are non-virtual
on the base (`MFP_hydro_gas.cpp:152-175`) and keep their p>0 semantics —
fine on the ideal corner where Sod lives. Virtualizing them and the hull-
aware floor logic is Stage 4 (W8/D10), where real cold-curve tables make
p>0 the wrong test.

### Flux-mode resolution update (completes W7)

- With a tabulated gas the `Gamma` prim slot carries γₑ, so **every existing
  solver already runs in `effective_gamma` semantics with zero solver
  changes** — the mode taxonomy (gamma / effective_gamma / eos, plan D2) is
  bookkeeping over the (gas type, `flux` key) pair, not a solver knob.
- Default resolution: `get_tag()=="tabulated"` + no explicit `flux` → the
  usual solver default (e.g. `HLLC`) reading γₑ from the slot. `flux =
  'HLLC_general_eos'` still aborts as unimplemented until Stage 5 (W10),
  where that separate Athena++-pattern solver file arrives and later becomes
  the tabulated default.
- Precedent note: this twin-run design matches Athena++'s `eos_riemann.py`
  regression test (ideal-gas *table* vs native adiabatic EOS on a shock
  tube) — independently converged designs, which is reassuring.

---

## W9-Sod — validation case `Exec/testing/EOS-Sod-Ideal/` (effort: M)

A **sibling** directory (the harness runs one case per directory with a
`check.py`, so it cannot live inside `EOS-Table/`); this plan stays here.
Precedents to copy from: `Double-Rarefaction/` (1D Riemann problem, recent
`check.py` with tolerance machinery) and `Eilmer-Gas/` (existing pattern for
exercising an alternate gas model).

Design — **twin-run comparison**, the sharpest cheap gate available:
1. `problem_definition.lua` parameterized by a `GAS_TYPE` variable
   (settable from the inputs/run script): identical Sod problem with
   `thermally_perfect` (γ=1.4) or `tabulated` + synthetic γ=1.4 table
   (generated by `eos_table_prep.py`, same artifact as Stage-2 tier 1).
   Left/right states chosen well inside the synthetic table's hull.
2. `run`: build (DIM=1, release), generate the table, run **both** variants,
   tee separate logs/plotfiles.
3. `check.py` gates, in order of diagnostic value:
   - **tabulated vs thermally_perfect**: L1 difference per field ≤
     interpolation-level tolerance (table resolution controls this — pick
     the grid so the bound is meaningfully tight; ~1e-4 relative is the
     starting guess, and the measured value becomes the committed threshold
     the first time the test passes — review amendment on guessed gates);
     this isolates *EOS-path* error with the numerics identical;
   - **tabulated vs analytic** exact Riemann solution (adapt the
     Double-Rarefaction checker): L1 ≤ the same tolerance used for the TPG
     run vs analytic (scheme error dominates — both variants must sit at
     the same distance from truth);
   - conservation check: total mass/energy drift ≤ machine-level.
4. Failure modes this catches that the one-zone test cannot: wrong γₑ slot
   filling, prim2cons/cons2prim asymmetry under real time stepping,
   reconstruction of the new slot contents, CFL misbehaviour.

## W14 — MHD guard (effort: S)

1. Comments in `Source/states/Eulerian/mhd/MFP_mhd.H`: at the class header
   and at the `Real gamma` member (`:33`), stating that the MHD state
   assumes a constant-gamma ideal gas throughout (`MFP_mhd.cpp:327,360,403,
   411,434`; `mhd/riemann/MFP_mhd_hlle.cpp`) and is **incompatible with
   `gas={type='tabulated'}` hydro states**; pointer to
   `doc/eos_implementation_plan.md` §1.
2. Config-time abort, placed after the state-construction loop in
   `read_config` (`MFP_config.cpp:112-140`): scan built states; if any is
   an MHD state AND any hydro state's `gas->get_tag() == "tabulated"` →
   `amrex::Abort` naming both offending states and the reason. Using the
   existing virtual `get_tag()` — no casts, no new interface.
3. Verification: manual negative test (two-state Lua config → expect clean
   abort with the message). Optional later: a MemCheck-style harness case
   asserting the nonzero exit; not required for stage exit.
4. Future path: the guard is not forever. `doc/eos_implementation_plan.md`
   §6 plans MHD-with-general-EOS along Athena++'s `general_mhd.cpp` lines
   (table sound speed inside the positive-definite fast-magnetosonic-speed
   formula; a separate `MFP_mhd_hlle_general_eos` solver). The abort message
   should point at that section so the guard documents its own lifting
   condition.

---

## Execution order

1. **W6 core** — ctor/loading + cons2prim/prim2cons/define_rho_p_T first
   (that's enough for ICs + stepping), then the remaining getters.
2. **Flux-mode resolution** — small; unblocks running at all with the
   tabulated gas.
3. **Smoke run** — Sod with tabulated gas, eyeball plotfile vs TPG before
   building the checker.
4. **W9** — twin-run case + `check.py`; iterate W6 until green.
5. **W14** — guard + comments (independent; any time before stage exit).
6. **Regression** — full suite vs `baseline/`; the tabulated model is
   opt-in so this must be clean, verify anyway.

## Stage-3 deliverables checklist

- [x] Verification reads done (define_rho_p_T convention, cons2prim-`false`
      contract, `MFP_PRIM_FLOOR` build gating) — before coding W6
- [x] `MFP_tabulated_gas.{H,cpp}` + factory registration
- [x] Reference-quantity ordering verified/asserted before table load
- [x] Flux-mode defaults: tabulated → existing solver with γₑ slot
      (`effective_gamma` semantics), `flux` override honoured,
      `HLLC_general_eos` still aborts
- [x] `Exec/testing/EOS-Sod-Ideal/` twin-run case green (vs TPG and vs
      analytic; guessed tolerances replaced by committed measured ones),
      conservation clean
- [x] Twin-run wall-time ratio (tabulated/TPG) measured, ≤ budget
- [x] MHD header comments + config-time abort (message pointing at plan §6),
      negative test observed
- [x] Full suite matches `baseline/` (verdict-identity)

### Stage-3 results notes (2026-07-08)

- **Verification-read findings.** (1) `define_rho_p_T`: positivity = given,
  priority (rho,p)→T then (p,T)→rho then (rho,T)→p — mirrored exactly.
  (2) `cons2prim`'s `return prim_valid(Q)` is NOT a graceful failure
  channel: `prim_valid` **aborts** on non-positive rho/p/T, and with
  `MFP_PRIM_FLOOR` (make flag `USE_PRIM_FLOOR`, default TRUE) the floors lift
  values first, so the bool is effectively always true. Inversion
  non-convergence therefore floors + proceeds rather than returning false.
  (3) Ordering confirmed in code: `MFP::update_ref()` (MFP_config.cpp:114)
  runs before state `init_from_lua()` (:164) — the ctor asserts anyway.
- **Units.** The `.eostab` file is CGS, Cerberus references are SI: the
  ctor converts refs to CGS before `nondimensionalise` (rho ×1e-3,
  u ×1e2, p ×10).
- **Twin-run measured numbers** (2048 cells, minmod/HLLC, 4 ranks): twin L1
  per field 1.5e-4 (rho) … 4.1e-4 (T) — committed tol 1e-3; vs exact
  Riemann: tabulated 1.804e-3 vs TPG 1.818e-3 (tabulated marginally
  *closer*); mass/energy drift 0; **advance-time ratio 2.67** (1.27 s vs
  0.48 s) vs budget 3× — the granular-getter inversion cost is real but
  bounded, thanks to the inverse-map seeds; the dedup design stays parked.
- The Sod case needed an explicit `actions = {fluxes = {type='CTU', ...}}`
  block — states define data, actions apply physics; a state without an
  action advances time but never changes (caught on first smoke run).
- W14 negative test observed: two-state config aborts with the message
  naming both states and pointing at plan §6.
