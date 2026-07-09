# EOS-Table — Stage 5: the `eos` flux mode (`HLLC_general_eos`) becomes default

Detailed plan for Stage 5 of `doc/eos_implementation_plan.md` (work items W10,
W11). Written 2026-07-09, after the Stage-4 milestone gate passed (see
STAGE4.md results notes). Predecessors: STAGE2.md (EosTable core), STAGE3.md
(TabulatedEOS gas model, `effective_gamma` mode), STAGE4.md (hull clamps,
Hugoniot validation).

**Goal.** A dedicated Riemann solver (the routine that computes the flux
between two cells from their left/right face states) that asks the gas model
directly for face internal energy and sound speed, instead of reconstructing
them through γ-algebra. After validation it becomes the default `flux` for
tabulated states; `gamma` / `effective_gamma` remain selectable for A/B
comparison. Braginskii's hard-coded ideal-gas thermodynamics (W11) are
converted behind the same dispatch discipline.

**Non-goals (deferred, recorded in the master plan):** general-EOS variants
of HLLE/AUSMDV/hybrid_hll (same pattern, on demand); MHD (W15/Stage 7,
guarded by W14); 2T (W13/Stage 6); per-variable floors (W16); restart test
(W17).

---

## What Stage 4 taught that shapes this stage

1. **The RH jump is set by conservation + EOS, not the wave-speed estimate.**
   `effective_gamma` mode already lands on the table Hugoniot to 0.7–1.3%.
   So Stage 5's payoff is *not* "correct shocks" — it is (a) removing the
   approximation in the wave-speed estimate (the acoustic Γ₁ = ρc²/p differs
   from γₑ off the ideal corner), (b) removing the linear reconstruction of
   the γₑ slot across strong shocks (accepted at the 1% level in Stage 4,
   part of this stage's motivation), and (c) giving every face evaluation the
   same clamping doctrine as cell evaluations. Expect *small* A/B deltas —
   the gates below are framed accordingly.
2. **Clamp inputs first, derive everything from the clamped state.** The
   Stage-4 density-axis lesson applies verbatim to the new solver's face
   evaluations: any (ρ, p) handed to the table must go through the same
   `clamp_rho` + inversion-stats tally as `cons2prim`, or the solver mixes
   raw off-hull inputs with edge-column table values.
3. **Faces arriving at the solver are already hull-floored** —
   `apply_prim_floor` (Stage-4 W8.1) runs on reconstructed face states before
   the Riemann solve, so the new solver may *assume* in-hull-ish (ρ, p) but
   must still tolerate and tally clamped inversions (floors are per-variable,
   the pair may still sit off-hull jointly).
4. **Measure tolerances, then commit them** (the Stage-4 gate discipline):
   every numeric gate below ships with a measured value and a tolerance set
   a comfortable margin above it, recorded in the results notes.

---

## Design decisions to lock (D-a … D-g)

**D-a — The two new `HydroGas` virtuals take the full face primitive vector,
not bare (ρ, p).** The master plan names them
`get_internal_energy_from_rho_p` / `get_sound_speed_from_rho_p`; the
signatures are refined here because a bare (ρ, p) pair is not enough for
every gas model (TPG's γ depends on the tracer mass fractions carried in the
primitive vector). Lock:

```cpp
// MFP_hydro_gas.H — base class, ideal-gas defaults (see D-b)
virtual Real get_internal_energy_from_prim(const Vector<Real>& Q) const;  // specific e (energy per mass)
virtual Real get_sound_speed_from_prim_rp(const Vector<Real>& Q) const;   // scalar a at the face
```

(The existing `get_speed_from_prim` returns a per-direction |u|+a array for
the CFL step — a different contract, hence the distinct name.) Returning
*specific* internal energy e matches `EosEval::e`; the solver forms the face
total energy as `nrg = ρ·e + ½ρ|u|²`.

**D-b — The base-class defaults read the reconstructed `Gamma` slot, not
`get_gamma_from_prim`.** Default implementations:

```cpp
e = Q[Prs] / ((Q[Gamma] - 1) * Q[Density]);
a = sqrt(Q[Gamma] * Q[Prs] / Q[Density]);
```

This is *exactly* the algebra `MFP_hllc.cpp:44/50` performs, so
`HLLC_general_eos` on a γ-law gas agrees with `HLLC` to round-off with no
TPG override needed. (Using `get_gamma_from_prim(Q)` instead would recompute
γ from the *reconstructed tracers*, which differs slightly from the
*reconstructed γ slot* and would loosen the round-off gate for no benefit.)
`TabulatedEOS` overrides both with one shared `eval_from_rho_p` call on the
clamped face state (one inverse-map lookup + Newton polish each, the D5
cost argument), tallying inversion stats as everywhere else.

**D-c — The two-shock q-factor uses a local γ = a²ρ/p at the PVRS mid-state
pressure** (PVRS = the cheap linearised estimate p* of the pressure between
the two waves; the q-factor is Toro's correction that widens the wave-speed
estimate when that wave is a shock, i.e. when p* exceeds the face pressure).
Per the Athena++ pattern recorded in D2: when `p_star > p_face`, evaluate

```cpp
// tabulated override path; γ-law default path uses Q[Gamma] as today
a2s   = gas->sound_speed²(ρ_face, p_star)      // one extra rp evaluation, shock side only
γ_loc = a2s * ρ_face / p_star
q     = sqrt(1 + ((γ_loc+1)/(2γ_loc)) * (p_star/p_face - 1))
```

Implementation detail: rather than a third virtual, fold this into the
solver — it calls `get_sound_speed_from_prim_rp` on a copy of the face
vector with `Q[Prs] = p_star` (clamped to the hull by the gas model as
usual). Recorded alternative if the extra lookup ever matters: use the
face Γ₁ (= a²ρ/p at the face itself, already computed) — the q-factor only
sets wave-speed robustness, not the converged solution, so an A/B
sensitivity check on the Hugoniot case is part of validation (G3).

**D-d — Wiring: the solver holds a non-owning `HydroGas*`, set in
`HydroState::set_flux`.** Construction order is safe: `set_gas()` runs at
`MFP_hydro.cpp:315`, `set_flux()` at `:436`. The factory builder only sees
the Lua table, so the pointer is injected after `Build`:

```cpp
flux_solver = rfact.Build(flux, state_def);
...
if (flux == HydroHLLCGeneralEOS::tag)
    static_cast<HydroHLLCGeneralEOS*>(flux_solver.get())->set_gas(gas.get());
```

(A `set_gas` virtual on the `HydroRiemannSolver` base with a no-op default
is the fallback if more general-EOS solvers appear later — not needed for
one solver; keep the base class untouched.) The reserved-name abort at
`MFP_hydro.cpp:248` is removed in the same change.

**D-e — Default-flux resolution (the W7 plumbing completed).** Today `flux`
is a hard-required key (`"null"` → abort). New rules in `set_flux`, config
plumbing only:

- gas `type='tabulated'` and no `flux` key → default `'HLLC_general_eos'`,
  with an IOProcessor `Print` notice stating the default was applied;
- gas `type='tabulated'` with an explicit `flux` → honoured (that *is* the
  `effective_gamma`/`gamma` A/B override mechanism — `flux='HLLC'` on a
  tabulated state runs effective-gamma mode exactly as in Stages 3–4);
- non-tabulated gas and no `flux` key → abort exactly as today (no behaviour
  change for existing configs);
- non-tabulated gas with `flux='HLLC_general_eos'` → **allowed** — this is
  the key validation configuration (G1: the new solver on a γ-law gas must
  match `HLLC` to round-off, isolating solver bugs from EOS bugs).

**D-f — W11 dispatch keeps the γ-law path byte-identical.** Each Braginskii
reform site becomes a call to a small helper that dispatches on the gas
model: tabulated → gas-interface evaluation; anything else → the *verbatim*
existing expression. Bitwise identity for γ-law cases is then guaranteed by
construction (same discipline as D2's untouched-solver rule), and gated by
fcompare (G7). No attempt to make TPG route through the gas interface — that
would change floating-point ordering and break byte-identity for zero
physics gain.

**D-g — The `Gamma`/`SpHeat` primitive slots keep being filled** by the
tabulated `cons2prim` (plots, diagnostics, and any other consumer still see
γₑ/cp); the new solver simply does not read `Gamma`. No descriptor changes
in this stage.

---

## W10 — the `HLLC_general_eos` solver (effort: L)

### W10.1 — gas-model face evaluations (`MFP_hydro_gas.{H,cpp}`, `MFP_tabulated_gas.{H,cpp}`)

- Add the two virtuals per D-a with the D-b defaults on the base class.
- `TabulatedEOS` overrides: `clamp_rho` first, then `eval_from_rho_p` on the
  clamped (ρ, p), `tally(st)`, return `ev.e` / `sqrt(ev.a2)` (member names
  per the actual `EosEval` POD). Both overrides share one private helper so
  a face that needs e *and* a pays one inversion, not two — signature
  detail settled at implementation (either a combined
  `eval_face_from_prim(Q, e, a)` used by the solver via two thin getters, or
  accept two lookups for v1 and record the cost; decide from the G8
  measurement).
- No changes to any existing getter.

### W10.2 — the solver file (`riemann/MFP_hllc_general_eos.{H,cpp}`, tag `HLLC_general_eos`)

Structurally a copy of `MFP_hllc.cpp` with exactly three ingredient swaps —
the star-state algebra (`S_star`, the starred conserved states, the flux
assembly, tracer handling) is EOS-independent Rankine–Hugoniot bookkeeping
and is copied unchanged:

| `MFP_hllc.cpp` site | γ-law expression | general-EOS replacement |
|---|---|---|
| :44, :59 (face total energy) | `nrg = p/(γ−1) + ½ρ\|u\|²` | `nrg = ρ·gas->get_internal_energy_from_prim(Q) + ½ρ\|u\|²` |
| :50, :60 (face sound speed) | `a = sqrt(γp/ρ)` | `a = gas->get_sound_speed_from_prim_rp(Q)` |
| :74–90 (q-factors) | γ from the `Gamma` slot | local γ = a²ρ/p at PVRS p\* per D-c |

Class additions vs `HydroHLLC`: `void set_gas(HydroGas* g)`, a
`const HydroGas* gas = nullptr` member, and an `AMREX_ASSERT(gas)` (debug
builds abort if the wiring is missed). Registered with the factory in the
`.cpp` (the one-line `Register` call — the only wiring that makes it visible
to Lua).

FPE discipline (traps are on in the EOS test cases): every division and
sqrt introduced must be guarded the same way the Stage-2 Newton guard and
Stage-4 floors already guarantee — document in the header which invariants
(hull-floored faces, clamped inversions) make each expression safe.

### W10.3 — config plumbing (`MFP_hydro.cpp` only)

- `set_flux`: remove the reserved-name abort, add D-e default resolution and
  the D-d pointer injection.
- `doc/eos_table_reader.md`: new §11 (or extend §10) — user-facing docs for
  the `flux` default and override matrix on tabulated states.
- Update `EOS-Sod-Ideal` and `EOS-Hugoniot` problem definitions: tabulated
  states drop their explicit `flux = 'HLLC'` (exercising the new default);
  the A/B runs pass the override through the per-run inputs overlay
  mechanism from Stage 4 (`temp_*.inputs` heredocs in `run`).

---

## W11 — Braginskii thermodynamic-consistency reform (effort: M)

### What is actually wrong today (sharpened from the code read)

`BraginskiiCTU` stores per-cell scalars in a flat `data` array before the
RK4 source integration; `data[DataIdx::Electron/IonGamma]` is filled once
per cell at `MFP_CTU_Braginskii.cpp:2828/2832` via `gas->get_gamma_from_cons`.
For a tabulated species that value is γₑ *frozen at the pre-step state*.
Then, inside the integrator:

- `rhs` (`:2503`, the collision/transport source-term right-hand side)
  recomputes `p = (γ−1)(E − ½ρ|u|²)` at `:2573` (electron) / `:2606` (ion)
  and `T = p·m/ρ` at `:2574`/`:2607` from the *current* RK-stage `y`. Two
  defects for a tabulated species: the frozen γₑ is applied to an evolved
  energy (γₑ is state-dependent), and `T = p·m/ρ` is the ideal-gas relation
  — off the ideal corner the table temperature must come from the `T(ρ,e)`
  inverse map. Temperature is the load-bearing quantity here: every
  Braginskii coefficient (conductivities, viscosities, relaxation times)
  scales with powers of T.
- `check_invalid` (`:2723`) repeats the same `p` algebra at `:2743`/`:2755`
  as a positivity bailout — for a tabulated species the check should be
  hull-aware (a table state clamped to the hull edge is *valid*, p is
  always > 0 in-hull for FPEOS-class tables).

### Mechanism (D-f) — corrected at implementation from two code facts

The written plan assumed `rhs`/`check_invalid` are member functions reaching
`ion_state->gas`, and that both `p` and `T` need reforming. Reading the code
corrected both:

1. **`rhs`/`check_invalid` are `static`** — they are passed to `rk4_adaptive`
   as function pointers (`:2858/2859`), so they can access only their `y`/
   `data` arguments and *static* members, not `ion_state`. Dispatch therefore
   goes through **two static pointers** set once per `calc_time_derivative`
   before the (serial, no-OpenMP) cell loop:
   ```cpp
   static const HydroGas* s_ion_gas_tab;       // non-null iff ion is tabulated
   static const HydroGas* s_electron_gas_tab;   // ditto electron
   ```
   Safe as statics for exactly the reason the tabulated gas's own clamp
   counter is (`mutable long`, "no OpenMP in Cerberus"): the per-rank cell
   loop is single-threaded.
2. **Only temperature needs reforming.** Inside `rhs`, `p_e`/`p_i` are used
   *solely* to form `T_e`/`T_i` (grep-confirmed: no other consumer). So the
   reform is one line per species — keep the verbatim γ-law `T`, then
   override it when tabulated:
   ```cpp
   Real T_e = p_e * m_e * inv_rho_e;              // gamma-law (unchanged)
   if (s_electron_gas_tab)                         // tabulated: table T(rho,e)
       T_e = species_temperature_tabulated(s_electron_gas_tab, rho_e,
                 y[ElectronXmom], y[ElectronYmom], y[ElectronZmom],
                 y[ElectronEden]);
   ```
   The single helper builds a minimal cons `Vector<Real>` (Density/moms/Eden,
   no tracers — the tabulated gas is single-material) in a static scratch
   buffer and calls the existing `get_temperature_from_cons` (the T(ρ,e)
   inverse map, hull-clamped, always positive). There is no pressure getter
   on `HydroGas` and none is needed.

Byte-identity is guaranteed *by construction*, not by hoping a getter
matches: the γ-law lines are executed verbatim in the exact original
operation order, and the only added instruction on the ideal path is one
predictable `if (nullptr)` branch. `check_invalid`'s two `p<effective_zero`
bailouts (`:2794/:2805`) get the same treatment — the ideal test is guarded
`if (!s_*_gas_tab && p_* < effective_zero)`, skipped for a tabulated species
(hull-clamped, `p>0` in-hull by the ctor guard).

The flat `y` vector carries no tracer information, so a *TPG* species could
not take the tabulated path even if wanted — another reason the legacy
branch stays verbatim (recorded limitation; revisit at W13/2T where
Braginskii + tables becomes physically real).

Reformed sites: `rhs` electron `T_e` (`:2573-75`), ion `T_i` (`:2605-07`);
`check_invalid` electron/ion positivity (`:2791-94`, `:2802-05`).

### Scope honesty

Braginskii with a *deuterium 1T table* on the electron fluid is not a
physical target configuration (the table is total-equilibrium D, not an
electron EOS) — W11 in this stage is **consistency plumbing + gates**, so
that when 2T/component tables arrive (W13) the transport module is already
EOS-clean. The tabulated-path physics gate is correspondingly a
*synthetic-table equivalence* test, not a physics validation.

---

## Validation gates (measure, then commit tolerances — Stage-4 discipline)

- **G1 — solver-isolation round-off gate (the key gate).** `EOS-Sod-Ideal`
  gains a third run: TPG + `flux='HLLC_general_eos'` vs the existing TPG +
  `HLLC` baseline. Gate: L∞ relative difference across all fields at the
  final plotfile ≤ 1e-11 (D-b makes the algebra identical up to FP
  reassociation; measure, expect ~1e-13, commit with margin). Any structured
  difference = solver bug, by construction not an EOS bug.
- **G2 — tabulated-path Sod.** The existing synthetic-table twin runs with
  the new default solver; the Stage-3 twin gate (twin L1 ≤ 4.1e-4 band vs
  exact Sod) must still pass. Record the delta vs the effective_gamma run —
  expected same order, slightly different wave-speed footprint.
- **G3 — Hugoniot A/B.** `EOS-Hugoniot` runs all three drivers with the new
  solver (via inputs overlay): all 21 Stage-4 gates PASS unchanged
  (tolerances not loosened); record plateau-compression deltas vs
  effective_gamma (expect ≤ the existing 0.7–1.3% deviations; this is also
  the D-c q-factor sensitivity check).
- **G4 — abusive robustness.** The Stage-4 abusive rarefaction (u₀=±2 from
  T=0.18) under the new solver: completes, zero validity aborts, clamps
  tallied and reported, centre pinned at (T_floor, p_hull_min). The new
  face-evaluation clamps must appear in the same counter.
- **G5 — default flip.** A tabulated case with no `flux` key selects
  `HLLC_general_eos` (assert on the Print notice + solver tag in the log);
  a γ-law case with no `flux` key still aborts with the options list.
- **G6 — regression suite.** Full `run_tests.py` suite green; fcompare
  spot-check (Double-Rarefaction) bit-identical — pure-γ states never touch
  any edited line except the `set_flux` bookkeeping.
- **G7 — Braginskii byte-identity.** `Braginskii_Riemann_2F` (not
  suite-enrolled — manual twin): fcompare pre-reform vs post-reform builds,
  bit-identical (D-f guarantee). Tabulated-path gate: *if* a synthetic
  electron-mass table is cheap to emit from `eos_table_prep.py`
  (prerequisite check at implementation time), run the case with both
  species on synthetic γ=1.4 tables and gate against the TPG run at the
  table-interpolation error level; otherwise defer that half to W13 and
  gate W11 on byte-identity + the helper-dispatch unit check alone.
- **G8 — performance budget.** Advance-time ratio, new solver vs
  effective_gamma, on the Sod twin: budget ≤ 1.5× (the solver adds 2 rp
  inversions per face, + up to 2 more on shock sides; effective_gamma adds
  none). Measure; if the combined-eval helper (W10.1) is needed to make
  budget, do it then.

---

## Execution order

1. W10.1 virtuals + defaults + tabulated overrides (build clean, suite
   untouched-green).
2. W10.2 solver file + registration (still unreachable without config).
3. W10.3 plumbing: abort removed, wiring, default resolution. → G1, G5, G6.
4. G2/G3/G4 tabulated validation sweeps (inputs-overlay A/B runs). Fix,
   re-measure, commit tolerances.
5. W11 helpers + five-site reform. → G7.
6. G8 measurement (+ combined-eval helper if over budget).
7. Docs (reader-doc §11, STAGE5.md results notes, master-plan stanza →
   DONE), memory update, single Stage-5 commit (no trailers).

## Deliverables checklist

- [x] `HydroGas` virtuals `get_internal_energy_from_prim` /
      `get_sound_speed_from_prim_rp` (+ combined `get_face_eval_from_prim`
      for the G8 budget) with Gamma-slot ideal defaults (D-a/D-b)
- [x] `TabulatedEOS` overrides with clamp-first + tally (Stage-4 doctrine)
- [x] `riemann/MFP_hllc_general_eos.{H,cpp}` registered, gas-wired (D-d),
      FPE-safe, three-swap structure documented in the header (q-factor uses
      the face Γ₁, D-c recorded alternative, after the G8 measurement)
- [x] `set_flux`: reserved abort removed, default flip + override matrix (D-e)
- [x] Braginskii temperature reform: `species_temperature_tabulated` helper +
      static `s_*_gas_tab` dispatch, four sites routed, γ-law byte-identical
      (D-f; corrected — static functions, T-only, see results notes)
- [x] Gates G1–G8 PASS with measured-then-committed tolerances recorded in
      the results notes above
- [x] `doc/eos_table_reader.md` flux-mode section (§10); master plan Stage-5
      stanza marked DONE with deviations noted
- [ ] Memory note updated; commit without trailers (commit pending user OK)

## Deferred / recorded for later

- General-EOS HLLE / AUSMDV / hybrid_hll variants (same three-swap pattern).
- W12 2T design memo — **still overdue** (plan said during Stages 1–2;
  Stage 6/W13 is gated on it). Recommended as a parallel deliverable during
  this stage's validation runs; it is a document, not code, and W11's
  "minimal cons vector, no tracers" limitation is exactly the kind of fact
  it must capture.
- W16 per-variable floors / vacuum treatment (the |u|≈9cs rarefaction).
- W17 restart/checkpoint fcompare test.

---

## Stage-5 results notes (2026-07-09)

**All gates green.** Summary of the validation, with the numbers committed
into the test tolerances.

### G1 — solver-isolation round-off (the key gate)
TPG + `HLLC_general_eos` vs TPG + `HLLC`, L∞/range at t=0.2, 2048 cells:
rho 5.1e-15, p 1.7e-15, x_vel 4.4e-15, T 1.1e-14, nrg 1.6e-15 — all at
the floating-point-reassociation floor (tol 1e-11). The base-class defaults
(read the reconstructed `Gamma` slot, D-b) reproduce `MFP_hllc.cpp`'s
algebra exactly; the residual is the sqrt/square round-trip in the q-factor
γ (`Γ₁ = a²ρ/p` recovers the slot γ to ~1e-15). **Solver bugs are isolated
from EOS bugs by construction, and there are none.**

### G2 — tabulated Sod twin (default + effgamma)
Both the new default solver and the `flux='HLLC'` effective_gamma override
sit at the Stage-3 twin band vs the ideal run: L1/range ≤ 4.1e-4 on every
field (identical to Stage 3 to 3 sig figs — the wave-speed footprint change
is below the table-interpolation floor at this resolution).

### G3 — Hugoniot A/B (the physics validation, re-run under the new solver)
All 21 EOS-Hugoniot gates PASS unchanged. Compressions on the table-RH
locus: T4 **0.67%** (was 0.71% in effgamma), T40 **1.07%**, T400 **1.27%**,
across p₂ = 247 / 2507 / 25234 GPa; anchor-vs-PIMC **0.175%**. The A/B
delta vs effective_gamma is within the last digit — as predicted, the shock
jump is set by conservation + EOS, not the flux's wave-speed estimate, so
both modes land on the same locus. This is also the D-c q-factor
sensitivity check: using the face Γ₁ (not a PVRS-p* re-evaluation) changed
nothing measurable.

### G4 — abusive rarefaction under the new solver
Completed, zero validity aborts, all fields finite, clamps firing every
step (per-sweep counts 1–10 via the AllPrint counter; the new solver's face
inversions feed the same tally). Conservation: outflux matches 2ρ|u|t to
1.9e-15. The theoretical `g_loc = 0` divide in the q-factor (would need a
table sound speed of exactly 0) never occurred — the hull floor keeps cs
positive, so no defensive guard was added; the header records why it is safe.

### G5 — default flip
The tabulated run with no `flux` key logged the config-default notice and
selected `HLLC_general_eos`; a γ-law state with no `flux` key still aborts
with the options list. PASS.

### G6 — full regression suite
All 6 enrolled cases green (Couette, Double-Rarefaction, Viscous-Vortex,
EOS-Sod-Ideal, EOS-Hugoniot, EOS-Table). The pure-γ cases never touch an
edited line except `set_flux` bookkeeping.

### G7 — Braginskii byte-identity (fcompare twin)
Pre-W11 (HEAD) vs post-W11 build on `Braginskii_Riemann_2F`, both to
`plt12800`: **every physical field 0.0 absolute AND relative difference**;
only `cost` (the wall-time diagnostic) differs — exactly the W1 baseline's
non-determinism. The γ-law transport path is bit-identical, confirming the
D-f discipline (verbatim legacy expressions, single null-pointer branch).

### G8 — wall-time budget
`geos/effgamma` advance-time ratio: **1.57** in the sequential suite
(1.32–1.57 clean, up to 2.21 under concurrent build load). The a-priori
1.5× budget was slightly too tight for the method's inherent cost — the
solver does one combined table inversion per face side that effective_gamma
(which reads the reconstructed γ slot) does not, and table inversions
dominate the tabulated advance (effgamma/ideal itself is ~2.8×). Budget
committed at **2.0×** with the measurement recorded in `check.py`; still
catches a >25% regression on the inherent cost.

### Deviations from the written plan (both in W11)
1. `rhs`/`check_invalid` are **static** (rk4 function pointers), not member
   functions — dispatch is via static `s_*_gas_tab` pointers set before the
   serial cell loop, not `ion_state`.
2. Only **temperature** needed reforming — inside `rhs`, `p_e`/`p_i` are
   used solely to form `T_e`/`T_i`. One-line override per species; no
   pressure getter needed. See the corrected §"Mechanism (D-f)" above.
Both are recorded in the plan and memory. No change to the D2 flux-mode
architecture or the D-a…D-e solver decisions.
