# EOS-Table — Stage 4: robustness + validation gate (`effective_gamma` mode)

Detailed plan for Stage 4 of `doc/eos_implementation_plan.md` (work items
W8 and W9-Hugoniot). Stage 3 wired the tabulated gas end-to-end and proved
it on a *synthetic ideal-gas* table (the Sod twin run). Stage 4 is the
**milestone gate**: it swaps the synthetic table for the **real conditioned
FPEOS deuterium table**, makes the clamp-and-flag machinery load-bearing so
an abusive strong shock survives instead of aborting, and validates the
shocked states against the published FPEOS Hugoniot locus.

> **Stage boundary (from the main plan):** *everything after this stage may
> be revised based on what Stages 1–4 teach.* Stage 5 (the `eos` flux solver)
> is not started until this gate is green and its lessons are folded back.

Exit criteria:
1. A strong-shock run on the **real FPEOS table** completes without aborting;
   cells pushed to the hull edge / into inversion non-convergence are
   **clamped and counted**, not fatal, and the counter is reported.
2. The post-shock states from a swept-strength shock series lie on the
   table's own Hugoniot locus (`qa/hugoniot_D_fpeos.txt`) to scheme
   tolerance — which in turn matched published PIMC to 0.08–0.26 % offline
   in Stage 1.
3. The full existing suite still matches the W1 `baseline/` verdict — the
   tabulated model is opt-in, so this must stay clean.

---

## Reframing after the Stage-3 verification reads

The main plan's D10/W8 is written around cold-curve tables where `p ≤ 0` is
physically valid and the `p > 0` floor is therefore *wrong*. Three findings
from Stage 3 sharpen what Stage 4 actually has to build:

1. **The v1 table has no cold curve.** FPEOS starts at T ≈ 1.35 eV; `p > 0`
   holds across the *entire* hull. So the literal cold-curve case
   (valid state at `p ≤ 0`) does **not** arise in v1. The pressing Stage-4
   robustness is not "p<0 is legal" but **hull-edge clamping and negative
   internal energy produced by the conservative update** (a strong
   rarefaction, or an over/undershoot, can drive `e_int = (Eden−ke)/ρ` below
   the table's `e_min(ρ)`). We build for *that*, and document the cold-curve
   generalization as deferred rather than writing speculative machinery.

2. **`prim_valid`/`cons_valid` abort and are non-virtual**
   (`MFP_hydro_gas.cpp:152,166`). `prim_valid` tests
   `Density/Prs/Temp ≤ 0 → amrex::Abort`; `cons_valid` tests
   `Density/Eden ≤ 0 → amrex::Abort` and **never inspects internal energy**.
   For FPEOS the `p>0` test never wrongly fires, so we do **not** virtualize
   them in v1 — we *prove* by the strong-shock run that they never trip, and
   note the virtualization as the cold-curve lifting point. What `cons_valid`
   misses (a valid `Eden>0` hiding a below-hull `e_int`) is caught one layer
   in, by the `cons2prim` inversion's hull clamp — which is the load-bearing
   guard this stage promotes.

3. **`apply_prim_floor` floors to an absolute `1e-14`**
   (`MFP_hydro_gas.H:88`), the documented dt-collapse pathology
   (ρ→1e-14 while p stays finite → sound speed ~1e6 → timestep collapses).
   The tabulated gas has a *principled physical floor already*: the hull's
   `ρ_min` and `p(ρ, T_min)`. Flooring to the hull edge instead of `1e-14`
   sidesteps that pathology **by construction, for this gas only** — no
   global per-variable-floor rework needed.

---

## W8 — hull-aware floors + load-bearing clamp/flag (effort: M)

Concretely, four pieces. Files: `MFP_tabulated_gas.{H,cpp}`,
`MFP_eos_table.{H,cpp}`, `MFP_hydro_gas.H` (make one method virtual),
`MFP_eulerian.cpp` (reporting hook only).

### W8.1 — `apply_prim_floor` → virtual, hull-physical for the tabulated gas
- Promote `HydroGas::apply_prim_floor` (`MFP_hydro_gas.H:88`) from a
  non-virtual base helper to `virtual`; the base body is unchanged (γ-law
  and Eilmer states keep the exact `effective_zero` floor — zero behaviour
  change, guarded by the baseline suite).
- `TabulatedEOS::apply_prim_floor` override: floor `Density` to the hull's
  `ρ_min` (nondimensionalised at load) rather than `1e-14`, and `Prs` to
  `p(ρ, T_min)` at the floored density (a single forward table eval). This
  keeps the sound speed finite and physical when the floor fires, so the
  timestep does not collapse — the resolution to the open absolute-floor
  issue, scoped to the table.
- The hull `ρ_min`/`T_min` are already members of the `EosTableView`; expose
  a `p_floor(ρ)` helper on the table so both the floor and the ctor's
  provenance dump use one code path.
- **Alternative considered — linear extrapolation off the hull edge**
  (Athena++'s policy: it silently extrapolates its interpolant beyond the
  grid). Instead of snapping a below-hull state to the edge value, continue
  `p`, `e` linearly along the edge derivatives (`dpdrho`, `dpdT` are already
  stored) for a bounded distance outside the hull. Pros: no artificial
  plateau at the edge, smoother behaviour for transients that graze the hull
  and come back. Cons: extrapolated states are *not* from the physics data
  (below a cold curve this can be badly wrong — the reason plan D5 chose
  clamp-and-flag), needs its own positivity guards (`extrapolated p`/`cv`
  can go negative), and hides hull-sizing problems the clamp counter would
  expose. **Decision for v1: clamp-and-flag stays the default**; linear
  extrapolation is recorded here as the candidate refinement if the Hugoniot
  runs show the edge plateau itself causing artefacts (e.g. a spurious wave
  reflected off the clamp discontinuity).

### W8.2 — hull-aware internal-energy clamp in `cons2prim` (the real guard)
- Today `cons2prim` computes `e_int = (Eden − ke)/ρ` then inverts for T.
  When `e_int` lands below `e_min(ρ)` (bottom of the hull) or the inversion
  does not converge, Stage 3 already floors-and-proceeds — Stage 4 makes
  that path *correct and observable*.
- **Plan-check against the code (2026-07-08)** — two findings that sharpen
  the original item:
  1. The inversion driver is *already* safe: `invert_1d`
     (`MFP_eos_table.cpp`) brackets Newton in `[T_min, T_max]` on the density
     column, and an unattainable `e` target clamps to the nearer end with
     `stats.flag = 1` (never aborts). So no pre-clamp is needed to protect
     Newton itself.
  2. **But the flag is discarded and γₑ can overflow.** `eval_from_rho_e`
     (`MFP_tabulated_gas.cpp:80`) creates a local `EosInvertStats` and drops
     it — the clamp signal never reaches `cons2prim`. Worse, the γₑ formula
     (`:121`) divides by `max(ρ·e_int, DBL_MIN)` using the **physical**
     `e_int`: a strong rarefaction that drives `e_int` negative produces
     γₑ ≈ 10³⁰⁸ in the `Gamma` slot the Riemann solvers consume — a latent
     overflow the Sod case never reaches but the abusive run will.
- Work, therefore:
  - clamp `e_int` to `[e_min(ρ), e_max(ρ)]` on the current density column
    **in `cons2prim`, before the inversion call** — not to protect Newton
    (already bracketed) but so that **γₑ, and every quantity derived from
    `e_int`, uses the same clamped value the returned (T, p) correspond
    to** (self-consistent primitive set; kills the overflow);
  - plumb `EosInvertStats` out of `eval_from_rho_e`/`eval_from_rho_p`
    (out-parameter or return) so `cons2prim` sees `flag ≠ 0`;
  - record the clamp/non-convergence in the W8.3 counter.
- This is where a strong shock's expansion fan or a conservative overshoot
  is caught — precisely the region `cons_valid`'s `Eden>0` test cannot see.

### W8.3 — make the `clamped` flag load-bearing (counter + report)
- The `clamped`/non-convergence flag is currently *set but consumed by
  nothing* (and, per the W8.2 check, discarded inside `eval_from_rho_e`
  before it even leaves the gas model). Wire it like `apply_prim_floor`'s
  existing reporting: a per-box counter accumulated across `cons2prim`
  calls, reduced and printed once per step under `verbosity ≥ 2` (e.g.
  `"[tabulated] clamped N cells to hull this step"`). No per-cell prints
  (hot loop).

#### Clamped-cell plot field (diagnostic; build with W9, not before)
Beyond the scalar counter ("*how many* cells were clamped this step"), a
**derived plotfile field** answers "*where and which way*" — which is what
actually guides debugging when a Hugoniot point sits off the locus:

- **What it shows.** Per cell, a small integer code from the last
  `cons2prim` on that cell: `0` in-hull converged; `1` e below the column
  minimum (clamped to `T_min` — rarefaction/undershoot side); `2` e above
  the column maximum (clamped to `T_max` — the driver is pushing past the
  table top); `3` Newton exhausted (`flag=2`, best-bracket value used);
  `+10` if the *density* itself was outside the hull (`locate` clamped the
  ρ index). Distinguishing the *direction* matters: code 1 says "extend the
  table downward / soften the IC", code 2 says "the sweep outran the table",
  code 3 says "conditioning problem in that region".
- **How it's plumbed (mechanism exists).** Cerberus already builds derived
  plot outputs from primitives via each state's plot-function machinery
  (`get_plot_output` in the hydro state; the Lua `plot` block selects
  fields). The natural v1 implementation is *stateless*: a plot-time
  function that **re-runs the classification** — recompute `e_int` from the
  cell's conserved state and compare against the hull's `e_min(ρ)/e_max(ρ)`
  column bounds (one bilinear read, no inversion needed). That avoids
  persisting a new cell-flag array through the update/regrid/checkpoint
  machinery for a purely diagnostic quantity; cost is paid only at plot
  cadence, only when the field is requested in the Lua `plot` block.
- **Limitation to note in the doc:** the stateless recompute reflects the
  state *at plot time*, so a transient clamp that fired mid-step and
  recovered is visible only in the step counter, not the field. Acceptable
  for the intended use (locating persistent off-hull regions in the
  Hugoniot sweep); a persisted flag array is the upgrade if transient
  events ever need mapping.
- **Gate tie-in:** the W9 abusive run's `check.py` can then assert not just
  `counter > 0` but that clamped cells are confined to the expected region
  (e.g. the expansion fan), turning the diagnostic into a test oracle.

### W8.4 — `get_positive_prim` / validity: verify-not-virtualize for v1
- `get_positive_prim` (`MFP_hydro.H:56`, the reconstruction-fallback
  bad-face set) hard-codes `{Prs, effective_zero}`. For FPEOS `p>0` holds
  hull-wide, so this never false-fires; **leave it unchanged in v1**.
- Document at the override (and in the reader doc) that a *cold-curve* table
  would require this set — and `prim_valid`/`cons_valid` — to test
  `e_int`/hull-membership instead of `p>0`; that virtualization is the
  cold-curve lifting task, not v1 work.
- **Verification gate:** the strong-shock run (W8 exit) must complete with
  **zero** `prim_valid`/`cons_valid` aborts. If either trips, the p>0
  assumption has met a real counterexample and virtualization is pulled
  forward — a decision point, not a silent workaround.

### Abusive-run robustness gate
A deliberately violent case (very high pressure-ratio shock tube, and/or a
double rarefaction that drives cells toward the hull's low-`e` corner) that
*would* off-hull without the guards. Assert: run completes to `t_end`;
clamp counter > 0 (the guard actually fired); no abort; conserved totals
still drift only at machine level. This exercises W8.2/W8.3 and the
Stage-2 Newton bisection fallback at the degenerate hull corner.

---

## W9-Hugoniot — validation case `Exec/testing/EOS-Hugoniot/` (effort: M)

A new sibling directory (one case per `check.py`, as with `EOS-Sod-Ideal/`).

### The physics gate, and why it works in `effective_gamma` mode
The Rankine–Hugoniot jump across a shock is fixed by **conservation of
mass/momentum/energy + the EOS** — not by the flux function's internal
wave-speed estimate. A conservative Godunov scheme enforces the jump
regardless of whether the solver used the exact sound speed or the
`effective_gamma` approximation; the approximation only affects the *width*
and robustness of the captured front, not the *state* it connects to. So the
Hugoniot locus can be validated now, in `effective_gamma` mode, **before**
the Stage-5 `eos` solver exists — an important sequencing point.

### Driving and measuring the locus
- **Drive:** a shock tube with a strong high-pressure driver, parameterized
  by a driver-strength variable swept over N values (like `GAS_TYPE` was
  swept in Stage 3). Each run produces one shocked state → one point on the
  principal Hugoniot. Uses only existing transmissive/reflecting BCs — no
  new inflow/piston BC. (A reflected-shock off a wall, to reach higher
  compression ratios, is an optional extension if the driver series does not
  reach the interesting FPEOS compression regime.)
- **Initial state:** deuterium at a chosen `(ρ0, T0)` inside the hull,
  matching the reference point of the published Hugoniot
  (`data/raw/hugoniot_MC2000_PRL85_1890.txt`, transcribed in Stage 1).
- **Measure:** from the last plotfile, extract the uniform post-shock region
  (ρ2, p2, u2) and the shock speed us (front position / time, or
  RH from the measured jump). One point per driver strength.

### `check.py` gates
1. **On-Hugoniot:** each measured post-shock `(ρ2, p2)` lies on the table's
   own Hugoniot `qa/hugoniot_D_fpeos.txt` (generated offline in Stage 1 from
   the *same* table) within scheme tolerance — this isolates solver/scheme
   error from EOS error, exactly as the Stage-3 twin isolated EOS-path error.
   Tolerance MEASURED-and-committed on first pass (no guessed gates).
2. **Physics anchor:** that table Hugoniot already matches published PIMC to
   0.08–0.26 % (Stage-1 result) — so passing (1) chains the running code to
   the literature. `check.py` re-asserts the table-vs-PIMC agreement as a
   guard against the table artifact drifting.
3. **Conservation:** total mass/energy drift ≤ machine level (as Stage 3).
4. **Robustness:** the strongest driver in the sweep doubles as the W8
   abusive-run gate (clamp counter reported, no abort).

### Precedent to copy
`EOS-Sod-Ideal/` (twin-run scaffolding, `run` structure, `check.py`
tolerance machinery, the exact-Riemann helper pattern) and
`Double-Rarefaction/` (strong-rarefaction robustness).

---

## Deferred / later goals (recorded here, scheduled later)

Tracked in the main plan's work-item table (W16, W17); neither gates Stage 4.

- **W16 — separate `rho_floor` / `p_floor` (global, all gas models).** The
  Stage-4 hull floor fixes the absolute-floor dt-collapse pathology *for the
  tabulated gas only*; γ-law states keep the single `effective_zero` for
  both density and pressure, and the documented failure (ρ floored to 1e-14
  while p stays finite → cs ~ √(p/ρ) explodes → dt collapses) remains
  reachable there. The general fix is per-variable floors: independent Lua
  keys `rho_floor` / `p_floor` (defaulting to `effective_zero` so absent
  keys change nothing), consumed by `apply_prim_floor`, the `cons2prim`
  floor blocks, and the `get_speed_from_*` guards — plus a decision on
  vacuum-cell treatment (zero the velocity when the density floor fires?).
  Touches every gas model and the MHD state, so it is its own change with
  its own regression pass, not a Stage-4 rider.
- **W17 — checkpoint/restart test with the tabulated gas.** Untested today.
  Expected to work: the gas model adds no new `StateData` (the table is
  reloaded from the `.eostab` file by every rank at config time, same as the
  Lua script itself), so a restart exercises only existing machinery. The
  test: extend a case's `run` script to stop at step N, restart from the
  checkpoint (`amr.restart`), run to `t_end`, and `fcompare` the final
  plotfile against an uninterrupted run — bit-identity expected (the W1
  determinism finding makes this a sharp gate). Worth doing before any long
  Hugoniot-style production runs rely on restart; natural home is a variant
  of the EOS-Hugoniot `run` once that case exists.

## Execution order

1. **W8.1** (virtual `apply_prim_floor` + hull floor) — smallest, unblocks
   running the real table without the 1e-14 pathology.
2. **W8.2/W8.3** (hull-e clamp + counter/report) — the load-bearing guard.
3. **Smoke:** one strong FPEOS shock, eyeball the profile + clamp report.
4. **W9-Hugoniot** case + `check.py`; sweep driver strengths; iterate W8.
5. **Abusive-run gate** (W8.4 verification: zero aborts, counter > 0).
6. **Regression:** full suite vs `baseline/` (verdict-identity + one
   fcompare spot-check).

## Stage-4 deliverables checklist

- [x] `HydroGas::apply_prim_floor` made `virtual`; base body byte-unchanged
      (destructor also made virtual — pre-existing UB, gas models are
      deleted through the base pointer)
- [x] `TabulatedEOS::apply_prim_floor` floors to hull `ρ_min` / `p_hull_min`
      (global in-hull pressure minimum, cached at load — cheaper than the
      per-column `p(ρ,T_min)` eval and sufficient; see results notes)
- [x] Off-hull inputs clamped self-consistently in `cons2prim` (`clamp_rho`
      up front + flagged-inversion `e_eff`; supersedes the planned
      pre-clamp of `e_int` — see results notes on why the density axis
      turned out to be the load-bearing clamp)
- [x] `clamped` counter accumulated + verbosity-gated per-step report
      (via `AllPrint` — see results notes)
- [x] Cold-curve generalization documented at `get_positive_prim` /
      `prim_valid` / `cons_valid` (deferred, not built); ctor aborts on a
      table with a non-positive in-hull pressure minimum
- [x] `Exec/testing/EOS-Hugoniot/` on the **real FPEOS table**: swept-driver
      shock series, post-shock states on the table-predicted RH locus
      (measured tols committed), conservation machine-zero
- [x] Abusive rarefaction run completes: clamp counter 10,966 > 0, **zero**
      `prim_valid`/`cons_valid` aborts
- [x] Full suite (6 cases) green; fcompare spot-check bit-identical

### Stage-4 results notes (2026-07-08)

- **The density axis, not the energy axis, was the load-bearing clamp.**
  The plan (and W8.2) focused on `e_int` falling below the hull; in
  practice `invert_1d` already clamps unattainable-`e` targets safely. The
  two REAL defects found and fixed:
  1. `eval_from_rho_e` discarded its `EosInvertStats`, and the γₑ formula
     divided by the *physical* `ρ·e_int` — the predicted ~1e308 `Gamma`
     overflow for negative `e_int` (fixed: flagged inversions use the
     table-consistent `ev.e`).
  2. **An off-hull density was silently clamped inside `locate()`** with no
     flag at all: `cons2prim` then mixed the raw ρ (e.g. 4.3e-3, below the
     1.15e-2 hull floor) with edge-column table values — inconsistent face
     states that drove `HydroHLLC` to an FPE (caught at `MFP_hllc.cpp:44`
     and `:143` by the trap-enabled abusive run). Fix: `clamp_rho()` at
     every gas-model entry point, BEFORE anything is derived, so
     velocities, `e_int`, γₑ, the evaluation column and `Q[Density]` all
     describe the same in-hull state. This also bounds `u = mx/ρ` and so
     prevents the absolute-floor dt-collapse pathology by construction.
- **`Print()` vs `AllPrint()`**: the clamp report initially never appeared
  — clamping cells belong to interior ranks and `amrex::Print()` emits from
  the IO rank only. The report uses `AllPrint()`. (The pre-existing face-
  floor report in `calc_fluxes` has the same rank-0-only behaviour; left
  untouched — noted for a possible cleanup.)
- **Abusive-run design**: a vacuum-forming double rarefaction (|u|≈9cs) is
  W16 territory (per-variable floors / vacuum velocity treatment) and still
  dies in dt-collapse for ANY gas model. The Stage-4 gate uses |u|=2·cs·3
  from T=0.18 (15% above the 15625 K floor): the centre pins at exactly
  (T_floor, p_hull_min) with ρ falling below the density hull — 10,966
  clamped evaluations, run completes, all fields finite, and the mass loss
  through the outflow boundaries matches the analytic boundary-flux
  prediction 2ρ|u|t to 1.9e-15.
- **Hugoniot numbers** (1024 cells, minmod/HLLC/RK2, 4 ranks): measured
  compressions 3.428/4.016/3.961 vs table-RH predictions 3.452/4.059/4.012
  at p₂ = 247/2507/25234 GPa → rel errs 0.71/1.07/1.26% (committed tol 2%;
  captured-shock plateau-median level). Pre-shock state closes through the
  table to 3.5e-7. Physics anchor (offline locus vs PIMC ≥250 kK): 0.175%.
  Conservation: ≤3.7e-12 with all waves in-domain. The `effective_gamma`
  sequencing argument held: the RH jump is set by conservation + EOS, and
  the solver's approximate wave speeds did not bias the locus at the 1%
  level.
- **Stop-time discipline**: the driver rarefaction head is much faster than
  γ-law intuition suggests (T400: >37 code units) — stop times were tuned
  empirically until first-vs-last sums hit machine zero (a leaked
  8.3e-9 at T4/0.07 was the rarefaction *head* just grazing the wall).
- `check.py` gates: 21, all PASS; W14-style negative checks not needed
  here (covered by Stage 3).
