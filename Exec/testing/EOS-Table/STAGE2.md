# EOS-Table — Stage 2: core machinery + config plumbing

Detailed plan for Stage 2 of `doc/eos_implementation_plan.md` (work items
W4, W5, W7). Stage 2 builds the table engine and proves it in isolation —
**no hydro coupling yet** (the `TabulatedEOS` gas model is Stage 3 / W6).
The one-zone test harness lives **in this directory**; its `check.py` +
`run` pair lands here at the end of this stage, at which point
`run_tests.py` starts picking it up automatically.

Exit criteria (from the main plan, amended after review):
1. rt→re→rt and rt→rp→rt round trips converge over the whole table hull,
   measured as the **relative residual in e (resp. p)** ≤ `ttol` — the
   quantities the Godunov update actually conserves. The T error is reported
   but only *monitored*: where cv→0 (degenerate corner) many temperatures
   produce indistinguishable energies, so recovering T tightly from e is
   ill-posed for any algorithm — a T-based criterion there would fail
   forever and tempt a global `ttol` loosening, which is the wrong fix;
2. derivative outputs (cs, dpde, dpdr_e) match centred finite differences
   of the interpolated surface;
3. the full existing test suite matches the W1 baseline
   **verdict-for-verdict**, plus an `fcompare` field-data comparison on one
   case (W1's corrected finding: all physical fields ARE bit-reproducible
   run-to-run; only the `cost` load-balance diagnostic and rank-to-file
   packing differ — so the fcompare gate is "zero error on every field
   except `cost`"). Stage-2 code touches no existing solver files (see
   W7), so this gate should be trivially clean.

---

## W4 — `EosTable` core (effort: M)

New files: `Source/states/Eulerian/hydro/gas/MFP_eos_table.{H,cpp}`.
CPU-only (plan D6): host `std::vector<Real>` storage, no GPU decorations.
No dependency on `HydroGas` — the class must be fully exercisable from the
self-test hook alone.

### Data structures

```c++
struct EosTableView {              // POD, passed by value, non-owning
    int  n_rho, n_T;
    Real lrho_min, dlrho;          // log10-uniform grid (axes reconstructed,
    Real lT_min,   dlT;            //  never stored — matches .eostab spec)
    const Real *p, *e;             // values
    const Real *dpdT, *dpdrho;     // smoothed derivative tables
    const Real *cv,   *dedrho;     //  (outputs only — see inversion note)
    // optional inverse maps (null when absent from the file; plan D5/D8)
    const Real *T_of_e = nullptr;  // T(rho, e) on (lrho x le) grid
    const Real *T_of_p = nullptr;  // T(rho, p) on (lrho x lp) grid
    Real le_min, dle; int n_e;
    Real lp_min, dlp; int n_p;
    const Real *hull;              // REQUIRED block (frozen spec): 1.0 =
                                   // inside source hull, 0.0 = filled cell
                                   // (FPEOS-D: 23% filled — ragged source)
    int idx(int i, int j) const { return i * n_T + j; }   // i=rho, j=T
};

struct EosEval {                   // one evaluation's outputs
    Real p, e, T, rho;
    Real cs, gam1, dpde, dpdr_e;   // Riemann contract
    Real cv, dpdT, dpdrho, dedrho; // for source terms / diagnostics
};

class EosTable {                   // owner: load, condition-check, nondim
    // std::vector<Real> storage for each block; EosTableView view() const;
    // Real small_temp, small_dens, ttol; int max_newton;
    // hull accessors: rho_min/max(), T_min/max() (code units after load)
};
```

### Functions

- `load(path)` — parse `.eostab` per the frozen W2 spec (README.md),
  including the inverse-map blocks + their `le:`/`lp:` axes and the
  **required `hull` block**. Validate hard: magic/version, dims > 1, finite
  values everywhere, cv > 0, block sizes, hull values ∈ {0,1}, hull
  non-empty. Echo the `e_shift` and `conditioning` provenance lines into
  the load report. Malformed input → `amrex::Abort` with file+line
  context. Every rank reads (plan D9).
- `nondimensionalise()` — divide once at load by the Cerberus reference
  quantities (`MFP::rho_ref`, `T_ref`, `prs_ref`, `u_ref` — `MFP.H:224`):
  rho/rho_ref, T/T_ref, P/prs_ref, e/u_ref², cv·T_ref/u_ref². Assert the
  refs are set (config ordering) before any load.
- `eval_rt(view, rho, T, EosEval&)` — bilinear interpolation in
  (log10 rho, log10 T) with hard hull clamping; then the Riemann-contract
  identities with the non-convexity floor:
  `dpde = dpdT/cv`, `dpdr_e = dpdrho − dpde·dedrho`,
  `cs² = max(dpdrho + (T/rho²)·dpdT²/cv, 0)`, `gam1 = rho·cs²/max(p,tiny)`.
- `invert_T_from_e / invert_T_from_p (view, rho, target, T_guess, flag&)` —
  **seed from the inverse map when present** (one bilinear lookup in
  `T_of_e`/`T_of_p` — Athena++-informed, plan D5/D8; lands the start point
  within interpolation error, so the driver typically takes 1–2 polish
  steps), else from `T_guess`. Then guarded Newton (plan D5): **slope =
  analytic d/dT of the bilinear interpolant itself** (NOT the smoothed
  cv/dpdT tables — those disagree with the value surface by construction and
  can stall the iteration near tolerance; they are outputs only). Bracket
  permanently bounded by the hull; any Newton step leaving the bracket
  becomes a bisection step; e(T)-monotonicity NOT assumed. Convergence is
  declared on the **e (resp. p) residual**, not the T increment.
  Non-convergence after `max_newton` → clamp to `small_temp`, set flag,
  never abort mid-step.
- `invert_rho_from_p(view, T, target, rho_guess, flag&)` — same driver,
  density axis (tp mode; IC fills).

### W4 acceptance (informal, pre-harness)

Compiles into the standard exe with zero references from live code paths
except the debug hook (below); no behavior change anywhere.

---

## W5 — One-zone self-test harness in this directory (effort: S–M)

### C++ hook

`eos_table_self_test(table_path, n_sweep)` registered as a free Lua function
exactly like the SDF hooks (`lua.set_function(...)`,
`MFP_ebgeometry_nodeshared.cpp:338` pattern), **gated `#ifdef AMREX_DEBUG`**
consistent with those hooks (release exes will not have it — the `run`
script builds DEBUG). Prints one `PASS`/`FAIL` line per check plus max-error
numbers; `check.py` parses them.

**Dimensional mode (review amendment):** the hook fires during Lua config
execution, and it is *unverified* whether `MFP::rho_ref` etc. are finalized
at that point (Stage 3 flags the same ordering question for the gas-model
ctor, but this hook hits it first). The hook therefore runs the table in
**dimensional (table-unit) mode, skipping `nondimensionalise()`** — every
check here is unit-agnostic (round trips, identities, hull behaviour). The
nondimensional path is first exercised where the references are guaranteed
set: the Stage-3 twin run.

Checks performed by the hook:
1. **Reader integrity** — echo dims/grid/provenance; per-block min/max.
2. **Round trips** — rt→re→rt and rt→rp→rt on an n_sweep × n_sweep grid of
   cells with `hull == 1` plus the hull *boundary* cells. "Hull" means the
   **mask block**, not the grid rectangle — the FPEOS source is ragged and
   23% of grid cells are nearest-value fill (where e is constant along T,
   so re inversion is degenerate by construction; those cells are exercised
   by check 4 instead). Pass criterion:
   max relative **e/p residual** ≤ `ttol`. Also reported: max relative T
   error (monitored, not gated), Newton iteration histogram **with a
   committed ceiling** (commit a hard number the first time the harness runs
   — e.g. fail if any in-hull query exceeds ~2× the synthetic table's p99 —
   rather than letting the harness enshrine whatever it measures),
   inverse-map seed hit rate, bisection-fallback count, non-convergence
   count (must be 0 inside the hull).
3. **Derivative identities** — cs², dpde, dpdr_e vs centred finite
   differences of the *interpolated* surface at cell centres; gam1
   consistency.
4. **Hull behavior** — queries outside the grid rectangle AND queries in
   masked (`hull == 0`) filled cells clamp and flag; no NaN/inf ever. The
   FPEOS tier exercises this naturally (76.9% coverage).
5. **Degenerate-corner stress** — dense sweep of the high-rho/low-T corner
   (small cv → ill-conditioned re inversion). Pass/fail is on the **e
   residual**; T-recovery error is recorded but unbounded by design there
   (see exit criterion 1). Fallback use is expected and counted, failures
   are not.

### Test tables (two-tier)

- **Tier 1 (the gate): synthetic ideal-gas table** generated by
  `eos_table_prep.py` (its `raw` reader path — also the Stage-1 fallback).
  Closed-form truth (p = (γ−1)ρe, cs² = γp/ρ) → *absolute* error thresholds:
  round trips to ttol, cs/dpde to interpolation-order tolerance. This tier
  must pass without FPEOS data existing, so Stage 2 never blocks on data
  acquisition.
- **Tier 2 (informative): conditioned FPEOS deuterium table** from Stage 1
  (`data/D_fpeos.eostab` — **exists, committed**, so this tier is enabled
  from day one): consistency-only checks (round trips on hull cells,
  identities, no non-convergence in-hull). Known stats to expect: 76.9%
  hull coverage; cv is floored in 51 hull cells (0.7%) and every filled
  cell; e spans ~6 decades after the recorded `e_shift`.

### Harness files (added to `EOS-Table/` at the end of this stage)

- `problem_definition.lua` — minimal single-hydro-state problem whose config
  section calls `eos_table_self_test('data/ideal_synthetic.eostab', 64)`
  (and tier 2 if present) at read-config time.
- `onezone.inputs` — AMReX section with `max_step = 0` (init, self-test
  during config, exit; no time stepping).
- `run` — same shape as `Collisions/run`: build
  `make -j8 DIM=1 USE_EB=FALSE AMREX_PARTICLES=FALSE DEBUG=TRUE` in
  `Exec/local/`, generate the synthetic table via `eos_table_prep.py`,
  `mpirun -n 1` the DEBUG exe, tee `run_log.txt`.
- `check.py` — parse `run_log.txt` for the PASS/FAIL lines and thresholds;
  nonzero exit on any failure. **This file's arrival enrols the case in
  `run_tests.py` — do not add it before the harness actually works.**

DIM=1 keeps the build cheapest and the EOS code is dimension-independent.

---

## W7 — Flux-mode plumbing (effort: S) — **revised: no solver edits**

The Athena++ review changed this item's shape. The `eos` flux mode will be a
**separate factory-registered solver**, `riemann/MFP_hllc_general_eos.{H,cpp}`
(Stage 5 / W10, following Athena++'s `hllc.cpp` GENERAL_EOS structure),
selected through the *existing* per-state Lua `flux` key — the same
mechanism that already picks `HLLC` vs `HLLE` vs `AUSMDV`
(`MFP_hydro.cpp:240` → `HydroRiemannSolverBuilder`, `MFP_hydro_riemann.H:27`).
`MFP_hllc.cpp` / `MFP_hydro_hlle.cpp` / `MFP_ausmdv.cpp` are therefore
**never edited**: the `gamma` path stays byte-identical by construction, and
the earlier verbatim-branch scaffolding idea is dropped. The three-mode
taxonomy (gamma / effective_gamma / eos, plan D2) survives as bookkeeping
over the (gas type, `flux` key) pair.

What remains in Stage 2 is small:
- config-time validation in `HydroState::set_flux()`: `flux =
  'HLLC_general_eos'` → clean `amrex::Abort("not yet implemented; see
  doc/eos_implementation_plan.md W10")` until Stage 5 lands;
- nothing else — defaulting-by-gas-type arrives in Stage 3 with the
  tabulated gas model, and the two `HydroGas` face-eval virtuals arrive
  with W10.

### Regression gate (retained)

Rebuild the standard (non-DEBUG) exe, rerun `run_tests.py`, and compare
every case's `check.py` **verdict** against the W1 baseline. Spot-check
with AMReX `fcompare` (`amrex/Tools/Plotfile`, built with plain `make`
there) on at least Double-Rarefaction: **zero error required on every
field except `cost`** (the load-balance timing diagnostic — the one field
that legitimately differs run-to-run; W1 corrected finding). Do NOT use
whole-file checksums — rank-to-file packing permutes identical data.
Since Stage 2 adds only new files plus the abort above, any deviation is a
genuine accident — fix before proceeding.

---

## Execution order

1. **W4** table core (biggest item; everything else consumes it).
2. **W5 hook + synthetic table** — validates W4 the same day it compiles;
   iterate W4 ↔ W5 until tier-1 passes.
3. **W7 config plumbing** — small; independent of W4/W5; gate on the
   verdict diff vs baseline.
4. Add `run` + `check.py` here (enrolling the case in the suite) only once
   tier-1 passes locally.
5. If Stage 1's FPEOS table exists by then, enable tier 2 and review its
   fallback/iteration statistics before calling the stage done.

## Stage-2 deliverables checklist

- [x] `MFP_eos_table.{H,cpp}` — view (incl. inverse-map members), owner,
      reader, eval, 3 seeded inversions
- [x] `eos_table_self_test` debug hook (SDF pattern; dimensional mode)
- [x] Synthetic ideal-gas `.eostab` generation path in `eos_table_prep.py`
      (incl. inverse-map blocks)
- [x] `problem_definition.lua`, `onezone.inputs`, `run`, `check.py` here
- [x] Tier-1 self-test green (2026-07-08: round trips ≤ 2.3e-15, identities
      ≤ 5e-16, fd-vs-blocks 3.3e-10, 0 non-convergences; iteration ceiling
      committed at 60 in check.py — observed max 2–3 on tier 1, 10 on
      tier 2)
- [x] `flux='HLLC_general_eos'` config-time abort (no existing-solver edits)
- [x] Full suite verdict-diff vs `baseline/` clean with Stage-2 code merged
      (+ fcompare field-identity on Double-Rarefaction vs a pre-Stage-2
      plotfile)

### Stage-2 results notes (2026-07-08)

- Tier 2 (FPEOS): all gated checks pass — 1759 in-hull round trips per
  mode, max residual 1.5e-15, worst case 10 iterations, 16–19 bisection
  fallbacks near the ragged hull edge, 0 non-convergences.
  `fd-vs-blocks` reports up to ~1.06 relative difference (PCHIP-conditioned
  slopes vs value-surface secants at fill-boundary cells) — reported, not
  gated, by design; revisit only if Stage-4 robustness work points here.
- Known soft spot: the FPEOS **degenerate-corner sweep finds 0 in-hull
  cells** (the high-rho/low-T corner of that table is entirely filled), so
  the corner stress currently bites only on tables whose hull reaches the
  corner. If a future table has an in-hull degenerate corner, the check
  arms itself automatically; no action needed now.
- The Newton driver's step guard (`|f| <= |slope|*(hi-lo)` before dividing)
  exists because `amrex.fpe_trap_*` turns the overflow from a
  denormal-slope division into a hard abort — found by the tier-2 FPEOS
  table on first run, invisible on synthetic data.
