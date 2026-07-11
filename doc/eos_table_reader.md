# The tabulated-EOS reader (`EosTable`)

Standalone documentation for the EOS table reader and how it is used in the
code. Status: covers the machinery as of Stage 2 of
`doc/eos_implementation_plan.md`; sections marked *(future)* describe how
later stages will consume it, and this document is updated as they land.

Related documents:
- `doc/eos_implementation_plan.md` — the engineering plan (design decisions
  referenced here as D1–D10, work items W1–W15).
- `Exec/testing/EOS-Table/README.md` — the frozen `.eostab` file-format
  specification and the Stage-1 data-conditioning results.

---

## 1. What it is, in one paragraph

Cerberus normally closes the fluid equations with an ideal-gas law
(p = (γ−1)ρe, one line of algebra). For warm-dense-matter problems the real
material behaviour is only available as *data* — a table of pressure and
energy over a grid of density and temperature, produced offline from
first-principles simulation results. `EosTable` is the class that reads such
a table from disk, holds it in memory, and answers the two kinds of question
a hydro code asks an equation of state:

1. **Forward:** "at this density and temperature, what are p, e, the sound
   speed, and the thermodynamic derivatives?" — one interpolation.
2. **Inverse:** "I know the density and the internal energy (or pressure) —
   what temperature is the material at?" — a small root-finding problem
   (finding where a curve crosses a target value), because the table is
   stored with temperature as an input, not an output.

## 2. Where the pieces live

| Piece | Location |
|---|---|
| Reader/evaluator (C++) | `Source/states/Eulerian/hydro/gas/MFP_eos_table.{H,cpp}` |
| Table *maker* (Python, offline) | `Exec/python_analysis/eos_tools/` package (`eos_table_prep.py` is its CLI shim; plan: `doc/eos_creation_plan.md`) |
| Format spec + conditioned tables | `Exec/testing/EOS-Table/` (`README.md`, `data/*.eostab`) |
| Self-test harness | `Exec/testing/EOS-Table/{problem_definition.lua, onezone.inputs, run, check.py}` |
| Lua registration | `Source/MFP_config.cpp` (`EosTable::register_with_lua`) |

The division of labour is deliberate (plan D4): **all data curation happens
offline in Python** — regridding scattered source points, smoothing
derivatives, flooring cv, building inverse maps — where it can be inspected,
plotted and version-controlled. The C++ reader parses exactly one canonical
format and does *no* conditioning; if a file is malformed it stops the run
immediately rather than computing with suspect data.

## 3. The `.eostab` file (what the reader reads)

A plain-text file (ASCII — human-readable, diff-able in git; tables are only
~1 MB). Authoritative spec: `Exec/testing/EOS-Table/README.md` §W2, mirrored
in the docstring of `eos_table_prep.py`. Structure:

```
EOSTAB 1                       <- magic line (file-type check) + version
material: D                    <- free-form provenance key/value lines
source: FPEOS ...                 (where the data came from, when, how)
units: cgs
e_shift: 1.23e+12              <- constant added to e so it is positive
conditioning: cv_floor=...     <- record of what the offline tool changed
grid: n_rho=96 n_T=96          <- grid dimensions
lrho: -2.7 3.2                 <- log10 axis ranges (axes are reconstructed
lT: 4.19 7.81                     from range+count, never stored as arrays)
le: ... n_e=96                 <- axes for the inverse maps
lp: ... n_p=96
block: p                       <- then the data blocks, each a header line
1.29e+09 4.71e+09 ...             followed by n_rho*n_T numbers
block: e
...
```

Nine blocks, all on the same (density, temperature) grid unless noted:

| Block | Meaning | Used for |
|---|---|---|
| `p`, `e` | pressure, specific internal energy | the value surfaces everything interpolates |
| `dpdT`, `dpdrho`, `cv`, `dedrho` | pre-smoothed derivative tables | sound speed and the other derivative outputs — **outputs only, never the Newton slope** (see §6) |
| `T_of_e`, `T_of_p` | pre-inverted maps T(ρ,e), T(ρ,p) on their own energy/pressure axes | starting guesses ("seeds") for the inversions |
| `hull` | 1.0 = real source data here, 0.0 = cell filled by the tool | flagging queries that land outside the trustworthy region |

Two format details worth understanding:

- **Everything is stored in ordinary CGS units** (grams, centimetres,
  Kelvin). Cerberus works internally in nondimensional units (all quantities
  scaled by reference values so they are O(1) numbers); the conversion
  happens **once at load time** (§5), so no evaluation ever multiplies by a
  unit factor.
- **The hull** exists because real source data is ragged: the FPEOS dataset
  has 5–16 temperature points per density, so when the tool regrids onto a
  rectangle, ~23% of cells have no data behind them. Those cells hold
  copies of the nearest real value (so arithmetic never sees garbage) and a
  0 in the hull mask (so the reader can *tell you* the answer is
  extrapolated).

## 4. The C++ objects and their lifecycle

Three small types (all in `MFP_eos_table.H`):

- **`EosTable`** — the owner. Holds the data in `std::vector<Real>` arrays
  (plain resizable memory on the CPU — plan D6 is CPU-only for v1). You
  `load()` it once, optionally `nondimensionalise()` once, then ask it for a
  view.
- **`EosTableView`** — a small struct of raw pointers + grid constants
  (POD, "plain old data": something that can be copied byte-for-byte with no
  hidden machinery). It does not own anything; it is the cheap handle that
  evaluation functions take. This owner/view split is the structure a future
  GPU port needs — a view can be copied to a device kernel wholesale —
  without any GPU code existing today.
- **`EosEval`** — the result of one evaluation: ρ, T, p, e, sound speed
  `cs`, `gam1` (Γ₁ = ρc²/p, the "effective gamma" the wave-speed estimates
  want), the derivative combinations the Riemann solvers need (`dpde`,
  `dpdr_e`), the raw derivatives, and a `clamped` flag (§7).

Typical use (this is exactly what the self-test does):

```c++
EosTable tab;
tab.load("data/D_fpeos.eostab");     // parse + validate; aborts on bad input
// tab.nondimensionalise(rho_ref, T_ref, prs_ref, u_ref);  // Stage 3 does this
EosTableView v = tab.view();

EosEval out;
EosTable::eval_rt(v, rho, T, out);   // forward: (rho, T) -> everything

EosInvertStats st;
Real T2 = EosTable::invert_T_from_e( // inverse: (rho, e) -> T
    v, rho, e_target, /*T_guess=*/-1.0, tab.ttol, tab.max_newton, st);
```

The evaluation functions are `static` (they belong to the class but need no
particular object — they act only on the view you pass them), which keeps
them free of hidden state and trivially testable.

## 5. Loading and units

`load(path)` parses the header, reads every block, and **hard-validates**:
magic line, grid dimensions > 1, every value finite, cv > 0 everywhere, hull
values exactly 0 or 1, required blocks present and full-length. Any failure
calls `amrex::Abort` (immediate stop of the whole run) with the file name
and the offending item — the philosophy is that a bad table should be
impossible to run with, not a warning you can scroll past. On success it
prints a short load report echoing the provenance lines (`material`,
`source`, `e_shift`, `conditioning`), so the log records exactly which data
produced the results.

Every MPI rank (each parallel process of the simulation) reads the file
independently — consistent with how Cerberus runs its Lua configuration on
every rank (plan D9). At ~1 MB per table this is harmless; a shared-memory
scheme like the STL/SDF work is a drop-in later if tables grow.

`nondimensionalise(rho_ref, T_ref, prs_ref, u_ref)` divides every block by
the appropriate combination of reference quantities and shifts the
logarithmic axes (dividing a quantity by a constant is a constant *shift* in
log space, so the uniform grid stays uniform). It may be called exactly once
— a second call aborts, which prevents double-scaling bugs. Note the
self-test never calls it: it runs in "dimensional mode" because it executes
during configuration reading, before the reference quantities are
guaranteed to exist.

One subtlety recorded in the file itself: **`e_shift`**. First-principles
energies are negative at low temperature (bound states), but much of the
hydro machinery implicitly prefers e > 0. The zero point of internal energy
is physically arbitrary, so the offline tool adds a constant to make the
whole table positive and records that constant in the header. It cancels
identically in anything the dynamics depends on — *provided every quantity
comes from the same table*, which is why the shift lives in the table file
and not in user configuration.

## 6. Forward evaluation: `eval_rt`

Interpolation is **bilinear in (log₁₀ρ, log₁₀T)** — the four grid values
surrounding the query point are blended linearly in each direction. Log
axes because the data spans ~6 decades in density; bilinear because it is
cheap, monotone within a cell (it cannot invent new maxima), and its slope
is available in closed form (needed by the inversions). The known trade-off
is that the interpolated surface is C⁰ (continuous values but slightly
kinked slopes at cell edges); a smoother bicubic scheme is a possible later
upgrade on the same file format.

After interpolating the six stored surfaces, the Riemann-contract outputs
are combined algebraically (plan D5):

```
dpde   = dpdT / cv                        (pressure response to energy)
dpdr_e = dpdrho − dpde · dedrho           (pressure response to density at fixed e)
cs²    = max(dpdrho + (T/ρ²)·dpdT²/cv, 0) (sound speed; floored at zero so a
                                           badly-conditioned corner degrades
                                           accuracy rather than producing an
                                           imaginary wave speed)
gam1   = ρ·cs²/p
```

## 7. Inverse evaluation: the guarded, seeded Newton

`invert_T_from_e` / `invert_T_from_p` (recover T along a fixed-density line)
and `invert_rho_from_p` (recover ρ along a fixed-temperature line — used for
initial conditions specified as "this pressure at this temperature") all run
the same driver:

1. **Seed.** Look up the pre-inverted map (`T_of_e`/`T_of_p`) — one
   bilinear read that lands within interpolation error of the answer. This
   is the idea adopted from Athena++'s EOS tables, and it is why inversions
   cost ~2–3 iterations in practice instead of ~tens. Without a map, the
   caller's guess or the axis midpoint is used.
2. **Bracket.** The root is confined to the table's temperature (or
   density) range. If the target value is not attainable anywhere on that
   line, there is no root: the driver returns the nearer endpoint and sets
   `flag = 1` — the caller decides what a clamped answer means.
3. **Newton iteration** (root refinement using the local slope), where the
   slope is the **analytic derivative of the bilinear interpolant itself**,
   *not* the smoothed `cv`/`dpdT` blocks. This matters: the smoothed blocks
   were conditioned offline and disagree slightly with the value surface by
   construction, so iterating with them can stall just above the tolerance;
   the interpolant's own slope is exactly consistent with the values being
   matched.
4. **Safeguards.** Any Newton step that would leave the bracket becomes a
   bisection step (halving the bracket — slow but unconditionally safe), and
   the step is size-checked *before* the division (`|f| ≤ |slope|·(hi−lo)`)
   so a near-zero slope can never cause a floating-point overflow — with
   `amrex.fpe_trap_*` enabled (a debugging option that turns floating-point
   anomalies into immediate crashes) that overflow would abort the run.
   This guard was added after the real FPEOS table triggered exactly that
   at a filled-cell boundary.
5. **Convergence is declared on the residual in e (or p)** — the quantity
   the conservative update actually conserves — **not** on the temperature
   increment. In degenerate matter cv → small means many temperatures give
   nearly the same energy: T is genuinely ill-determined there, but the
   energy residual can still be driven to tolerance, which is what the
   dynamics needs.
6. **Failure is a flag, never an abort.** `flag = 2` after `max_newton`
   iterations returns the best bracketed value; mid-timestep the correct
   response is a floor and a diagnostic, not a crash (plan D5). The
   `EosInvertStats` struct reports iterations, bisection count, whether the
   seed was used, and the flag — the self-test aggregates these into its
   pass/fail lines.

Measured behaviour (Stage-2 harness, see `EOS-Table/STAGE2.md`): synthetic
ideal-gas table — 2–3 iterations, residuals ≤ 2.3e-15; FPEOS deuterium —
worst case 10 iterations with 16–19 bisection fallbacks near the ragged hull
edge, zero non-convergences over 1759 in-hull points per mode.

## 8. The `clamped` flag and the hull

`eval_rt` sets `out.clamped = true` in two situations: the query point lies
outside the table's (ρ,T) rectangle (coordinates are clamped to the edge —
the table never extrapolates, plan D5), or the containing cell is marked 0
in the hull mask (the value returned is arithmetic on *filled* placeholder
data, not physics).

Since Stage 4 (W8) the clamping is **load-bearing** in the gas model, with
one important subtlety learned the hard way: the low-level `locate()` clamps
an off-grid *density* silently (no flag — only an unattainable inversion
*target* sets `EosInvertStats::flag`). The gas model therefore clamps the
density itself (`TabulatedEOS::clamp_rho`) at every entry point BEFORE
deriving anything from it, so velocities (`u = m/ρ`), internal energy, the
effective gamma and the evaluation column all describe the same in-hull
state. Skipping that produced inconsistent Riemann-face states and a
floating-point crash in HLLC on the Stage-4 abusive run. Flagged inversions
additionally swap the physical `e_int` for the table-consistent `ev.e` in
the γₑ formula (a negative post-rarefaction `e_int` would otherwise drive
γₑ through the `max()` guard to ~1e308).

Every clamp — density or inversion-target — is tallied in a per-gas counter
and reported once per primitive sweep under `verbosity >= 2` (via
`amrex::AllPrint()`, because the clamping cells usually belong to a non-IO
rank whose `Print()` output would be dropped):

```
[fluid] EOS hull clamp applied to 142 evaluations
```

`apply_prim_floor` (the pre-Riemann face guard) floors density and pressure
to the **hull edge** (`ρ_min`, smallest in-hull `p`) instead of the absolute
`effective_zero = 1e-14` used by the analytic gas models — a floored face
then has a finite, table-consistent sound speed, which is what prevents the
time-step collapse an absolute density floor allows. Constructor-time guard:
a table whose in-hull pressure minimum is not positive (a cold-curve table)
is rejected with a pointer to the deferred D10 hull-membership work.

## 9. How it is used in the code today

Two live entry points: the `TabulatedEOS` gas model (§9a) and the debug
self-test hook (§9b).

### 9a. The `TabulatedEOS` gas model (Stage 3, W6)

`Source/states/Eulerian/hydro/gas/MFP_tabulated_gas.{H,cpp}` is a `HydroGas`
backend (the same plug-in interface the ideal-gas and Eilmer models
implement) built on `EosTable`. Configured per state in Lua:

```lua
gas = {
    type   = 'tabulated',
    table  = 'data/D_fpeos.eostab',  -- path relative to the run directory
    mass   = 1.0,                    -- code units, TPG array conventions
    charge = 0.0,
    -- optional: names, ttol (default 1e-10), max_newton (default 100)
}
```

How the interface maps onto the reader:

- **`cons2prim`** (conserved → primitive, every cell every step): subtract
  kinetic energy, `invert_T_from_e`, then `eval_rt` fills pressure and
  temperature. The `Gamma` slot gets the energy-consistent effective gamma
  γₑ = 1 + p/(ρe), so the **unchanged** Riemann solvers reconstruct the
  exact face energy from p/(γₑ−1); the `SpHeat` slot carries cp from the
  general-EOS identity cp = cv + (T/ρ²)(∂p/∂T)²/(∂p/∂ρ). Floors mirror the
  ideal-gas model (`MFP_PRIM_FLOOR`).
- **`prim2cons`**: `invert_T_from_p` then `eval_rt` for e; total energy =
  ρe + kinetic.
- **`define_rho_p_T`** (initial conditions): same "positive means given"
  convention as the ideal model — (ρ,p) given → invert for T; (p,T) given →
  `invert_rho_from_p`; (ρ,T) given → direct `eval_rt`.
- **`get_speed_from_cons/prim`** (CFL time step): the **true table sound
  speed** — exact even though the flux mode is `effective_gamma`.
- **Units**: the constructor calls `nondimensionalise` with the `MFP`
  reference quantities converted from SI to the table's CGS (ρ ×10⁻³,
  u ×10², p ×10), after asserting the references are set.
- **Tracers are thermodynamically passive** (one table closes the state);
  mass/charge arrays exist so the plasma source terms work unchanged.
- **MHD is guarded**: combining an MHD state (hard-coded constant γ) with a
  tabulated hydro state aborts at configuration time with a message pointing
  at the lifting plan.

Validation: `Exec/testing/EOS-Sod-Ideal/` runs the same Sod shock tube with
the ideal-gas model and with the tabulated model on a synthetic γ=1.4 table;
the two agree to 1.5–4×10⁻⁴ per field (table-interpolation level), both sit
equally close to the exact Riemann solution, conservation is exact, and the
tabulated run costs 2.7× the ideal run's stepping time (inversions are
cheap because of the inverse-map seeds).

### 9b. The debug self-test hook (Stage 2, W5)
`EosTable::register_with_lua` (called from `MFP::read_config`,
`Source/MFP_config.cpp`) registers the Lua function

```lua
eos_table_self_test('data/ideal_synthetic.eostab', 48)
```

**only in DEBUG builds** (`#ifdef AMREX_DEBUG` — same pattern as the SDF
geometry self-tests; in a release executable the function simply does not
exist and calling it is a Lua error). It loads the named table and runs five
checks — reader integrity, e- and p-round-trips over the hull, derivative
identities, hull/out-of-range behaviour, degenerate-corner stress — printing
one `EOSTAB-SELFTEST[...] PASS/FAIL` line each, which
`Exec/testing/EOS-Table/check.py` parses. The `EOS-Table` case runs in the
standard test suite (`run_tests.py`) like any other.

To run it by hand:

```sh
cd Exec/testing/EOS-Table
sh run          # builds DIM=1 DEBUG, runs the one-zone case, checks the log
```

To build a fresh table (details in `EOS-Table/README.md`):

```sh
# synthetic ideal gas (closed-form truth, used as the test gate)
python3 ../../python_analysis/eos_table_prep.py synthetic \
    --out data/ideal_synthetic.eostab
# condition the FPEOS source into a deuterium table + QA plots
python3 ../../python_analysis/eos_table_prep.py fpeos \
    --src data/raw/FPEOS/H_EOS_09-18-20.txt \
    --out data/D_fpeos.eostab --qa qa \
    --hug-ref data/raw/hugoniot_MC2000_PRL85_1890.txt
```

## 10. How the Riemann solver uses it (`flux` modes, Stage 5 / W10)

A tabulated state can run under two Riemann-solver (the flux routine
between two cells) configurations, selected by the per-state Lua `flux`
key:

```lua
states = {
  fluid = {
    type = 'hydro',
    gas = { type = 'tabulated', table = 'data/D_fpeos.eostab', ... },
    -- no 'flux' key            -> 'HLLC_general_eos'  (the default)
    -- flux = 'HLLC_general_eos'-> same, explicitly
    -- flux = 'HLLC'            -> effective_gamma mode (Stages 3-4 path)
  },
}
```

- **`HLLC_general_eos` (default).** A dedicated solver
  (`riemann/MFP_hllc_general_eos.{H,cpp}`) asks the gas model for the face
  specific internal energy and sound speed through one combined evaluation
  per face (`HydroGas::get_face_eval_from_prim` — for a tabulated gas, a
  single rp inversion answers both). Its two-shock q-factor uses the local
  acoustic gamma Γ₁ = a²ρ/p at the face, reusing that same sound speed, so
  no extra table calls. On a gamma-law gas the base-class defaults reduce
  every one of these to exactly `MFP_hllc.cpp`'s algebra — the new solver
  agrees with `HLLC` to round-off, which is also the validation gate that
  isolates solver bugs from EOS bugs (run TPG with
  `flux = 'HLLC_general_eos'` to reproduce it).
- **`HLLC` (effective_gamma mode).** The unmodified solver reads the
  reconstructed `Gamma` slot, which the tabulated `cons2prim` fills with
  γₑ = 1 + p/(ρe). Face energy is exact; the wave-speed estimate is
  approximate (Γ₁ ≠ γₑ off the ideal corner) and the γₑ slot is linearly
  reconstructed across shocks. Kept selectable for A/B comparison and for
  wall-time-sensitive runs — it adds no per-face table calls.

Non-tabulated states are unaffected: they still require an explicit `flux`
key, and their solvers are byte-identical to pre-Stage-5 builds. All face
evaluations inside the general-EOS path go through the same hull clamping +
counter as §8 (clamp the inputs first, derive everything from the clamped
state).

Wall-time expectation: the general-EOS solver pays one combined table
inversion per face side per stage that effective_gamma mode does not; on
the 2048-cell Sod twin that measured 1.57x the effective_gamma advance
time (gate `wall-geos` in `Exec/testing/EOS-Sod-Ideal/check.py`, budget
2.0x — noise-sensitive at these short runtimes).

## 11. *(future)* Two-temperature

- **Stage 6:** two-temperature (separate ion/electron tables) extends the
  file format and this reader; the design memo is W12.
