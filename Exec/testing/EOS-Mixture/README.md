# EOS-Mixture — Dalton mixture validation (Stage 9, W28)

Validation for the `tabulated_mixture` gas model
(`Source/states/Eulerian/hydro/gas/MFP_mixture_gas.{H,cpp}`) under the Dalton
(partial-pressure) mixing rule. Plan: `doc/eos_mixture_dalton_plan.md` §11;
design background: `doc/eos_wide_range_and_mixtures.md` §5.

## What runs

Five Sod shock tubes on identical numerics (1024 cells, minmod, RK2,
`HLLC_general_eos` via the no-`flux`-key config default), driven by `sh run`:

| run | gas | alpha | isolates |
|---|---|---|---|
| `single` | `tabulated`, one γ=1.4 table | — | baseline |
| `mix_pure1` | mixture, two IDENTICAL γ=1.4 tables | 1 everywhere | pure-cell short-circuit (component 0) |
| `mix_pure0` | same | 0 everywhere | pure path through the DERIVED last fraction |
| `mix_uniform` | same | 0.3 everywhere | the full Dalton solve in every cell |
| `mix_binary` | γ=1.4 + γ=5/3 tables | step at x=0 | the real two-material tube (M3/M4/M5) |

The synthetic tables are closed-form ideal-gas tables
(`eos_table_prep.py synthetic`), for which Dalton is analytically exact —
that is what makes the gates below sharp.

## Gates and measured values (first passing run, 2026-07-12)

- **M1 pure round-off** (`pure{1,0}-roundoff-*`): the pure-cell path is
  required to be formula-identical to `TabulatedEOS`; measured
  Linf/range = **0.0 exactly** (bitwise) on rho/p/u/T/nrg for both runs.
- **M1b uniform mix** (`uniform-twin-*`): alpha=0.3 vs single measures pure
  bilinear-interpolation error at the partial densities (Dalton exact
  analytically). Measured L1/range ≤ **2.5e-4** (tol 1e-3).
- **M2 self-test** (`selftest-*`, printed by the ctor when the gas def
  carries `self_test = <n>`): driver round-trips over an alpha-grid ×
  hull-sweep (max residual **8.6e-11**, ttol=1e-10, 0 non-converged),
  pure-parity **bitwise**, frozen sound speed vs an isentropic FD of the
  mixture pressure (max cs² error **1.2e-2**, tol 5e-2), and the drop_tol
  crossing (p jump **3.1e-12** — kink-free thanks to the dilute ideal-limit
  scaling of below-hull partial pressures).
- **M3 binary tube**: stable; alpha ∈ [0, 1] exactly; **45** mixed contact
  cells (bounds [1, 60]); positivity.
- **M4 conservation**: total mass/energy and BOTH per-component masses
  (Σ rho·alpha tracer) drift ≤ **3.3e-16**.
- **M5 wall-time**: binary/single advance-time ratio **2.39** (budget 2.5
  for N=2; noise-sensitive at ~2 s runtimes — same caveat as the
  EOS-Sod-Ideal wall gates).
- **W24 default flux**: the log must show the config default selecting
  `HLLC_general_eos` for the mixture gas (the `needs_general_eos_solver()`
  predicate).

## Notes

- The pure round-off gates are deliberately posed on IDENTICAL ideal-gas
  tables: on a general nonlinear table Dalton does not reproduce the
  single-material run (shared-volume rule), so that exact-reproduction gate
  belongs to Amagat (plan W27/W28a).
- Dilute condensed components under Dalton land below their density hull by
  construction; the model handles that with the ideal-gas low-density limit
  (p and dpdT scaled by rho_k/rho_hull_min at the hull-edge row) plus
  per-component clamp accounting. See the header of `MFP_mixture_gas.H`.
