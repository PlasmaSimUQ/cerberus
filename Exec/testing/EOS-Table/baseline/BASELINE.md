# W1 baseline record — 2026-07-07

- **git**: ffc2135 (`feature_general_eos`), working tree carried unrelated
  local modifications to `Double-Rarefaction/{problem_definition.lua,run}`
  (recorded in `run1_raw.log` header) — the baseline is of the tree *as run*.
- **Build config**: per-case self-builds in `Exec/local/` — Couette: DIM=2
  USE_EB=TRUE AMREX_PARTICLES=FALSE; Double-Rarefaction: DIM=3 USE_EB=FALSE
  AMREX_PARTICLES=FALSE USE_PRIM_FLOOR=TRUE (the make flag defining
  `MFP_PRIM_FLOOR`; **default TRUE**); Viscous-Vortex: DIM=2. gcc 11.4.0,
  OpenMPI (see log header), AMReX 22.06, 20-core host, cases run 8–10 ranks.
- **Suite content**: exactly 3 cases have `check.py` and are executed by
  `run_tests.py`: Couette, Double-Rarefaction, Viscous-Vortex.
  **No single-rank case exists in the suite** (review amendment #2's
  "deliberately deterministic single-rank case" is unavailable as-is).

## Verdicts (both runs)

| case | run 1 | run 2 |
|---|---|---|
| Couette | PASS | PASS |
| Double-Rarefaction | PASS | PASS |
| Viscous-Vortex | PASS | PASS |

`run_tests.py`: "Completed with 0 failed tests" in both runs.

## Determinism verdict

Final-plotfile sha256 checksums (`run{1,2}_plotfile_checksums.txt`)
**differ between the two runs for all three cases** → the suite is NOT
bit-reproducible run-to-run under multi-rank MPI (reduction/atomic order).
Consequence, as anticipated by the plan amendment: **all later regression
gates are verdict-identity gates**; plotfile bit-comparison is only
meaningful if a single-rank case is added to the suite (candidate: the
Stage-2 one-zone EOS-Table case, which runs `mpirun -n 1` — but it has no
plotfiles; a bit-compare gate would need a dedicated 1-rank hydro case).

## Files

- `run1_raw.log`, `run2_raw.log` — full harness output (incl. builds)
- `run{1,2}_<case>_log.txt` — per-case run logs
- `run{1,2}_plotfile_checksums.txt` — final-plotfile checksums
