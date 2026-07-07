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

## Determinism verdict (corrected 2026-07-08)

Final-plotfile sha256 checksums (`run{1,2}_plotfile_checksums.txt`) differ
between the two runs for all three cases — but a follow-up probe
(Double-Rarefaction run twice, plotfiles compared with AMReX `fcompare`,
built in `amrex/Tools/Plotfile` via `make`) shows **every physical field is
bit-identical (absolute error 0 on all 17 fields)**. The raw-byte checksum
differences come from exactly two benign sources:
1. the `cost` field — a load-balancing diagnostic storing wall-clock time
   per cell, i.e. timing noise by construction;
2. rank-to-file packing (`Cell_H` / `Cell_D_*` layout), which permutes
   identical numbers between files.

Consequence for later gates: the bulk regression gate is still
**verdict-identity** (cheap, covers everything), but **field-level
bit-identity via `fcompare` is available and should be used as the
spot-check** — gate = all fields except `cost` show zero error. Whole-file
checksums are NOT a usable gate; do not resurrect them.

Caveat: verified on Double-Rarefaction (3D, 10 ranks); assumed for the
other cases until an fcompare spot-check says otherwise.

## Files

- `run1_raw.log`, `run2_raw.log` — full harness output (incl. builds)
- `run{1,2}_<case>_log.txt` — per-case run logs
- `run{1,2}_plotfile_checksums.txt` — final-plotfile checksums
