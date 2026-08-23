# air_lowT_extension — eref295 member 5031 re-emitted with its T floor at ~17.7 K

Plan: `doc/eos_air_lowT_extension_plan.md` (design doc, not committed).
Built 2026-08-23. Decisions: D1-a native-spacing prepend, D2-a uniform
ideal-thermal branch, D3-a new rows hull 0 (constructed), D4-a `--T-floor`
stage, D5-a NEW file (the shipped `../data/dry-air_5031_s301_eref295` is
kept), D6-a seam kink accepted, D8-a CoolProp check as a script here.

## What this is

A mixture is bracketed in T by the HIGHEST native floor among its
retained components; in the eref295 set that floor belonged to SESAME
5031 (100 K) while every other member reaches lT 1.25 (17.78 K).
`data/dry-air_5031_s301_eref295_Tf1p25.eostab(.gz)` is the SAME emission
as the shipped eref295 table of that member (same arguments, same shared
shift S = 2.0017551224e13 erg/g, same 295 K reference state) plus 44 rows
prepended BELOW the native 100 K floor on the native log spacing
(h = 0.0170803): T floor 100 K → 17.72 K (lT 1.24847; target 1.25), 428 T
nodes instead of 384. Mixture cells retaining this member now bracket on
[17.78 K (set floor), 3.4815e8 K (this member's ceiling)] instead of
[100 K, 3.4815e8 K].

Construction (stage 7c, `eos_tools/lowT.py`), per density column with
anchor row (T0 = 100 K, p0, e0), m = 28.97 amu, R_s = k/m:

    e = e0 - 5/2 R_s (T0 - T)                 gas cv, rotations classical
    p = max(p0 - rho R_s (T0 - T), p0 T/T0)   anchored slope / scaled row

dpdT is continuous and > 0 in the extension; both branches are <= p0 so
the T-enforcement can never lift a native row; rho-monotone wherever the
anchor row is; cv constant → linear T(rho,e) inversion. Measured on 5031:
`T_ext_vapor=0 anchor=4136 scaled=21208` — the 100 K row is 1–3 % SUB-ideal
in the gas phase (real-gas second virial), so the gas columns take the scaled
branch, which preserves the source's own compressibility factor at the
anchor (dpdT = p0/T0 ≈ 0.97 rho R_s); 94 condensed columns take the
anchored branch.

## Gates (make_air_lowT.sh, 2026-08-23 — all PASS)

- G-a native byte identity vs the shipped table: lrho identical; p, e,
  hull, dpdrho, dedrho identical on ALL native rows; dpdT and cv identical
  on all native rows except the 100 K seam row, where the PCHIP stencil
  legitimately changes from end-point to interior (report: dpdT max rel
  0.98 / median 1.2e-2; cv max 0.86 / median 3.3e-2 — liquid columns with
  large native dpdT). lT native nodes reproduce to 2.3e-12 (the header
  stores 12-digit endpoints; nodes are reconstructed). TRAP found and
  fixed on the way: `monotonise_T`'s strictifying ramp (1e-12·j·|F|) is
  indexed from the first row, so prepending rows shifted every native node
  by ~4e-11 relative; the CLI now enforces the native block on its own
  index origin (flag-free emission stays byte-identical to the shipped
  table — T1 re-verified 2026-08-23: blocks identical, generator line only).
- G-b gauge: `e_shift` header identical; e(rho_ref = 2.361e-4 g/cc, 295 K) =
  2.0017551224e13 unchanged; min(e) > 0 under the forced S (2.0014e13);
  6-member design-point spread 1.216081e-10 code (unchanged).
- G-c conditioning: G1 428/428 isotherms monotone, 0 loops, nonpos 0;
  monoT/monoRho moves 0; c_cav 1.6941e4 cm/s UNCHANGED (the data-derived
  default now measures the coldest NATIVE isotherm — another trap closed);
  thermal floor 41296 floored (22112 native + 19184 new: the scaled-branch
  gas cells sit 3 % under kT/m and are floored to it, in-hull semantics
  moot since hull 0), 7442 demoted (1942 native + 5500 new dome cells,
  sub-thermal → banded → cavitated-response floor dpdrho >= c_cav² =
  kT0/m on those 5500 cells, i.e. up to 5.6× the ideal stiffness at
  17.7 K — block-level, conservative, same policy as the native band);
  hull coverage 88.93 % (native 99.1 % × 384/428 — the new rows are all
  hull 0 by decision).
- G-e physics (REPORT-ONLY, `qa_coolprop_check.py`, Lemmon et al. 2000 via
  CoolProp fluid `Air`, 60–100 K, rho <= 0.5 g/cc; CSV in qa/):
  vapor side (rho < 1e-3): p median 0.8 %, max 2 % (the 0.8 % is 5031's
  own offset from Lemmon at 100 K, inherited); energy drop from 100 K
  median 0.2–0.3 % (gas cv matches) — except at 60 K where cells between
  ~6e-4 and 1e-3 g/cc are two-phase in reality (p max 2.8×, de max 84 %).
  Dense side (1e-3–0.5 g/cc): construction — Lemmon's two-phase pressure
  is 4× (80 K) to 84× (60 K) lower than the scaled-row branch (median),
  energy drop 70–90 % larger (latent heat not modelled); at 100 K itself
  the anchor row already differs from Lemmon by 10 % median (SESAME vs
  Lemmon dome placement).
- G-f solver (`sh run` from `..`, 2026-08-23, DEBUG DIM=1, all 13 tables):
  **check.py PASS**; new file: reader dims 576x428 hull_frac 0.889,
  roundtrip-e n=2048 max_res 3.4e-15 iters_max 11 / p99 9, 0 bisections,
  nonconv 0; roundtrip-p iters_max 6 / p99 4; corner PASS (max_Terr
  5e-12); identities 4.4e-16; fd-vs-blocks 0.115 (report). MIXEOS
  self-tests and the 1D reproducers with the rail counter moved to
  17.7 K are user-side.

## Limitations (read before using)

1. Below 59.75 K at every density, and below 100 K above ~1e-3 g/cc, the
   rows are CONSTRUCTION, not data — no validated EOS for this material
   exists there.
   They are hull 0; the solver uses the hull flag only for its in-hull
   pressure-minimum floor (unchanged by this file).
2. Metastable ideal vapor: no condensation, dome, or latent heat below
   100 K. An adiabatically expanding ambient parcel (rho ~ T^2.5) never
   reaches the dome, so this is exact for it; compressed-then-cooled gas
   is not represented.
3. Condensed density: gas cv (real ~2×), Grüneisen ~0.4 (real ~2), no
   solid phases; pressure is cold-curve dominated so p is within ~35 % at
   1 g/cc and < 5 % above ~1.5 g/cc.
4. The 100 K row (incl. its Maxwell flat and demoted cells) is inherited
   down every column; the seam carries a dpdT kink at condensed density.
5. The rail MOVES (100 K → 17.7 K); count it there, and watch the hull-0
   fraction of cells retaining this member as the honesty metric.

## Files

- `make_air_lowT.sh` — fixed-S driver, gates inline, gzip; run from here.
- `qa_coolprop_check.py` — G-e report; writes `qa/coolprop_check_*.csv`.
- `data/`, `qa/`, `*.log` — gitignored (SESAME-derived / generated).

sha256 (2026-08-23 emission):
- `dry-air_5031_s301_eref295_Tf1p25.eostab.gz` a41bf08c226299ae1ba131352e3592c634c4cd8629f6c889c1782b23e854ebf6
- `dry-air_5031_s301_eref295_Tf1p25.eostab`    c97793d5554d086fe8cf30f2ef2993d87ce6004ac2ef62c2486e54151d8aab1d
