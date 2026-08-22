# EOS-Sesame — SESAME ASCII2 extraction: tables, gates, records

The run directory for the SESAME extraction work
(`doc/eos_sesame_plan.md`; user guide: `Quickstart.md`). Holds the
emitted Track-S tables, their QA plots, and the standing one-zone
self-test harness that `run_tests.py` picks up.

## Directory contract

```
EOS-Sesame/
  README.md               this record
  Quickstart.md           user guide (selection, approximations, commands)
  data/                   committed .eostab tables (native-resolution, ~0.2-3 MB)
  qa/                     QA plots per table
  problem_definition.lua  one-zone self-test config (DEBUG-only Lua hook)
  onezone.inputs          AMReX inputs (max_step = 0)
  run / check.py          harness: DEBUG build + self-test + gate
```

The source library `sesame-unc.ascii2` (244 MB, LANL distribution,
2026-02-05 release, git 465e025f) is **not committed**; sha256
`209d5629...bad7a5ee` (full value in the `sesame-unc` entry of
`../EOS-Table/data/raw/sources.yaml`), echoed as a 12-char prefix in
every emitted table's `source:` header line.

## Generating commands (the committed tables)

One-shot regeneration of the full six-table set (five Track-S at
source-native axes + the Track-P mixture-bracket Ti table), run from
this directory:

```sh
SES=/mnt/c/Users/ktap0992/Downloads/sesame-unc/sesame-unc.ascii2   # not committed; sha256 in sources.yaml
PREP="python3 ../../python_analysis/eos_table_prep.py"

# Track S: source-native axes, auto table selection (311 -> 301)
for m in 3337 5267 7834 5251 2963; do
    $PREP sesame --src $SES --mat $m --qa qa
done

# Track P: Ti casing on the shared mixture T-bracket (H3), c_cav pinned
$PREP sesame --src $SES --mat 2963 --lT 1.25 9.0 384 \
    --lrho -5.30 1.69897 576 --c-cav 5.34 \
    --out data/ti-beta-21s_2963_trackP.eostab --qa qa
gzip -k9f data/ti-beta-21s_2963_trackP.eostab   # committed form is the .gz

sh run    # acceptance: DEBUG one-zone self-test on all six tables
```

The common-gauge mixture set (solid + interior + air on one shared
energy zero — see the eref295 record below) has its own one-batch
driver:

```sh
SES=... bash make_eref295_set.sh      # probe -> shared S -> emit + QA + spread gate
SES=... bash make_eref295_solids.sh   # Al + diamond onto the same gauge (S fixed)
```

Regeneration gate (SS1 precedent): byte-identical to the committed
tables modulo the `generator:` provenance line.

## WS5 record (2026-07-27) — first extraction set, all G1-clean

All five tables emitted at source-native resolution and pass the emitted
G1 gate (zero decreasing p intervals in ρ and T; p > 0 everywhere;
`tension_clip=0` by construction). Routes exercised per material:

| table | source | raw G1 classification | conditioning routes | band cells | hull | c_cav (km/s, source) |
|---|---|---|---|---|---|---|
| `copper_3337_s311` | 311, 117x91 | 48 flat / 43 mono — the exact-tie-line (Maxwell-constructed) trap, zero loops | 48 crossover | 1781 | 86.2 % | 3.37 (cold slope at 9.8 g/cc) |
| `deuterium_5267_s301` | 301, 208x97 | 90 mono / 7 flat, **46 decreasing-in-T intervals** (worst 0.96 rel) | 7 crossover + monoT repair | 295 (+274 monoT) | 95.0 % | 0.30 (cold slope at 0.20 g/cc) |
| `diamond_7834_s301` | 301, 112x79 | **43 loops, 1043 nonpositive tension cells** (worst drop 1.95 rel) — the full 301 pathology | 3 Maxwell-resolved + 40 deep-cold direct ramps | 1118 | 92.9 % | 8.93 (cold slope at 3.8 g/cc) |
| `hydrogen_5251_s301` | 301, 49x24 (P,U only — no Helmholtz array) | 18 mono / 6 flat | 6 crossover | 14 | 99.0 % | 0.15 (201: √(B_S/ρ₀)) |
| `ti-beta-21s_2963_s311` | 311, 130x80 | 35 flat / 45 mono | 35 crossover | 1556 | 85.2 % | **5.34** (cold slope at 5.8 g/cc — inside the Ti bulk-speed literature window without being an input) |

Notes against the plan:
- **D-c confirmed by measurement**: both 311 tables (Cu 3337, Ti 2963)
  shipped exact tie-line flats — 48 and 35 isotherms respectively —
  demonstrating 311 is not a "clean" path; and both 301 dome-carrying
  tables confirmed the loop/tension pathology (diamond) or the
  decreasing-in-T defect (D 5267) the requirements doc predicted.
- **Deep-cold Maxwell degeneracy confirmed**: only 3 of diamond's 43
  folded isotherms had a numerically resolvable equal-area tie line; the
  other 40 (vapor pressure below double-precision) took the direct-ramp
  route — the WS2b design's underflow caveat, observed.
- The c_cav fallback initially sampled the coldest-isotherm slope *at*
  ρ₀, which for 2963 landed inside the softened band and returned
  0.54 km/s (circular); fixed to the first non-band cell at/above ρ₀
  (5.34 km/s). Recorded as a pipeline trap.
- regrid LOO (leave-one-out) medians 6e-5…1.6e-3 per material; maxima
  (0.2…38) are confined to crossover-kink/hull-0 territory. Native
  resolution retained for v1; raise `--n-rho/--n-T` if in-hull LOO or
  the future Hugoniot gates demand it (plan §3.1).
- `ti-beta-21s_2963_s311` is a **Track-S pipeline-QA artifact**, not the
  production Track-P casing table (which needs the H3 bracket + R4
  extensions; plan WS2e).

## WS2e record (2026-07-28) — Track-P tables

**Design (recorded in `doc/eos_sesame_plan.md` §6c):** the Track-P Ti
table is a **SESAME-2963 single-source backbone** on the shared mixture
bracket — no solid-model seam (the WDM band that was uncontrolled in Ti
v0 is now real SESAME data; fewer seams = R3 trivially satisfied, no H5
constants). The band below SESAME's 72.5 K floor is the constant-in-T
extension `regrid_onto`'s refuse-and-fill semantics produce, hull-0.

**`data/ti-beta-21s_2963_trackP.eostab(.gz)`** — 576x384,
lρ [-5.30, 1.699] (4.9e-6…50 g/cc), lT **[1.25, 9.0]** — the H3 bracket
endpoints exact (17.782794100389228 K, 1e9 K); `--c-cav 5.34` pinned
(the native-grid Track-S value; the fine-grid fallback samples a
different cell — recorded trap). Generation:

```sh
python3 ../../python_analysis/eos_table_prep.py sesame --src $SES \
    --mat 2963 --lT 1.25 9.0 384 --lrho -5.30 1.69897 576 \
    --c-cav 5.34 --out data/ti-beta-21s_2963_trackP.eostab --qa qa
```

**New conditioning element found by the V-a audit — the thermal floor.**
The first emitted Track-P table carried 26,560 in-hull cells softer than
c_cav, of which ~4k sat below even the ideal-gas isothermal stiffness
k_B·T/m̄ (tie-line residue the band mapping missed; min c ~ 1e-15 cm/s —
in-hull dust, the exact class this work exists to kill, invisible at
Track-S resolution). Fix (now in the pipeline for ALL tables):
`dpdrho ≥ k_B·T/m̄` everywhere; sub-thermal cells join the band (hull-0 +
c_cav response); cells between thermal and c_cav *outside* the band are
genuine near-critical softening and are kept. `thermal_floored=` joins
the conditioning line; Track-S tables regenerated (1547/5690/1124/160/
1628 cells for Cu/D/C/H/Ti-alloy).

**Final V-a audit (Track-P):** p strictly monotone in ρ and T, p > 0;
zero sub-thermal cells (audit-mass caveat: 2963's Ā = 50.75, not the
2962 value); in-hull min c = 0.73 km/s (genuine near-critical, above
its thermal floor); compressed branch ρ > ρ₀ untouched (hull-1 fraction
0.919, min c 3.53 km/s); hull coverage 70.5 % (the extension band and
dome are hull-0 honestly).

**Self-test (all six tables, check.py PASS):** the Track-P table is the
best-behaved table in the stack — roundtrip-p residual 3.8e-15 with
**zero bisections** and T-recovery error 5.6e-11 (vs ~0.9 monitored on
coarse tables); fd-vs-blocks 0.27 (vs Ti v0's 49.1).

**D side:** `splice_deuterium` gained the monotonise-in-ρ pass + the G1
gate (raises on failure); `data/../EOS-Table/data/D_spliced.eostab.gz`
regenerated — G1 clean (the 5/40 rippled isotherms of v1 cleared;
monoRho_p=197 cells, max 2.3e-3 — the predicted ripple scale), **all
SPLICE-QA gates PASS with SS3-identical physics** (alignment constancy
0.0090/0.0154/0.0107, convexity 0/126,725, Hugoniot peak 4.400, PIMC
anchor median 1.19 %), C++ self-test robustness checks PASS
(fd-vs-blocks 28.1 reported, better than v1's 49.1). pytest 64 pass
(CoolProp installed into `cerberus_python`, so the old skip now runs).

**Initialization envelope (Track-P Ti, measured 2026-07-28).** Global
min sound speed anywhere on the table = 0.815 km/s (genuine
near-critical cells; band cells carry 5.34, extension ≥ thermal) — the
guard-2 degenerate-speed threshold is unreachable at ANY (ρ,T), so
EOS-driven HLLE diversion is impossible by construction. The in-hull
condensed-branch edge sits slightly ABOVE ρ₀: ρ ≥ 5.02 g/cc for
100–500 K (p ≈ 1.6–2.8 GPa there), relaxing to 4.37 g/cc by 4 kK. The
exact ambient state (4.93 g/cc, 293 K, 1 bar) is a band cell (hull-0,
p→0.69 GPa, c = 5.34 km/s): ambient pressure is fundamentally
grid-unresolvable near ρ₀ (Δρ/ρ ≈ 2.9 %/cell ⇒ ~2–3 GPa/cell at the
cold bulk modulus). Initialize either at the band cell via (ρ,T)
(finite, monotone, clamp-counted) or slightly compressed at ρ ≥ 5.02
(in-hull, ~2 GPa pre-stress — small vs multi-10-GPa drives). QA-plot
note: the pre-fix Track-P runs' QA tag collided with the s311 tables'
(same `<name>_<mat>_s<table>` base); fixed (tag now follows the output
filename) and all 18 plots regenerated and verified unique.

**Deferred / user-side:** the ℓ=2 and ℓ=3 discriminating runs and the
triple's MIXEOS gates (the air table + TiD case live outside this
repository); the Amagat-on-Ti variant additionally needs W27, which is
**not on this branch** (`mixing_rule='amagat'` still aborts). ρ-grid
extension below SESAME's 4.9e-6 g/cc floor toward the air table's 1e-8
(full R4) would need a vapor-branch model joined per R1 — deferred with
rationale, since SESAME's native floor already meets the recorded 1e-5
compromise.

## T3 record (2026-08-02) — cold-extended casing table

**`data/ti-beta-21s_2963_coldext.eostab(.gz)`** — the Track-P backbone
plus the `doc/ti_splice_plan_v2.md` cold extension (`--cold-extend`):
the sub-hull cold fill replaced by a solid model fitted to 2963's own
306 cold curve (Vinet fit over ρ ∈ [4, 12] g/cc, n=42, rms 4.7e-2:
ρ₀(0K)=4.974, B0=109.2 GPa vs raw-306 slope 110.6, B0p=3.856, Slater θ₀
543 K), capped at 0.9× the 411 melt line, blended into SESAME over
≥10-cell tanh half-bands (R3), bridged across SESAME's own gap at
expanded ρ, tension foot repaired by a second crossover pass. Extension
restricted to the fit window ρ ≤ 12 (measured validity: the H5 e-offset
is constant to ~3.5 kT inside it, drifting to 1.2e12 erg/g by 50 g/cc).
H5: one constant offset cE = −1.26e7 erg/g, std/kT = 1.006 over 690
overlap cells (gate ≤ 5). Generation: the §3b Quickstart command with
`--cold-extend`; committed form is the `.gz`.

Measured against the trackP before-state (`region_survey.py --health`):

| isochore | metric | trackP | coldext |
|---|---|---|---|
| 6.00 g/cc | sub-hull flat fraction | 90 % (28/31) | **0 %** |
| 6.00 g/cc | max adjacent dpdT log10 jump | 7.55 | **1.55** |
| 6.00 g/cc | dp/de at the T floor | 1.5e-4 | **6.8** |
| 4.93 g/cc | sub-hull flat fraction | 35 % | 14 % |
| 4.93 g/cc | dp/de at the T floor | 3.5e-6 | 4.2 |
| 4.40 g/cc | sub-hull flat fraction | 48 % | 27 % |
| 4.40 g/cc | dp/de at the T floor | 2.2e-8 | 3.8e-2 |

**Honest limits (recorded, not hidden):** (1) at ρ < ρ₀ the gap columns
remain anchor-driven — a per-column T-re-grade was tried and removed
because `monotonise_rho`'s cummax against the vapor-side fill overwrites
it (2-D monotonicity squeeze); the plan-v2 decisive-gate target
dp/de ~ 1e2 at (4.40, 30 K) is *not* reached (3.8e-2; six orders above
the before-state) and cannot be without the deferred R1 vapor-branch
model. (2) The blend's p-mismatch near ρ₀ (H5 shifts e only; max 54 %
relative at small p) plus the honest hull policy shrink the in-hull
edge at 4.40 from 3776 K to 5232 K and at 4.93 from 850 K to 2370 K —
demoted cells are altered cells. (3) The exact-flat cells that remain
below ~72 K at compressed ρ are largely physical (a Debye solid at
T ≪ θ has dp/dT → 0). The §2.2 operating-point lever (initialise at
ρ ≥ 5.02) composes with this table exactly as recorded in plan v2.

## eref295 record (2026-08-10) — common-energy-reference set + real air

**`data/{ti-beta-21s_2963_coldext_eref295, deuterium_5267_s301_eref295,
dry-air_5031_s301_eref295}.eostab(.gz)`** — the three mixture members
(solid / interior / air) on ONE shared energy gauge, per
`HANDOFF_common_energy_reference.md` + its 2026-08-10 addendum. The
defect being removed: `MixtureEOS::prepare_weights` drops trace
components without adjusting the energy target, so independently-gauged
tables inject a spurious `w_j·(e_j − <e>)` of order the *gauge spread* —
measured 32.7 code units, ~2e6× Ti's cold-range thermal span, railing
inversions at the T-axis bottom at `drop_tol = 1e-4`. (The old shifts
were dominated by each table's own `1e-3·span` positivity margin — pure
per-table arbitrariness.)

Recipe (`make_eref295_set.sh`, one batch, generator 90dd55f):
`--e-ref-state <rho_fill> 295` per material (solid 4.1856, interior
0.14775, air 2.361e-4 g/cc — e sampled off the finished surface with
the C++ reader's bilinear convention, then subtracted), then ONE shared
`--e-shift S`, `S = max deficit + 2e13 erg/g = 2.0017551224e13`
(4.85 code — the addendum's O(1)–O(10)-code target band; forced shifts
hard-fail on `min(e) ≤ 0`, no silent top-up). Post-reference deficits:
Ti −9.03e8, D2 −1.755e10, air −3.29e9 erg/g. **Achieved design-point
gauge spread: 4.29e-11 code units (target < 0.15) — PASS.** Every
table's hull minimum lands at 4.84–4.85 code; `e_ref_state:` header
records the anchor.

**The air table is the first real air EOS in the set** (replaces the
synthetic γ=1.4 table): SESAME **5031 dry air** (D. Sheppard 2018,
N₂ 0.7809 / O₂ 0.2195 / Ar 0.0096, built to fix 5030's numerical
issues; 5030 is not in this library). **Native hull only — no flank
extension** (deliberate; supersedes the handoff's lrho/lT prescription;
addendum §5 endorses no-fill-first): 576×384 over ρ 1e-7..15 g/cc,
T 100..3.4815e8 K exactly. Source is exceptionally clean: ONE
non-monotone isotherm (100 K, near the ~133 K critical point) → 1
Maxwell construction, 2 band cells; zero monotonise interventions; LOO
p rel median 3.3e-4; hull coverage **99.1%**. Ambient verification
(2.361e-4 g/cc, 295 K): in-hull, p = 0.1993 bar, dpdrho/(p/ρ) = 1.0034,
**cs = 341 m/s**. Health: zero flat isochore segments in the whole
table, dpdT > 0 everywhere.

**Molecular thermal floor (`--floor-mass-amu`)**: the D-g floor used the
201 table's per-ATOM ā, which over-floors cold molecular fluids by the
association factor (air ~2×) and would have demoted the ambient itself
to hull-0 (the shipped `deuterium_5267_s301` carries exactly this
artifact: 28% of cells). With the molecular mass (air 28.97, D₂ 4.028)
an ideal molecular gas sits ON the bound, so only cells < 0.95·kT/m are
demoted: air 1942 (0.9%), D2 12580 (5.7% at the new resolution). D2 is
also raised 208×97 → 576×384 (native bracket lρ −10..3, lT 0.602..9).

Gates this side: T1 re-baseline (flag-free re-emission at head
content-identical to the shipped coldext table); Ti gauge-mechanism gate
(p/hull/dpdrho/dpdT byte-identical, e one constant to 1.6e-11, cv/dedrho
1e-11 scale-normed); G1 clean ×3; `sh run` **check.py PASS all ten
tables** — air roundtrips at max_res 2.7e-15 with iters_max 11 / p99 9 /
1 bisection / nonconv 0 (the shared-pedestal cost on the collapsed cold
`le` range: mild), D2 9.0e-11 / iters_max 23 / nonconv 0.

Caveats (standing): (1) tables + runs regenerate TOGETHER — the gauge
change invalidates old checkpoints, plotfiles, and any stored-energy
baseline; do not mix eref295 and non-eref295 tables in one mixture set.
(2) An explicit `flux = 'HLLC'|'HLLE'|'AUSMDV'` on a tabulated state
rides the gauge-dependent Γ slot → different fluxes under any
re-reference; the mixture deck's default `HLLC_general_eos` is
gauge-clean.
(3) Air-retaining cells narrow the mixture T-bracket to
[100, 3.4815e8] K (retained-set bracketing; air-free cells keep
[17.78, 1e9]); watch a 100 K-rail counter separately from hull pins in
the 1D reproducers. (4) Under Amagat every retained component evaluates
at the mixture (p,T): expect and *attribute* air ceiling pins (the
addendum's 2.91 GPa / 667 K state wants air at ~15 g/cc) — they are the
instrument that decides whether a v2 high-ρ armor flank is ever needed,
not a failure of this work. (5) Ti's reference state sits on constructed
bridge fill (hull weight 0) — deterministic and solver-consistent, but
any future bridge change silently moves Ti's gauge; the `e_ref_state:`
header is the check. (6) `MFP_CTU_Braginskii.cpp` `check_invalid` tests
ABSOLUTE energy against a constant (the only rk4 rejection test) — out
of scope here, but it will bite re-referenced TiD Braginskii cases;
follow-up item. Acceptance still pending: MIXEOS self-tests + the
run-side 1D reproducers at production drop_tol (zero Riemann fallbacks;
pin/rail counters recorded).

## eref295 solids record (2026-08-11) — aluminum + diamond on the same gauge

**`data/{aluminum_3720_coldext_eref295,
diamond_7834_s301_eref295}.eostab(.gz)`** — two alternate solid members
emitted onto the *existing* eref295 gauge (`make_eref295_solids.sh`):
`S = 2.0017551224e13` adopted verbatim (never recomputed — recomputing
would re-gauge the shipped set; the driver cross-checks the D2 header),
with the one remaining gate that post-reference minima stay in the
O(1)–O(10)-code band: **Al 4.840, diamond 4.846 code — PASS**. Both
referenced at their 201 ρ₀ at 295 K (Al 2.70, diamond 3.515 g/cc).
Beryllium was the original request and is **not possible from this
library**: sesame-unc.ascii2 carries only Be conductivity (22021/22024)
and combined melt/shear (32020/32025, 411/431) materials — no 301/311
pressure–energy surface (the classic 2020/2021/2023 tables are not in
the unclassified release).

**Aluminum 3720** (Crockett 2003, 311, Helmholtz-consistent) gets the
full Ti-style cold extension on the shared mixture bracket
(lT 1.25..9 ×384, lρ −5.56..1.4314 ×576 = native min..10×ρ₀). This
required generalising two formerly Ti-hard-coded knobs, both opt-in and
byte-preserving for the default path (Ti probe reproduces
min_full = −9.0300113996e8 exactly): `--fit-rho` (Vinet window; the
default (4, 12) extrapolates at Al's ρ₀ = 2.70 → ambient residual
+3.99 GPa, gate FAIL; the chosen (2.2, 8.0) keeps the tension foot that
pins v₀ and gives −0.31 GPa, rms 1.06e-1, ρ₀K = 2.755, B₀ = 64 GPa vs
raw-306 slope 74) and `--z-cond` (Sommerfeld valence, Al = 3). H5
offset constancy 3.08 kT (gate 5.0); overlap p-mismatch max 7.28
(report-only; Ti ships 5.45e-1 — Al's blend edge is rougher). Reference
sits on constructed fill (hull 0), Ti caveat (5) applies identically.

**Diamond 7834** (Crockett 2006, 301) has no 411 melt table → no cold
extension possible; it gets the D2/air treatment: native hull
(ρ 3.52e-6..7.03e4 g/cc, T 72.53..1.16e9 K), raised 112×79 → 576×384.
Source is rough (43 looped isotherms, 1043 nonpositive cells → 3
Maxwell + 40 ramp constructions, 1118 band cells) but emits G1-clean;
hull 86.3%; LOO p rel median 8.0e-5. Its reference point has hull
weight 0.92 (mostly real data). Diamond-retaining mixture cells narrow
the T-bracket to [72.53 K, 1.16e9 K] — same class as air's 100 K rail.

`--c-cav` pinned for both by the Ti convention (T-floor isotherm slope
at the first column above the crossover tension foot, ~1.1–1.2×ρ₀):
Al **5.82** km/s at 3.04 g/cc, diamond **12.7** km/s at 3.83 g/cc. The
data-derived default measures ON the repaired foot and lands 2.7 (Al) /
1.06 (diamond) km/s — far below the bulk sound speeds (~5.3 / ~11.2);
pinning keeps the cavitated-response floor physical and
grid-independent.

Caveats (1)–(6) of the eref295 record apply unchanged; the five-table
family {Ti, D2, air, Al, diamond} shares one gauge, so any mixture
subset is legal — but never together with old-gauge tables.

## Quiet-foot re-emission (2026-08-11) — solid ICs at ~ambient pressure

All three solid tables re-emitted (`doc/eos_quiet_foot_plan.md`; D2/air
byte-identical, S and the gauge unchanged to the last digit, five-table
spread 6.86e-11 code). Motivation: a solid initial condition
(fluid-represented) must not over-expand — but the emitted 295 K
isochores at solid density read 1.45 (Ti), 96.7 (Al!), 0.58 (diamond)
GPa where raw data and physics say ~bar scale.

Two mechanism classes fixed:
- **Ramp-level artifact (all solids)**: the raw 311s pin the whole
  sub-solid region — vacuum THROUGH ρ₀ — at one placeholder value
  (Al 7.9e-82, Ti 2.0e-254 GPa); H1 replaces the flat with a ramp
  whose log-linear span parks anchor-scale pressure at ρ₀. Fix:
  `--p-foot` floor-hugging shape (see flag table in Quickstart) —
  opt-in, legacy byte-identical when absent.
- **Solid-model spike below the spinodal (Al, latent for any coldext
  material)**: a single-cell +96.3 GPa value at 1.943 g/cc — the
  Slater theta's hard B-floor kink at the Vinet spinodal — propagated
  table-wide by the crossover envelope (NOT the Maxwell construction;
  ramp-only reproduced the identical plateau). Fixes, all
  unconditional in `cold_extend_stage`: spinodal validity floor
  (`rho_min=` in the cond string), construction right-envelope clamp
  (hull-0 cells capped below the smallest genuine p at larger ρ on
  their isotherm), ramp-only crossover-2, and a hull guard that
  refuses to emit if pristine cells move > 50% (measured ladder:
  legitimate seam/dome ripple 5–15%, the defect 3400%).

**ρ₀-aligned lattices (same-day follow-up)**: the first emission left
GPa-scale readings *at* ρ₀ although every stored node ≤ ρ₀'s cell was
quiet — the lattice had no node at ρ₀, so the reader's linear-in-value
bilinear blend across the single knee cell (1 bar → ~1 GPa genuine
compression at the first node above ρ₀K) dominated any sample inside
that gap. Fix: the Ti/Al `--lrho` upper endpoints are now chosen so
log10(ρ₀) is EXACTLY a lattice node (Ti hi = 1.7038353223, node 492;
Al hi = 1.4278989136, node 493 — see driver comments). Off-node sample
densities inside the knee cell still read the blend — that smear is
irreducible at 2.84%/cell and is precisely what the piecewise-uniform
axis proposal would sharpen.

Results at 295 K (`foot_achieved=` in each header): Ti 0.198 bar at
its 4.1856 fill and 1 bar at ρ₀ = 4.93 exactly; Al **1 bar at ρ₀ =
2.70 exactly** (was 96.7 GPa); diamond 5798 bar unchanged-and-correct
(its ρ₀ node is genuine +0.32 GPa source data; fill at ≤3.50 g/cc for
a quiet IC). All sub-ρ₀ cold isochores hug the foot to ~300 K then
track the raw thermal rise; above ρ₀K the compression branch matches
raw. Gates: pytest 82, T1 (D2/air), check.py PASS ×12, minima band +
spread (1.22e-10 code) PASS. Deck fill densities should still be read
off each table's own 295 K isotherm; ρ₀ itself is now a safe choice
for Ti and Al.

## Harness

`sh run` — builds the DIM=1 DEBUG executable, runs the one-zone config
(`max_step = 0`), and gates the `EOSTAB-SELFTEST` lines via `check.py`.
Tier-2 policy for all five tables (real conditioned data): reader,
roundtrip-e/p, identities, hull, corner gated; `fd-vs-blocks` reported
not gated — the cavitated-response floor (plan §3.0) makes the
derivative blocks deliberately stiffer than the value surface inside the
crossover band.

Self-test record (2026-07-27, first harness run): **check.py PASS, all
five tables** — reader / roundtrip-e / roundtrip-p / identities / hull /
corner all PASS per table; roundtrip residuals ≤ 1.0e-10 with
`nonconv=0` everywhere (iters_max ≤ 28, bisection fallbacks ≤ 67);
identities ≤ 5.1e-16. `fd-vs-blocks` reported 1.38–16.3 per table (the
cavitated-response band, by construction). `max_Terr ≈ 0.9` monitored —
the known flat-in-T / residual-only-gates class. Noted: the 311-derived
tables carry astronomically small deep-cold vapor pressures from the
source data itself (Cu min p 5.7e-235, Ti 2.0e-244 erg/cc), stretching
the inverse-map p-axis over ~250 decades; the guarded Newton absorbed it
(seeded, ≤ 22 bisections). If seed quality ever degrades on such tables,
a `--p-trim` floor for the inverse-map axis is the recorded lever.
