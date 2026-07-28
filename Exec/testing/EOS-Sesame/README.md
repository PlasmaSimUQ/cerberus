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
