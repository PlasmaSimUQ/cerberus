# Quickstart — SESAME → Cerberus EOS tables (`sesame` subcommand)

How to turn a material from a LANL SESAME ASCII2 library file into a
Cerberus `.eostab` table, what the tool does to the original data on the
way, and how to check the result. Design record: `doc/eos_sesame_plan.md`;
binding conditioning requirements:
`doc/eos_table_conditioning_requirements.md`.

All commands below are run from this directory
(`Exec/testing/EOS-Sesame/`) in an environment with numpy/scipy/matplotlib
(e.g. `conda activate cerberus_python`).

```sh
SES=/path/to/sesame-unc.ascii2     # the LANL ASCII2 library file
PREP="python3 ../../python_analysis/eos_table_prep.py"
```

The library file itself is **not** committed (244 MB, LANL distribution);
its identity is pinned by sha256 in every emitted table's `source:` line.

---

## 1. Selecting a material

List everything in the file, or filter by name / ID:

```sh
$PREP sesame --src $SES --list                  # all 752 materials
$PREP sesame --src $SES --list --grep copper
$PREP sesame --src $SES --list --grep titanium
```

Each row shows the material ID, the name from its 101 comment table, and
the SESAME tables it carries. **Only materials with a 301 or 311 table can
produce an EOS closure.** Read the table IDs like this:

| table | contents | usable? |
|---|---|---|
| 301 | total EOS: P, U(, A) on a (ρ,T) grid | yes — the normal source |
| 311 | total EOS, Maxwell-constructed (equilibrium tie lines built in) | yes — preferred when present |
| 303/304/305 | ion+cold / electron / ion component EOS | accepted with a warning: a *partial* EOS, physically incomplete as a single-table closure (future two-temperature feed) |
| 306 | T=0 cold curve (1-D) | not a closure by itself; consumed by `--cold-extend` (Vinet fit input, approximation 12) |
| 401/411/412 | vapor dome / melt boundaries | 411 consumed by `--cold-extend` (melt cap); 401/412 not used |
| 5xx/6xx | opacity / conductivity | not EOS data |

ID families: 3xxx-9xxx EOS proper; 1xxxx opacity; 2xxxx conductivity;
3xxxx melt/shear — the same element appears in several families, so e.g.
"copper" lists 3332…3337 (EOS) *and* 23331… (conductivity). You want the
EOS-family entry.

Known contents of the unclassified 2026 release (see
`doc/eos_sesame_plan.md` §1.1): Cu 3336/3337, diamond (HDC) 7830/7834,
D 5263/5266/5267, H 5250/5251, Ti **alloys only** (2961/2962,
2963 = Ti-Beta-21S; there is **no pure-Ti EOS** in the file), **no
tritium** (mass-scale a deuterium table if needed).

## 2. Producing a table

```sh
# one material, defaults (table auto-select: 311 if present, else 301):
$PREP sesame --src $SES --mat 3337 --qa qa

# the current four-material set used by this case:
for m in 3337 5267 7834 2963; do
    $PREP sesame --src $SES --mat $m --qa qa
done
```

Output lands in `data/<name>_<mat>_s<table>.eostab` (e.g.
`copper_3337_s311.eostab`); `--out` overrides. `--qa qa/` writes the
standard QA plots (isotherms, cv heatmap, hull map). Runtime ≈ 20 s per
material (single streaming pass of the library file).

Options:

| flag | meaning (default) |
|---|---|
| `--table N` | force a SESAME table ID (auto: 311 → 301) |
| `--n-rho/--n-T` | target grid counts (source counts — the floor, not the recommendation; raise if the LOO/QA numbers say so, see §4) |
| `--T-min K`, `--rho-min g/cc` | trim the source grid (T=0 and ρ=0 rows are always dropped — log axes cannot hold zero) |
| `--c-cav km/s` | override the cavitated-response sound speed (default: √(B_S/ρ₀) from the 201 table, else the coldest-isotherm slope at the first non-band cell at/above ρ₀) |
| `--fit-rho LO HI` | `--cold-extend` Vinet fit window in g/cc (default 4 12 — 2963's measured window; scale to ~0.8–2.4× the material's own ρ₀, keeping the 306 tension foot inside or v₀ is unpinned) |
| `--z-cond N` | `--cold-extend` conduction-electron count for the Sommerfeld term (default 4 = Ti; Al 3) |
| `--T-floor LOG10T` | extend BELOW the native T floor down to 10^LOG10T K by prepending rows on the native log spacing (native nodes untouched): anchored ideal-thermal branch off the floor row — e = e₀ − 5/2 (k/m)(T₀−T), p = max(p₀ − ρ(k/m)(T₀−T), p₀·T/T₀) — hull 0; needs `--floor-mass-amu` (molecular). A mixture is bracketed by the highest native floor among its retained components; this lowers one member's floor (first use: the eref295 member with a 100 K native floor → 17.7 K, `air_lowT_extension/`); see `doc/eos_air_lowT_extension_plan.md` for what it is NOT (metastable ideal vapor, gas cv at condensed density, no solid phases) |
| `--p-foot [BAR]` | quiet-foot target (1.0 if given bare; absent = legacy): crossover bridges hug this pressure across clipped/tension spans instead of the full-span log ramp that parks anchor-scale GPa at ρ₀; steep rise confined to the final cell before the data anchor; 1e-3 T-tilt keeps dpdT > 0. The achieved value at the `--e-ref-state` density is recorded as `foot_achieved=` (genuine data is never overridden — diamond's ρ₀ node carries real +0.32 GPa and reports honestly) |
| `--out path` | output path |

## 3. What the tool does to the original data — every approximation

Read this section before trusting a table. The one-line philosophy
(requirements doc): *steep monotone everywhere, never exactly flat, never
folded* — solvers break on flats and folds, so equilibrium coexistence
structure is deliberately traded away. In pipeline order:

1. **Zero rows dropped.** SESAME grids start at ρ=0 and T=0; the `.eostab`
   format is log-uniform in both axes, so those rows cannot be
   represented. The T=0 isotherm (cold curve) is *not* in the emitted
   table — the coldest surviving isotherm is the table's floor (printed
   at ingest).
2. **Units.** SESAME GPa → erg/cc (×1e10), MJ/kg → erg/g (×1e10),
   Mg/m³ ≡ g/cc. Stored dimensional CGS; Cerberus nondimensionalises at
   load.
3. **The Helmholtz array is ignored** (where present). Only P and U are
   used; thermodynamic consistency of derived quantities is *not*
   re-imposed from a free energy.
4. **Two-phase regions are replaced** (the load-bearing approximation).
   Every isotherm is classified (G1 scan): monotone / exact tie-line flat
   (311-style Maxwell data) / folded (301-style van der Waals loops,
   possibly with negative-pressure tension cells).
   - Loops: an equal-area **Maxwell construction** locates the
     equilibrium tie line where it is numerically resolvable; at
     deep-cold temperatures the vapor pressure underflows double
     precision and the band is bridged directly instead.
   - Every tie line/flat/tension band is then replaced by a **strictly
     monotone geometric ramp** in ρ between the branch endpoints. p is
     never pinned flat and never ≤ 0 anywhere in the emitted table
     (`tension_clip=0` always).
   - Energy across a replaced band follows the **two-phase lever rule**
     (linear in specific volume between the band edges) — latent-heat
     content is preserved.
   - Consequences: exact coexistence (flat-P phase equilibrium, correct
     dome shape, tension/metastable states) is **not representable**;
     shock *endpoints* and energy budgets are preserved; split-shock /
     slow near-dome expansion physics is approximate by construction.
     This is the documented UC1-d trade (overdriven-shock scope, no
     strength model in Cerberus).
5. **Monotone-in-T enforcement** on p and e at every density (the
   single-branch inversion contract). SESAME data violating it (e.g. 46
   decreasing-in-T intervals in D 5267) is raised to the running maximum.
6. **Hull honesty.** Every cell whose value was *replaced* (crossover
   band, tension cells, significantly-monotonised cells) is flagged
   `hull=0`. At runtime such cells evaluate finite, table-consistent
   values but are counted and reported as clamped ("EOS hull clamp"
   lines at verbosity ≥ 2). Typical hull coverage: 85–99 % depending on
   how much dome the material's grid covers.
7. **Resampling to log-uniform axes** (the format requires it; axes are
   reconstructed from range+count). SESAME's hand-clustered grid points
   (extra resolution near melt/ambient/ionization) are lost; the
   leave-one-out (LOO) regrid error is printed — median is typically
   1e-4–1e-3, the max sits at crossover kinks (hull-0 territory). If
   in-hull LOO is too large, raise `--n-rho/--n-T`.
8. **e-shift.** First-principles/chemical-model energies are negative at
   low T; a constant (recorded as `e_shift:` in the header) makes e > 0.
   Physically inert if every quantity comes from the same table — but a
   `tabulated_mixture` compares energies ACROSS tables (its component-drop
   renormalisation), so a mixture set must share ONE gauge: see §3c.
9. **cv floor.** cv = ∂e/∂T is floored at 1e-3 · (3/2 k_B/m̄) (m̄ from the
   201 table); floored-cell count recorded.
9b. **Thermal stiffness floor.** The emitted `dpdrho` block is floored
    at k_B·T/m̄ everywhere — no single-phase fluid is isothermally softer
    than an ideal gas, so cells below it are tie-line residue the band
    mapping missed; they are floored, joined to the band, and hull-0
    (`thermal_floored=` count in the provenance line). Genuine
    near-critical softening above the thermal floor is left untouched.
    **Molecular fluids need `--floor-mass-amu`** (air 28.97, D₂ 4.028):
    the 201 ā is per *atom*, so k_B·T/m_atom over-floors the cold
    molecular region by the association factor (air ~2×) and would
    demote the ambient itself to hull-0. With the flag, an ideal
    molecular gas sits exactly ON the bound, so only cells below
    0.95·k_B·T/m are demoted; marginal cells are floored but stay
    in-hull (`floor_mass_amu= floor_floored= floor_demote_tol=` in the
    provenance line).
10. **Cavitated-response stiffness (decided D-e).** On band cells the
    emitted `dpdrho` derivative block is floored at c_cav² — the
    material's own bulk sound speed (√(B_S/ρ₀), e.g. Cu 3.37, diamond
    8.93, Ti-Beta-21S 5.34 km/s). The reported sound speed in the dome
    never falls toward zero (a cavitated solid still carries elastic
    waves), which is what keeps Riemann faces off the degenerate-speed
    guard. Cost: the derivative blocks deliberately disagree with the
    finite-difference of the p surface *inside the band* (confined,
    reported-not-gated). `dpdT` and `cv` are untouched, so inversion
    behaviour is unchanged.
11. **Alloy caveat.** 2963 is Ti-Beta-21S (an alloy): Z̄ = 23.17602,
    Ā = 50.74763, ρ₀ = 4.93 from its 201 table are the *alloy* constants
    and are what the emitted `composition:` line carries. Do not expect
    pure-Ti transition pressures from it. (An earlier revision of this
    note quoted Z̄≈21.1/Ā≈45.9 — those belong to material **2962**, a
    different unnamed Ti alloy; see `doc/ti_splice_plan_v2.md` B1.)
12. **Cold extension (opt-in, `--cold-extend`).** Replaces the sub-hull
    cold fill with a solid model fitted to the material's own 306 cold
    curve (Vinet + Slater–Debye ions + Sommerfeld electrons), capped at
    0.9× the 411 melt line (never evaluated in the liquid), blended into
    SESAME over ≥10-cell tanh half-bands where they overlap, and bridged
    (linear in log T, hull-0, declared) across SESAME's own gap where
    they don't. Restricted to the Vinet fit window (ρ ≤ 12 g/cc for
    2963 — measured validity, see `doc/ti_splice_plan_v2.md` T3); one
    constant e-offset aligns the energy zeros (H5, gated on constancy);
    p is never shifted (overlap mismatch reported). A second crossover
    pass repairs the solid model's tension foot (ρ < ρ₀) with the
    standard monotone ramp + lever-rule e. All model-supplied and bridge
    cells are hull-0. The fit window and valence are per-material knobs
    (`--fit-rho`, `--z-cond`): the defaults are 2963's, and applying them
    to a material with a different ρ₀ makes the ambient-pressure gate
    fail on extrapolation (measured on Al 3720: +3.99 GPa at ρ₀ under
    the default window; −0.31 GPa with `--fit-rho 2.2 8.0`).

13. **Sub-floor T extension (opt-in, `--T-floor LOG10T`).** Prepends rows
    BELOW the native T floor on the native log spacing (native nodes are
    never resampled): per density column, anchored on the floor row
    (T₀, p₀, e₀), e = e₀ − 5/2 (k/m)(T₀−T) and p = max(p₀ − ρ(k/m)(T₀−T),
    p₀·T/T₀), m the `--floor-mass-amu` molecular mass, hull 0. Exact for
    an ideal vapor, cold-curve-dominated at condensed density, merely
    safe (positive, strictly monotone, continuous dpdT, constant cv)
    between — there is NO condensation, latent heat, dome or solid phase
    below the anchor row. The floor term is p₀·T/T₀ rather than ρkT/m
    deliberately: in two-phase columns an ideal-gas floor would exceed p₀
    and the T-enforcement would lift native rows. First used on the
    eref295 member whose 100 K native floor set the mixture bracket
    (100 K → 17.7 K, `air_lowT_extension/`); records
    `T_floor= T_ext_rows= T_ext_lT= T_ext_model=ideal T_ext_vapor=
    T_ext_anchor= T_ext_scaled=` in the header.

Every count above is recorded in the emitted header's `conditioning:`
line, e.g.:

```
conditioning: cv_floor=... cv_floored=424 monotonised=355
  maxwell=constructed:0/flats:48/ramp:0 band_cells=1781
  c_cav=3.3741e+05 resp_floored=1471 tension_clip=0
  monoT_p=808(6.5e+02) monoT_e=0(9.0e-11) monoRho_p=0(1.2e-10)
```

(`flats` = shipped 311 tie lines crossed over; `ramp` = deep-cold
degenerate isotherms bridged directly; `tension_clip` is **always 0** —
if you ever see otherwise, the table predates this pipeline.)

## 3b. Track-P mode: prescribed axes (mixture-set tables)

Tables destined for a `tabulated_mixture` state must share the exact
T-bracket (the constructor intersects component T-hulls bitwise). Pass
the axes explicitly:

```sh
# the Ti solid table on the shared mixture bracket [17.78 K, 1e9 K]:
$PREP sesame --src $SES --mat 2963 --lT 1.25 9.0 384 \
    --lrho -5.30 1.69897 576 --c-cav 5.34 \
    --out data/ti-beta-21s_2963_trackP.eostab --qa qa
```

Cells outside the SESAME source span (e.g. below its 72.5 K T-floor)
are nearest-filled (constant-in-T) and marked hull-0 — an *extension*,
not data. Pin `--c-cav` for reproducibility across grids (the fallback
samples a grid-dependent cell; record the value you pin and why).

```sh
# the cold-extended solid table (approximation 12; plan v2 T3-T5):
$PREP sesame --src $SES --mat 2963 --lT 1.25 9.0 384 \
    --lrho -5.30 1.69897 576 --c-cav 5.34 --cold-extend \
    --out data/ti-beta-21s_2963_coldext.eostab --qa qa
```

## 3c. Common-energy-reference sets (`--e-ref-state` / `--e-shift`)

Each table's `e_shift` is an independent gauge constant, and the
`tabulated_mixture` component-drop renormalisation takes weighted
*differences of energies across tables* — so a mixture set built from
independently-gauged tables carries an O(10 code units) spurious energy
term at any finite `drop_tol` (measured: inversions rail at the T-axis
bottom). The fix is one shared gauge, applied at generation time (the
header `e_shift` is documentary; the C++ reader never applies it):

1. `--e-ref-state RHO T` — subtract `e(RHO, T)` (each material at its own
   fill density, one common T, sampled off the finished surface with the
   C++ reader's bilinear convention) — the common physical zero;
2. `--e-shift S` — ONE shared positivity constant for the whole set,
   chosen from `--probe-shift` output so every table's minimum lands
   O(1)–O(10) code units above zero (never epsilon-positive: the
   log₁₀(e) inverse axis, the γₑ overflow guard, and dead-cell
   resurrection all assume e = 0 sits comfortably below every hull).
   A forced shift that leaves `min(e) <= 0` hard-fails — there is no
   silent per-table top-up, by design.

`make_eref295_set.sh` runs the whole recipe for the solid/interior/air set
(probe pass → shared S → emit + QA → design-point spread check, target
< 0.15 code units). The emitted `e_ref_state:` header line records the
anchor; e reads exactly `e_shift` at that state.

Extending an existing set: new members must adopt the set's S *verbatim*
(never recompute it — that would re-gauge the whole set); the only gate
left is that their post-reference minima still land in the O(1)–O(10)
band under the fixed S. `make_eref295_solids.sh` does this for the
aluminum (3720, cold-extended — note the material-scaled `--fit-rho`
Vinet window and `--z-cond` valence) and diamond (7834, native-hull)
solid alternates. A member may also be *re-emitted* on the same gauge
with its axes extended (e.g. `air_lowT_extension/make_air_lowT.sh`,
`--T-floor 1.25`): the gauge couples the set only through S and each
table's own 295 K reference, so as long as the native nodes are
byte-identical (gated) the other members need not be re-emitted.

## 4. Checking what you produced

The tool self-gates: the emitted surface must pass the **G1 scan** (zero
decreasing p intervals in ρ and in T, p > 0 everywhere) or the run aborts.
Beyond that:

- read the two `G1[...]` lines (raw vs emitted) and the `crossover:` /
  `LOO` lines against §3;
- look at the QA plots in `qa/`;
- run the C++ one-zone self-test on every table in `data/` (DEBUG build;
  see `run` in this directory once WS4 lands): reader integrity,
  e/p round-trips, derivative identities, hull behaviour.

## 5. Using the table in Cerberus

```lua
states = {
  fluid = {
    type = 'hydro',
    gas = {
      type  = 'tabulated',
      table = 'data/copper_3337_s311.eostab',
      mass  = 1.0, charge = 0.0,
      -- optional: ttol (1e-10), max_newton (100)
    },
    -- flux defaults to 'HLLC_general_eos'
  },
}
```

Reference quantities must be set (the gas model nondimensionalises the
table at load). See `doc/eos_table_reader.md` §9–10.

## 6. What these tables are NOT

- **Not the mixture-set (Track P) tables.** Tables destined for the
  Ti/D/air mixture must share the exact common T-bracket and the low-ρ/
  low-T extensions (plan §3 stage 5, requirement H3) — that pipeline
  (WS2e) builds *on* these extractions but is separate. Mixing
  arbitrary Track-S tables in one `tabulated_mixture` state will abort
  at the T-hull intersection check.
- Not a two-temperature closure (304/305 extraction exists for the
  future 2-T track only).
- Not valid below the coldest kept isotherm or outside the source ρ
  span — the runtime clamps (never extrapolates) and counts. With
  `--T-floor` the rows below the source floor exist but are construction
  (§3 item 13, hull 0): the rail moves, it does not disappear.
- Not a substitute for the published SESAME data in any context where
  the §3 replacements matter (dome equilibrium, tension, melt-line
  structure).
