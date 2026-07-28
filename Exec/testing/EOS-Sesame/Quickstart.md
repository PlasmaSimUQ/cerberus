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
| 306 | T=0 cold curve (1-D) | no — cannot populate a (ρ,T) table |
| 401/411/412 | vapor dome / melt boundaries | not used (v1) |
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
   Physically inert if every quantity comes from the same table.
9. **cv floor.** cv = ∂e/∂T is floored at 1e-3 · (3/2 k_B/m̄) (m̄ from the
   201 table); floored-cell count recorded.
9b. **Thermal stiffness floor.** The emitted `dpdrho` block is floored
    at k_B·T/m̄ everywhere — no single-phase fluid is isothermally softer
    than an ideal gas, so cells below it are tie-line residue the band
    mapping missed; they are floored, joined to the band, and hull-0
    (`thermal_floored=` count in the provenance line). Genuine
    near-critical softening above the thermal floor is left untouched.
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
11. **Alloy caveat.** 2963 is Ti-Beta-21S (an alloy): Z̄≈21.1, Ā≈45.9 from
    its 201 table are the *alloy* constants and are what the emitted
    `composition:` line carries. Do not expect pure-Ti transition
    pressures from it.

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
# the Ti casing table on the shared mixture bracket [17.78 K, 1e9 K]:
$PREP sesame --src $SES --mat 2963 --lT 1.25 9.0 384 \
    --lrho -5.30 1.69897 576 --c-cav 5.34 \
    --out data/ti-beta-21s_2963_trackP.eostab --qa qa
```

Cells outside the SESAME source span (e.g. below its 72.5 K T-floor)
are nearest-filled (constant-in-T) and marked hull-0 — an *extension*,
not data. Pin `--c-cav` for reproducibility across grids (the fallback
samples a grid-dependent cell; record the value you pin and why).

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
  span — the runtime clamps (never extrapolates) and counts.
- Not a substitute for the published SESAME data in any context where
  the §3 replacements matter (dome equilibrium, tension, melt-line
  structure).
