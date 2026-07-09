# Lua input options: tracer fluids, reconstruction fallback, effective zero

This is a short guide to three related options in the Lua problem definition
(`problem_definition.lua` or the Lua section of the inputs file) for `hydro`
states. A worked example using all three is
[`Exec/testing/Double-Rarefaction/problem_definition.lua`](../Exec/testing/Double-Rarefaction/problem_definition.lua).

## 1. Tracer fluids

A `hydro` state carries mass fractions *alpha* (one value per cell per tracer)
that track sub-components of the fluid. For a state with *N* components only
*N−1* alphas are stored — the last component's fraction is always computed as
`1 - sum(alpha)`.

Tracers are set through the `alpha` key of the state's `value` table (the
initial conditions). Each entry may be a constant or a function of position,
like any other `value` entry. There are two ways to use them:

### Passive tracers (single-species gas)

Keep the gas properties scalar and give `alpha` a value or a table. The number
of tracers is taken from the table length and the gas properties are duplicated
internally, so the tracers are purely passive "colour" fields — useful for
tagging regions, tracking mixing, and diagnostics.

``` lua
fluid = {
  type = 'hydro',
  gas = {
    type = 'thermally_perfect',
    mass = 1.0, charge = 0.0, gamma = 1.4,
  },
  value = {
    rho = 1.0, x_vel = 0, y_vel = 0, z_vel = 0, p = 0.4,
    alpha = tracer_fn,                 -- one tracer, or:
    -- alpha = {fn_1, fn_2, fn_3},     -- several tracers
  },
}
```

### Active multi-component gas

Give the gas `mass`, `charge` and `gamma` as equal-length lists (they must all
have the same length, or the run aborts). This defines *N* species and requires
exactly *N−1* `alpha` entries. Here the alphas are physically active: the
effective mass, charge and ratio of specific heats of each cell are
mixture-averages over the local fractions.

``` lua
fluid = {
  type = 'hydro',
  gas = {
    type = 'thermally_perfect',
    names  = {'deuterium', 'tritium', 'electron'}, -- optional labels for output
    mass   = {1.0, 1.5, 0.01},
    charge = {1.0, 1.0, -1.0},
    gamma  = {1.4, 1.4, 5/3},
  },
  value = {
    rho = 1.0, x_vel = 0, y_vel = 0, z_vel = 0, p = 0.4,
    alpha = {alpha_D, alpha_T},  -- N-1 entries; electron fraction = 1 - sum
  },
}
```

If `alpha` is omitted the tracers default to 0. Boundary condition blocks
follow the same convention: `alpha` in a BC may also be a table with one entry
per tracer.

## 2. Reconstruction fallback

High-order reconstruction schemes (how cell-average data is extrapolated to
cell faces before the flux calculation) can overshoot at steep gradients and
produce unphysical face values — negative density, pressure, or tracer
fractions. The optional per-state key `reconstruction_fallback` names a
lower-order scheme that is re-applied *only* to the individual face values that
come out unphysical:

``` lua
fluid = {
  type = 'hydro',
  reconstruction = 'O6',              -- primary (high order, unlimited)
  reconstruction_fallback = 'minmod', -- used only where the primary fails
  ...
}
```

Notes:

- The feature is **off by default**; the presence of the key enables it.
- The fallback must be a *lower-order* scheme than the primary (its stencil
  must fit inside the primary's), and cannot be `null`. Violations abort at
  startup with an explanatory message.
- The guarded quantities are density and pressure (checked against
  `effective_zero`, see below) and every tracer alpha (checked against exactly
  0, since `alpha = 0` is a legitimate value). For MHD states, density and
  pressure are guarded.
- The check and fallback are applied per face and per component, so untroubled
  regions are untouched and a run that never triggers the fallback is
  bit-identical to one without the key.
- If the fallback value is *still* unphysical, the face value degrades to the
  cell-centre value (first order).
- With `verbosity = 2` or higher, a summary line reports how many face values
  were repaired in each grid patch.

A face-positive reconstruction does not strictly guarantee a positive cell
average after the update (see Zhang & Shu, J. Comput. Phys. 229 (2010)
8918–8934), so tiny negative residuals can survive; the fallback is a
robustness measure, with the floors below as the hard guard.

## 3. Effective zero

`effective_zero` is an optional per-state key (default `1e-14`, in
non-dimensional code units) giving the smallest value of density and pressure
the state treats as physically meaningful:

``` lua
fluid = {
  type = 'hydro',
  effective_zero = 1e-14,
  ...
}
```

It is used as:

- the positivity threshold for density and pressure in the reconstruction
  fallback guard (tracers use 0 instead, as above);
- the floor applied to density, pressure and temperature during the
  conserved-to-primitive conversion;
- the floor on density and pressure inside the wave-speed (time step)
  calculation;
- the clamp applied to the left/right face states just before each Riemann
  solve.

The floor layers (everything except the fallback guard) are compiled in only
when the build flag `USE_PRIM_FLOOR=TRUE` is set — this is the default; setting
it to `FALSE` appends `.NOFLOOR` to the executable name and prints a warning at
build time.

Choose `effective_zero` relative to your problem's non-dimensional scales: it
should be far below any physically meaningful density or pressure in the run.
Note that near-vacuum states floored at a very small density while pressure
remains finite acquire an enormous sound speed, which can collapse the time
step — if a run stalls after flooring, revisit the floor value relative to the
problem scales.
