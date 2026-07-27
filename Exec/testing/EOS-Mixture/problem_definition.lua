-- ======== SOD SHOCK TUBE: DALTON MIXTURE VALIDATION (Stage 9, W28) =========
--
-- The Sod problem run in five modes by the run script (GAS_MODE is set in
-- each .inputs file just before this file loads), designed to isolate
-- mixing-rule bugs from table bugs (doc/eos_mixture_dalton_plan.md sect 11):
--
--   single       TabulatedEOS, one gamma=1.4 synthetic table   (baseline)
--   mix_pure1    tabulated_mixture, TWO IDENTICAL gamma=1.4 tables,
--                alpha = 1 everywhere (all mass in component 0): every cell
--                takes the pure-cell short-circuit, which is required to be
--                FORMULA-IDENTICAL to the single-table path -> the whole run
--                must match 'single' to round-off (gate M1, pure path)
--   mix_pure0    as mix_pure1 with alpha = 0: all mass in the DERIVED last
--                component (1 - sum alpha) -> catches last-slot indexing bugs
--   mix_uniform  as mix_pure1 with alpha = 0.3 everywhere: every cell takes
--                the full Dalton solve. Dalton is analytically EXACT for
--                identical ideal-gas components (sum alpha_k rho R T = rho R T),
--                so the difference vs 'single' measures pure table-
--                interpolation error at the partial densities (gate M1b,
--                interpolation-level tolerance)
--   mix_binary   TWO DIFFERENT tables (gamma=1.4 left material, gamma=5/3
--                right material), alpha step at the diaphragm: the real
--                two-material tube. Gates: stability, alpha in [0,1], a
--                mixed contact layer, per-component mass + total energy
--                conservation, wall-time budget (M3/M4/M5)
--
-- The mixture runs carry `self_test` (plan M2): the ctor sweeps driver
-- round-trips, pure-cell bitwise parity, the frozen sound speed vs an
-- isentropic finite difference, and the drop_tol crossing, printing
-- MIXEOS-SELFTEST lines that check.py greps from the run logs.
--
-- Units as in EOS-Sod-Ideal: table rho-hull [1e-4, 10] g/cc, T [1e3, 1e7] K;
-- with ref_density = 1000 kg/m^3 and ref_temp = 1e5 K the Sod states sit
-- comfortably inside every hull.

-- === REFERENCE QUANTITIES ===

ref_length = 1.0
ref_density = 1000.0        -- kg/m^3 -> code rho=1 is 1 g/cc
ref_mass = 3.34358e-27      -- kg (deuteron mass, matching the tables)
ref_temp = 1e5              -- K

-- === SETTINGS ===

verbosity = 1
cfl = 0.5
time_integration_scheme = 'RK2'

-- === PROBLEM ===

local mix_components_identical = {
  { table = 'data/mix_a.eostab', name = 'mat_a', mass = 1.0, charge = 0.0 },
  { table = 'data/mix_a.eostab', name = 'mat_b', mass = 1.0, charge = 0.0 },
}

-- B-tier component lists (nonlinear tables; doc/eos_amagat_plan.md sect 5).
-- Paths relative to this case dir; run_amagat unpacks the committed .gz.
local spliced_dir = '../EOS-Table/data'
local b1_components = {
  { table = spliced_dir .. '/D_spliced.eostab', name = 'd_a', mass = 1.0, charge = 0.0 },
  { table = spliced_dir .. '/D_spliced.eostab', name = 'd_b', mass = 1.0, charge = 0.0 },
}
-- B2 pair: condensed-hulled D + dilute air. Ti was the intended heavy
-- component but Ti_spliced v0 isotherms are NON-MONOTONE in rho (vapor-wedge
-- fill vs tension clip: a 6-decade sawtooth), which makes the Amagat inner
-- rho-inversion ill-posed — the Ti variant is BLOCKED on the R0 monotone
-- reconditioning (doc/eos_amagat_plan.md section 7 finding, 2026-07-26).
local b2_components = {
  { table = spliced_dir .. '/D_spliced.eostab', name = 'light', mass = 1.0,  charge = 0.0 },
  { table = spliced_dir .. '/air_g14.eostab',   name = 'air',  mass = 14.38, charge = 0.0 },
}

local mix_components_binary = {
  { table = 'data/mix_a.eostab', name = 'left_g14',  mass = 1.0, charge = 0.0 },
  { table = 'data/mix_b.eostab', name = 'right_g53', mass = 1.0, charge = 0.0 },
}

function alpha_step(dat)
  if dat['x'] < 0.0 then
    return 1.0
  else
    return 0.0
  end
end

alpha_def = nil
flux_def = nil

if GAS_MODE == 'single' then
  gas_def = {
    type = 'tabulated',
    table = 'data/mix_a.eostab',
    mass = 1.0,
    charge = 0.0,
  }
elseif GAS_MODE == 'mix_pure1' then
  gas_def = {
    type = 'tabulated_mixture',
    mixing_rule = 'dalton',
    components = mix_components_identical,
    self_test = 12,
  }
  alpha_def = 1.0
elseif GAS_MODE == 'mix_pure0' then
  gas_def = {
    type = 'tabulated_mixture',
    mixing_rule = 'dalton',
    components = mix_components_identical,
  }
  alpha_def = 0.0
elseif GAS_MODE == 'mix_uniform' then
  gas_def = {
    type = 'tabulated_mixture',
    mixing_rule = 'dalton',
    components = mix_components_identical,
  }
  alpha_def = 0.3
elseif GAS_MODE == 'mix_binary' then
  gas_def = {
    type = 'tabulated_mixture',
    mixing_rule = 'dalton',
    components = mix_components_binary,
    self_test = 12,
  }
  alpha_def = alpha_step
-- Amagat twins (W27 plan Tier A, doc/eos_amagat_plan.md section 5):
-- identical numerics/ICs, only the mixing rule differs. amg_pure1 carries
-- the ctor self-test so the Amagat rp sweep runs on the identical-table
-- pair (ideal-gas: analytically exact reference).
elseif GAS_MODE == 'amg_pure1' then
  gas_def = {
    type = 'tabulated_mixture',
    mixing_rule = 'amagat',
    components = mix_components_identical,
    self_test = 12,
  }
  alpha_def = 1.0
elseif GAS_MODE == 'amg_pure0' then
  gas_def = {
    type = 'tabulated_mixture',
    mixing_rule = 'amagat',
    components = mix_components_identical,
  }
  alpha_def = 0.0
elseif GAS_MODE == 'amg_uniform' then
  gas_def = {
    type = 'tabulated_mixture',
    mixing_rule = 'amagat',
    components = mix_components_identical,
  }
  alpha_def = 0.3
elseif GAS_MODE == 'amg_binary' then
  gas_def = {
    type = 'tabulated_mixture',
    mixing_rule = 'amagat',
    components = mix_components_binary,
    self_test = 12,
  }
  alpha_def = alpha_step
-- B1: identical NONLINEAR (D_spliced) tables — the W28a gate, exact only
-- under Amagat; the Dalton twin is the control that measures the
-- shared-volume limitation on a real table
elseif GAS_MODE == 'b1_single' then
  gas_def = {
    type = 'tabulated',
    table = spliced_dir .. '/D_spliced.eostab',
    mass = 1.0,
    charge = 0.0,
  }
elseif GAS_MODE == 'b1_dalton' then
  gas_def = {
    type = 'tabulated_mixture',
    mixing_rule = 'dalton',
    components = b1_components,
  }
  alpha_def = 0.3
elseif GAS_MODE == 'b1_amg' then
  gas_def = {
    type = 'tabulated_mixture',
    mixing_rule = 'amagat',
    components = b1_components,
  }
  alpha_def = 0.3
-- B2/C2: Ti + air dilute-condensed stress (inputs set RHO_L... + ALPHA_TI)
elseif GAS_MODE == 'b2_dalton' then
  gas_def = {
    type = 'tabulated_mixture',
    mixing_rule = 'dalton',
    components = b2_components,
  }
  alpha_def = ALPHA_TI
elseif GAS_MODE == 'b2_amg' then
  gas_def = {
    type = 'tabulated_mixture',
    mixing_rule = 'amagat',
    components = b2_components,
    self_test = 8,  -- the heterogeneous REAL pair (differing hulls)
  }
  alpha_def = ALPHA_TI
-- C1: dilute ideal binary — both rules analytically exact, must agree
-- A4b: the TPG ANALYTIC reference for the binary tube — the active
-- multi-component thermally-perfect gas mixes gamma by mass-fraction-
-- weighted specific heats (MFP_thermally_perfect_gas.cpp:52-84), which is
-- EXACTLY what both Dalton and Amagat reduce to for ideal gases: an
-- interpolation-free ground truth. flux is pinned to HLLC_general_eos so
-- the SOLVER is identical across all three runs (round-off-identical to
-- HLLC on a gamma-law gas, STAGE5.md G1) and the comparison isolates the
-- closure.
elseif GAS_MODE == 'tpg_binary' then
  gas_def = {
    type = 'thermally_perfect',
    names = { 'left_g14', 'right_g53' },
    mass = { 1.0, 1.0 },
    charge = { 0.0, 0.0 },
    gamma = { 1.4, 1.6666666666666667 },
  }
  alpha_def = alpha_step
  flux_def = 'HLLC_general_eos'
elseif GAS_MODE == 'c1_dalton' then
  gas_def = {
    type = 'tabulated_mixture',
    mixing_rule = 'dalton',
    components = mix_components_binary,
  }
  alpha_def = 1.0e-6
elseif GAS_MODE == 'c1_amg' then
  gas_def = {
    type = 'tabulated_mixture',
    mixing_rule = 'amagat',
    components = mix_components_binary,
  }
  alpha_def = 1.0e-6
else
  error("GAS_MODE must be one of 'single', 'mix_pure1', 'mix_pure0', " ..
        "'mix_uniform', 'mix_binary', 'amg_pure1', 'amg_pure0', 'amg_uniform', 'amg_binary', ... " ..
        "(set in the .inputs file)")
end

-- B/C-tier state overrides: inputs may set RHO_L/P_L/RHO_R/P_R/ALPHA_TI in
-- the mfp.lua block BEFORE the dofile; defaults reproduce the classic tube
RHO_L = RHO_L or 1.0
P_L = P_L or 1.0
RHO_R = RHO_R or 0.125
P_R = P_R or 0.1
ALPHA_TI = ALPHA_TI or 1.0e-4

function rho0(dat)
  if dat['x'] < 0.0 then
    return RHO_L
  else
    return RHO_R
  end
end

function p0(dat)
  if dat['x'] < 0.0 then
    return P_L
  else
    return P_R
  end
end

states = {
  fluid = {
    type = 'hydro',
    gas = gas_def,
    flux = flux_def,  -- nil (key absent) except the tpg_binary reference
    reconstruction = 'minmod',
    -- no 'flux' key: the config default must select HLLC_general_eos for
    -- both the tabulated and the tabulated_mixture gas (the W24 predicate)
    value = {
      rho = rho0,
      x_vel = 0,
      y_vel = 0,
      z_vel = 0,
      p = p0,
      alpha = alpha_def,  -- nil for 'single' -> key absent -> no tracers
    },
  },
}

actions = {
  fluxes = {
    type = 'CTU',
    corner_transport = true,
    states = { 'fluid' },
  },
}
