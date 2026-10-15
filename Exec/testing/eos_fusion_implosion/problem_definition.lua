-- =============================================================================
-- 2D cylindrically converging implosion, single fluid, tabulated SESAME EOS
-- =============================================================================
--
-- The eos_fusion_implosion series, selected by flags set from the inputs
-- file (defaults = step 1):
--   step 1  TWO_FLUID = false, BURN = false   hydro only, tabulated mixture
--   step 2  TWO_FLUID = true                  + electron fluid, Maxwell field,
--                                             plasma5 coupling (Thermal-RMI
--                                             plasma settings, as in
--                                             Exec/testing/fusion_implosion)
--   step 3  TWO_FLUID = true,  BURN = true    + D-T target, He4 ash, `fusion`
--                                             action; T and He4 are isotope-
--                                             scaled stand-ins of the D table
--                                             (make_isotope_tables.py)
-- DRIVER_VELOCITY / DRIVER_TEMPERATURE (SI) override the prescribed driver;
-- CASE_DIR locates the tables when run from a step subdirectory.
--
-- The hydrodynamic layout is that of Exec/testing/fusion_implosion (itself
-- from 2D_RMI_implosion_test), closed by the `tabulated_mixture` gas model on
-- real SESAME-derived tables of the common-gauge eref295 set
-- (Exec/testing/EOS-Sesame/data, see that README):
--
--   deuterium   5267  table 301   rho 1e-10..1e3 g/cc,  T 4 K..1e9 K
--   diamond     7834  table 301   rho 3.5e-6..7e4 g/cc, T 72 K..1.2e9 K
--
-- The mixture's shared temperature bracket is the intersection of the
-- component hulls, 72 K .. 1e9 K. (Dry air 5031 was the background in the
-- first version; at 0.1 g/cc and 5000 K its pressure is inside the table's
-- interpolation noise, giving cells of either sign, and it would have capped
-- the bracket at 30 keV. The background is now a dilute deuterium gas.)
--
-- Radial layout (R = distance from the origin), materials by zone:
--
--   R < r_inner                  : deuterium (step 3: 50/50 D-T) target, rho_target, T_ambient
--   r_inner <= R < r_outer       : diamond shell,    rho_shell,  T_ambient
--   r_outer <= R <= r_background : diamond driver, compressed, hot, moving in
--   R > r_background             : deuterium gas background, rho_background, T_ambient
--
-- Every zone is initialised by (rho, T) and the tables supply the pressure:
-- with condensed-matter tables a pressure-specified initial state is
-- ill-conditioned (GPa per 100 K), and the sub-GPa pressure mismatches
-- between the zones at the ambient temperature are negligible against the
-- driver (TPa).
--
-- Mixing rule and ambient state (found by trial, 2026-10-12):
--   * `amagat`, not `dalton`: under Dalton a trace of diamond in an air cell
--     evaluates the diamond table at a partial density below its hull and
--     clamps, which poisoned the whole background (negative pressures at
--     t = 0, blow-up at the driver edge). Amagat evaluates each component
--     at its own density on its own hull. `drop_tol`/`pure_tol` discard
--     round-off traces of a component before any evaluation.
--   * T_ambient = 5000 K, not 295 K: at 295 K the air table's pressure
--     (~1-10 MPa) is below its own interpolation noise (~50 MPa), giving
--     cells of either sign; at 5000 K the thermal pressure dominates.
--   * rho_background = 0.01 g/cc deuterium gas (700:1 driver/background):
--     the D table is the best-conditioned of the set (95 % hull, smooth);
--     see the note on air above.
--   * Steps 2-3 run on a UNIFORM 128^2 grid (amr.max_level = 0). With AMR,
--     the conservative prolongation of (rho, rho E) across the driver /
--     background jump at each regrid produces cells whose specific energy
--     lies below the table floor of the eref295 gauge (e ~ 40 % low at the
--     right rho and T): an ideal gas just turns cold there, a gauged table
--     returns a negative pressure. Step 1's weaker driver stayed above the
--     floor; the 300 km/s, 1e6 K driver of steps 2-3 does not.
--
-- The driver is a PRESCRIBED state (density ratio, temperature, inward
-- velocity), not the ideal-gas shock jump of the parent cases: with a real
-- EOS the jump is the table's Hugoniot, which the hydro computes itself.
--
-- Two-fluid closure (steps 2-3): the SESAME 301 tables are TOTAL
-- (ion + electron + cold) equations of state, used here to close the heavy
-- `fluid` state; the `electron` state is a thermally-perfect gas with the
-- Thermal-RMI reduced electron mass, singly ionised everywhere (n_e = sum
-- Z_k n_k with Z_k = 1 for D, T, C, air and 2 for He4). The electron
-- thermal energy is therefore counted in both fluids at high temperature;
-- the proper split is the SESAME 303/304 two-temperature closure
-- (MFP_eos_tab2T.H draft), out of scope for this series.
--
-- The r_inner contact carries an azimuthal sinusoidal perturbation of
-- amplitude `a_perturb` and arc-length wavelength `lambda_perturb`; the
-- inward-running shock crosses it and deposits vorticity -> Richtmyer-Meshkov
-- growth, on top of the Bell-Plesset convergence effect. Only the first
-- quadrant is solved (x = 0, y = 0 symmetry planes), which requires an EVEN
-- azimuthal mode number.

-- ======== CASE SWITCHES ==========
if TWO_FLUID == nil then
  TWO_FLUID = false
end
if BURN == nil then
  BURN = false
end
case_dir = CASE_DIR or '.'
driver_velocity_ms = DRIVER_VELOCITY or 5.0e4    -- m/s
driver_temperature_K = DRIVER_TEMPERATURE or 1.0e5 -- K

-- ======== REFERENCE QUANTITIES ==========
-- code mass in deuteron masses, code length 1e-4 m, code rho = 1 is 1 g/cc
-- (the tables' own density unit). The plasma settings are Thermal-RMI's:
-- `lightspeed = 20` (non-dimensional) with the physical reference light
-- speed, which fixes T_ref = (c / lightspeed)^2 m_ref / kB (~5.4e10 K); the
-- code reads this as `ref_temp` must then not be set. Temperatures below
-- are given in kelvin and converted with T_K.
kB, e_charge, c_SI = 1.38064852e-23, 1.602176634e-19, 299792458.0
ref_length = 1.0e-4                       -- m: r_inner = 1 is a 100 um target
ref_density = 1000.0                      -- kg/m^3 at rho = 1 (1 g/cc)
ref_mass = 3.34358e-27                    -- kg, deuteron
lightspeed = 20.0                         -- Thermal-RMI value (code units)
T_ref = (c_SI / lightspeed) ^ 2 * ref_mass / kB -- K, as derived by the code
u_ref = c_SI / lightspeed                 -- m/s
T_K = 1.0 / T_ref                         -- code T per kelvin

-- ======== MATERIALS ==========
-- masses in deuteron masses; charges = ionisation state used by the plasma
-- coupling (singly ionised approximation, see header; He4 fully ionised)
sesame_dir = case_dir .. '/../EOS-Sesame/data/'
standin_dir = case_dir .. '/data/'
amu_D = 2.014101
materials = {
  { table = sesame_dir .. 'diamond_7834_s301_eref295.eostab', name = 'C', mass = 12.011 / amu_D, charge = 1.0 },
}
if BURN then
  table.insert(materials, { table = standin_dir .. 'tritium_standin_from_D5267_eref295.eostab', name = 'T',   mass = 3.016049 / amu_D, charge = 1.0 })
  table.insert(materials, { table = standin_dir .. 'helium4_standin_from_D5267_eref295.eostab', name = 'He4', mass = 4.002602 / amu_D, charge = 2.0 })
end
table.insert(materials, { table = sesame_dir .. 'deuterium_5267_s301_eref295.eostab', name = 'D', mass = 1.0, charge = 1.0 })
-- component order = alpha order; D (last) is the implicit 1 - sum(alpha):
-- it fills both the target (with T in step 3) and the background

-- target composition by mass: D only, or 50/50 D-T by number (step 3);
-- D is the implicit component so only T needs a tracer
m_T = 3.016049 / amu_D
target_alpha = {}
if BURN then
  target_alpha = { T = m_T / (1.0 + m_T) }
end
-- Thermal-RMI reduced electron mass (deuteron masses). Finding
-- (2026-10-11, step 2 at 300 km/s): the electron pressure front leaves the
-- driver at the electron sound speed (~4 code units, decoupled from the
-- ions), reaches the outflow boundary by t ~ 0.13 and the electron fluid
-- runs away there (T_e > 10 keV everywhere by t = 0.5, flow speeds above
-- the reduced c). Not resolved; the mass ratio is not to be changed.
mass_electron = 0.01

-- ======== PROBLEM PARAMETERS ==========

r_inner = 1.0      -- radius of the perturbed material interface
r_outer = 1.1      -- inner edge of the driver
r_background = 1.3 -- outer edge of the driver

rho_target = 0.2          -- g/cc, deuterium (ice-like)
rho_shell = 3.515         -- g/cc, diamond at its reference density
rho_background = 0.01     -- g/cc, deuterium gas background (see header)
T_ambient = 5000.0 * T_K  -- every zone except the driver (see header)

-- prescribed driver: diamond compressed, heated and moving inwards
driver_density_ratio = 2.0
driver_temperature = driver_temperature_K * T_K -- K -> code
driver_velocity = driver_velocity_ms / u_ref    -- m/s -> code, inward

-- interface perturbation.  `lambda_perturb` is the wavelength measured as
-- arc length along the unperturbed r_inner circle; it is snapped to the
-- nearest EVEN azimuthal mode -- an integer number of wavelengths has to fit
-- around the circle, and cos(n*theta) is only mirror-symmetric about the x = 0
-- plane (theta -> pi - theta) when n is even, which the quarter domain needs.
lambda_perturb = 2 * math.pi * r_inner / 12
a_perturb = 0.02

-- width over which the contact is smeared (avoids a stair-stepped circle);
-- cells inside it are genuine D/C mixtures handled by the Dalton solve
interface_transition = 0.02

-- ======== DERIVED QUANTITIES ==========

-- nearest even mode >= 2 (see the note on the quarter domain above)
mode_perturb = 2 * math.max(1, math.floor(math.pi * r_inner / lambda_perturb + 0.5))
lambda_effective = 2 * math.pi * r_inner / mode_perturb

rho_driver = rho_shell * driver_density_ratio

-- planar estimate of the time for the driver front to reach the axis
t_focus = r_outer / driver_velocity

-- geometry sanity check
if not (r_inner < r_outer and r_outer < r_background) then
  error('require r_inner < r_outer < r_background')
end

-- ======== HELPERS ==========

-- Lua 5.3+ dropped math.tanh
function tanh(x)
  if x > 30.0 then
    return 1.0
  end
  if x < -30.0 then
    return -1.0
  end
  local e = math.exp(2.0 * x)
  return (e - 1.0) / (e + 1.0)
end

-- atanh(0.85): the blend below reaches 92.5% of the far-field value at
-- s = +/- width/2
blend_const = 2.0 * 1.2561551

-- smooth heaviside about s = 0: -> L for s << 0, -> R for s >> 0
function blend(s, L, R, width)
  if math.abs(L - R) < 1e-14 then
    return L
  end
  local t = tanh(blend_const * s / width)
  return 0.5 * ((L + R) - t * (L - R))
end

-- radius of the perturbed contact at the azimuthal angle of (x, y)
function contact_radius(x, y)
  local theta = math.atan(y, x)
  return r_inner + a_perturb * math.cos(mode_perturb * theta)
end

function radius(dat)
  local x, y = dat['x'], dat['y']
  return math.sqrt(x * x + y * y), x, y
end

-- ======== INITIAL CONDITIONS ==========

function density(dat)
  local r, x, y = radius(dat)
  if r > r_background then
    return rho_background
  end
  if r >= r_outer then
    return rho_driver
  end
  -- smoothed target/shell contact
  return blend(r - contact_radius(x, y), rho_target, rho_shell, interface_transition)
end

function temperature(dat)
  local r = radius(dat)
  if r >= r_outer and r <= r_background then
    return driver_temperature
  end
  return T_ambient
end

-- zone fractions: target (blended into the shell at the contact) and
-- shell + driver; the background is what remains
function target_fraction(dat)
  local r, x, y = radius(dat)
  if r >= r_outer then
    return 0.0
  end
  return blend(r - contact_radius(x, y), 1.0, 0.0, interface_transition)
end

function shell_fraction(dat)
  local r, x, y = radius(dat)
  if r > r_background then
    return 0.0
  end
  if r >= r_outer then
    return 1.0
  end
  return blend(r - contact_radius(x, y), 0.0, 1.0, interface_transition)
end

-- mass fraction of each explicit-tracer material: C in the shell/driver,
-- target materials by target_alpha (T in step 3), He4 absent; D is implicit
function material_fraction(name)
  if name == 'C' then
    return shell_fraction
  elseif target_alpha[name] then
    return function(dat)
      return target_alpha[name] * target_fraction(dat)
    end
  end
  return function(dat)
    return 0.0
  end
end

alpha_funcs = {}
for k = 1, #materials - 1 do
  alpha_funcs[k] = material_fraction(materials[k].name)
end

-- ion number density and (singly ionised) electron density for the plasma
function n_ion(dat)
  local rho = density(dat)
  local sum = 0.0
  local n = 0.0
  for k = 1, #materials - 1 do
    local a = alpha_funcs[k](dat)
    sum = sum + a
    n = n + a / materials[k].mass
  end
  n = n + (1.0 - sum) / materials[#materials].mass
  return rho * n
end

function n_electron(dat)
  local rho = density(dat)
  local sum = 0.0
  local n = 0.0
  for k = 1, #materials - 1 do
    local a = alpha_funcs[k](dat)
    sum = sum + a
    n = n + materials[k].charge * a / materials[k].mass
  end
  n = n + materials[#materials].charge * (1.0 - sum) / materials[#materials].mass
  return rho * n
end

function electron_density(dat)
  return mass_electron * n_electron(dat)
end

function electron_pressure(dat)
  return n_electron(dat) * temperature(dat) -- p = n T in code units
end

-- only the driver moves; the target, the shell and the background are at rest
function velocity_x(dat)
  local r, x, y = radius(dat)
  if r < r_outer or r > r_background or r < 1e-14 then
    return 0.0
  end
  return -driver_velocity * x / r
end

function velocity_y(dat)
  local r, x, y = radius(dat)
  if r < r_outer or r > r_background or r < 1e-14 then
    return 0.0
  end
  return -driver_velocity * y / r
end

-- ======== SETTINGS ==========

verbosity = 1
cfl = 0.3

time_integration_scheme = 'strang'

force_dt = 0

skin_depth = 1.0
beta = 1.0

-- ======== DEFINE STATES ==========

-- x = 0 and y = 0 are symmetry planes of the problem: reflecting.
-- the outer boundaries sit in the uniform background: outflow.
hydro_bc = {
  x = { lo = { fill_hydro_bc = 'symmetry' }, hi = { fill_hydro_bc = 'outflow' } },
  y = { lo = { fill_hydro_bc = 'symmetry' }, hi = { fill_hydro_bc = 'outflow' } },
}

states = {

  fluid = {
    type = 'hydro',
    gas = {
      type = 'tabulated_mixture',
      mixing_rule = 'amagat', -- see header
      drop_tol = 1.0e-3,      -- drop a component below this mass fraction
      pure_tol = 1.0e-6,      -- single-table path above 1 - pure_tol
      components = materials,
    },
    reconstruction = 'minmod',
    -- no 'flux' key: the config default selects the general-EOS HLLC
    refinement = { name = 'hydro_gradient', rho = 0.1 },
    value = {
      rho = density,
      x_vel = velocity_x,
      y_vel = velocity_y,
      T = temperature, -- pressure from the tables
      alpha = alpha_funcs,
    },
    bc = hydro_bc,
  },
}

if TWO_FLUID then
  states.electron = {
    type = 'hydro',
    gas = {
      type = 'thermally_perfect',
      mass = mass_electron,
      charge = -1.0,
      gamma = 5 / 3,
    },
    reconstruction = 'minmod',
    refinement = { name = 'hydro_gradient', rho = 0.1 },
    flux = 'HLLC',
    value = {
      rho = electron_density,
      x_vel = velocity_x,
      y_vel = velocity_y,
      p = electron_pressure,
    },
    bc = hydro_bc,
  }

  states.field = {
    type = 'field',
    reconstruction = 'O6',
    flux = 'RankineHugoniot',
    bc = {
      x = { lo = { fill_D_bc = 'symmetry', fill_B_bc = 'symmetry' } },
      y = { lo = { fill_D_bc = 'symmetry', fill_B_bc = 'symmetry' } },
    },
  }
end

-- ======== DEFINE ACTIONS ==========

hydro_states = TWO_FLUID and { 'fluid', 'electron' } or { 'fluid' }

actions = {

  hydro_fluxes = {
    type = 'CTU',
    corner_transport = true,
    states = hydro_states,
  },
}

if TWO_FLUID then
  actions.em_fluxes = {
    type = 'CTU',
    corner_transport = true,
    states = { 'field' },
  }

  actions.plasma = {
    type = 'plasma5',
    solver = 'explicit',
    states = { 'fluid', 'electron', 'field' },
  }

  actions.divergence_cleaning = {
    type = 'elliptic',
    projection = 1,
    state = 'field',
  }
end

-- D + T -> He4 + n: fuel accounting through the `fluid` tracers (the mixture
-- composition is thermodynamically active), alpha energy deposited locally
-- with the Fraley ion/electron split when the electron fluid exists
if BURN then
  burn_species = { 'D', 'T', 'He4' }
  if TWO_FLUID then
    table.insert(burn_species, 'electron')
  end
  actions.burn = {
    type = 'fusion',
    states = hydro_states,
    species = burn_species,
    fuel = { 'D', 'T' },
    product = 'He4',
    alpha_fraction = 1.0,
    alpha_partition = TWO_FLUID and 'fraley' or 'fixed',
    ion_fraction = 1.0,
  }
end
