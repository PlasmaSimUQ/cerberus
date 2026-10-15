-- =============================================================================
-- 2D cylindrically converging two-fluid (ion / electron) DT implosion with burn
-- =============================================================================
--
-- The hydrodynamic layout is that of Exec/testing/2D_RMI_implosion_test; the
-- plasma model, numerics and non-dimensional plasma parameters are those of
-- Exec/testing/Thermal-RMI (five-moment two-fluid + Maxwell, plasma5 coupling,
-- elliptic divergence cleaning, reduced light speed and electron mass).
--
-- Radial layout (R = distance from the origin):
--
--   R < r_inner                  : light  target, rho_i = density_inner (= 1)
--   r_inner <= R < r_outer       : heavy  shell,  rho_i = density_shell (= 2)
--   r_outer <= R <= r_background : post-shock gas driving inwards (M = shock_mach)
--   R > r_background             : ambient low-density background, at rest
--
-- Both fluids are 50/50 D-T by number, fully ionised, with T_e = T_i and
-- n_e = n_i everywhere initially; the He4 ash fraction starts at zero. The
-- `fusion` action (D + T -> He4 + n, NRL reactivity) consumes D and T into
-- He4 through the ion-state tracers and deposits the alpha energy locally,
-- split ion/electron by the Fraley partition. The neutron energy is lost.
--
-- The r_inner contact carries an azimuthal sinusoidal perturbation of
-- amplitude `a_perturb` and arc-length wavelength `lambda_perturb`.  The
-- inward-running shock crosses that perturbed contact and deposits vorticity
-- -> Richtmyer-Meshkov growth, on top of the Bell-Plesset convergence effect.
--
-- The shock at R = r_outer runs INWARDS into the (stationary) heavy shell.
-- Both fluids share the velocity and temperature, so the post-shock state is
-- the normal-shock jump of the COMBINED gas (p_i + p_e, rho_i + rho_e)
-- applied to each fluid's density and pressure separately.
--
-- Only the first quadrant is solved: x = 0 and y = 0 are symmetry planes and
-- carry reflecting boundary conditions for the fluids and the fields, the
-- outer x/y boundaries are outflow.  It requires the perturbation to be
-- mirror-symmetric about both planes, i.e. an EVEN azimuthal mode number.

-- ======== CASE SWITCHES ==========
-- set from the inputs file before dofile(); defaults: burn on, no losses
if BURN == nil then
  BURN = true
end
if BREMS == nil then
  BREMS = false
end

-- ======== REFERENCE QUANTITIES ==========
-- Thermal-RMI plasma settings: `lightspeed = 20` (non-dimensional) with the
-- physical reference light speed, so T_ref = (c / lightspeed)^2 m_ref / kB
-- (~5.4e10 K) is fixed by the code (`ref_temp` must then not be set). The
-- fusion action converts to keV through T_ref, so this is only a question of
-- what code T = 1 means; temperatures below are given in keV via T_keV.
kB, e_charge, c_SI = 1.38064852e-23, 1.602176634e-19, 299792458.0
ref_length = 1.0e-4                       -- m: r_inner = 1 is a 100 um target
ref_density = 1000.0                      -- kg/m^3 at rho = 1 (n_ref ~ 3e29 m^-3)
ref_mass = 3.34358e-27                    -- kg, deuteron
lightspeed = 20.0                         -- Thermal-RMI value (code units)
T_ref = (c_SI / lightspeed) ^ 2 * ref_mass / kB -- K, as derived by the code
T_keV = 1.0e3 * e_charge / (kB * T_ref)   -- code T per keV

-- ======== PROBLEM PARAMETERS ==========

r_inner = 1.0      -- radius of the perturbed material interface
r_outer = 1.1      -- radius of the incoming shock
r_background = 1.3 -- outer edge of the driving (post-shock) shell

density_inner = 1.0 -- ion rho, R < r_inner
density_shell = 2.0 -- ion rho, r_inner < R < r_outer
T_ambient = 0.1 * T_keV -- 0.1 keV, uniform ahead of the shock (both fluids)

-- ambient background beyond the driver, at rest, in pressure equilibrium with
-- the unshocked target/shell (so the only discontinuity at R = r_background
-- is in density): hot and thin, it sets the hydro CFL time step.
density_background = 0.1 -- ion rho, R > r_background

shock_mach = 10.0 -- shock Mach number w.r.t. the shell state (combined gas)

gam = 5 / 3

-- DT composition: masses in deuteron masses, 50/50 by number, no ash
m_D, m_T, m_He4 = 1.0, 1.497877, 1.987289
alpha_D = m_D / (m_D + m_T)
alpha_T = m_T / (m_D + m_T)
m_mix = 0.5 * (m_D + m_T)       -- ion mass per particle
mass_electron = 0.01            -- Thermal-RMI reduced electron mass

-- interface perturbation.  `lambda_perturb` is the wavelength measured as
-- arc length along the unperturbed r_inner circle; it is snapped to the
-- nearest EVEN azimuthal mode -- an integer number of wavelengths has to fit
-- around the circle, and cos(n*theta) is only mirror-symmetric about the x = 0
-- plane (theta -> pi - theta) when n is even, which the quarter domain needs.
lambda_perturb = 2 * math.pi * r_inner / 12
a_perturb = 0.02

-- width over which the contact is smeared (avoids a stair-stepped circle)
interface_transition = 0.02

-- ======== DERIVED QUANTITIES ==========

-- nearest even mode >= 2 (see the note on the quarter domain above)
mode_perturb = 2 * math.max(1, math.floor(math.pi * r_inner / lambda_perturb + 0.5))
lambda_effective = 2 * math.pi * r_inner / mode_perturb

-- per-fluid state from the ion density and the common temperature:
-- n_i = rho_i / m_mix, n_e = n_i (Z = 1), p_s = n_s T, rho_e = m_e n_e
function n_ion(rho_i)
  return rho_i / m_mix
end
function electron_rho(rho_i)
  return mass_electron * n_ion(rho_i)
end

-- pre-shock (the heavy shell) COMBINED state
rho_pre = density_shell + electron_rho(density_shell)
p_pre = 2 * n_ion(density_shell) * T_ambient
a_pre = math.sqrt(gam * p_pre / rho_pre)

-- normal shock jump for a shock of Mach `shock_mach` running into the shell
p_ratio = 1 + ((2 * gam) / (gam + 1)) * (shock_mach ^ 2 - 1)
rho_ratio = 1 / (1 - (2 / (gam + 1)) * (1 - 1 / shock_mach ^ 2))
p_post = p_pre * p_ratio
rho_post = rho_pre * rho_ratio
a_post = math.sqrt(gam * p_post / rho_post)

-- lab-frame speed imparted to the gas by the shock (positive = along the
-- shock's direction of travel, i.e. inwards here)
u_post = shock_mach * a_pre
  - a_post * math.sqrt(((gam - 1) * shock_mach ^ 2 + 2) / (2 * gam * shock_mach ^ 2 - (gam - 1)))

-- background temperature from pressure equilibrium with the shell
T_background = T_ambient * density_shell / density_background

-- planar estimates (the real shock strengthens as it converges, so these are
-- conservative): the shock reaches the axis at t_focus, and the rarefaction
-- launched at r_background catches it at t_unloaded.  Want t_unloaded > t_focus.
u_shock = shock_mach * a_pre -- inward shock speed
u_rarefaction = u_post + a_post -- inward speed of the rarefaction head
t_focus = r_outer / u_shock
t_unloaded = (r_background - r_outer) / (u_rarefaction - u_shock)

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

-- ======== INITIAL CONDITIONS ==========

-- ion density by zone (the shock jump applies to the pre-shock SHELL state)
function ion_density(dat)
  local x, y = dat['x'], dat['y']
  local r = math.sqrt(x * x + y * y)

  if r > r_background then
    return density_background
  end
  if r >= r_outer then
    return density_shell * rho_ratio
  end

  -- smoothed light/heavy contact
  return blend(r - contact_radius(x, y), density_inner, density_shell, interface_transition)
end

function electron_density(dat)
  return electron_rho(ion_density(dat))
end

-- temperature by zone, common to both fluids
function temperature(dat)
  local x, y = dat['x'], dat['y']
  local r = math.sqrt(x * x + y * y)

  if r > r_background then
    return T_background
  end
  if r >= r_outer then
    return T_ambient * p_ratio / rho_ratio
  end
  return T_ambient
end

function ion_pressure(dat)
  return n_ion(ion_density(dat)) * temperature(dat)
end

function electron_pressure(dat)
  return n_ion(ion_density(dat)) * temperature(dat) -- n_e = n_i
end

-- only the driver shell moves; the target, the shell and the background are
-- all initially at rest.  Both fluids share the velocity.
function velocity_x(dat)
  local x, y = dat['x'], dat['y']
  local r = math.sqrt(x * x + y * y)

  if r < r_outer or r > r_background or r < 1e-14 then
    return 0.0
  end
  return -u_post * x / r
end

function velocity_y(dat)
  local x, y = dat['x'], dat['y']
  local r = math.sqrt(x * x + y * y)

  if r < r_outer or r > r_background or r < 1e-14 then
    return 0.0
  end
  return -u_post * y / r
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

  ion = {
    type = 'hydro',
    gas = {
      type = 'thermally_perfect',
      mass = { m_D, m_T, m_He4 },
      charge = { 1.0, 1.0, 2.0 },
      gamma = { gam, gam, gam },
      names = { 'D', 'T', 'He4' },
    },
    reconstruction = 'minmod',
    refinement = { name = 'hydro_gradient', rho = 0.1 },
    flux = 'HLLC',
    value = {
      rho = ion_density,
      x_vel = velocity_x,
      y_vel = velocity_y,
      p = ion_pressure,
      alpha = { alpha_D, alpha_T }, -- He4 is the implicit last component
    },
    bc = hydro_bc,
  },

  electron = {
    type = 'hydro',
    gas = {
      type = 'thermally_perfect',
      mass = mass_electron,
      charge = -1.0,
      gamma = gam,
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
  },

  field = {
    type = 'field',
    reconstruction = 'O6',
    flux = 'RankineHugoniot',
    bc = {
      x = { lo = { fill_D_bc = 'symmetry', fill_B_bc = 'symmetry' } },
      y = { lo = { fill_D_bc = 'symmetry', fill_B_bc = 'symmetry' } },
    },
  },
}

-- ======== DEFINE ACTIONS ==========

actions = {

  hydro_fluxes = {
    type = 'CTU',
    corner_transport = true,
    states = { 'ion', 'electron' },
  },

  em_fluxes = {
    type = 'CTU',
    corner_transport = true,
    states = { 'field' },
  },

  plasma = {
    type = 'plasma5',
    solver = 'explicit',
    states = { 'ion', 'electron', 'field' },
  },

  divergence_cleaning = {
    type = 'elliptic',
    projection = 1,
    state = 'field',
  },

}

-- D + T -> He4 + n: fuel accounting through the ion tracers, alpha energy
-- deposited locally with the Fraley ion/electron split (T_e from the
-- electron fluid), neutron energy lost.
if BURN then
  actions.burn = {
    type = 'fusion',
    states = { 'ion', 'electron' },
    species = { 'D', 'T', 'He4', 'electron' },
    fuel = { 'D', 'T' },
    product = 'He4',
    alpha_fraction = 1.0,
    alpha_partition = 'fraley',
  }
end

-- optically-thin bremsstrahlung loss from the electron fluid
if BREMS then
  actions.brems = {
    type = 'bremsstrahlung',
    states = { 'ion', 'electron' },
    species = { 'D', 'T', 'He4', 'electron' },
  }
end
