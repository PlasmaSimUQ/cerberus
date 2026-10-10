-- 0-D (uniform, periodic, no fluxes) checks of the local ICF source terms:
--   CASE = 'brems'    bremsstrahlung cooling of a 50/50 DT plasma
--   CASE = 'fusion'   DT burn with fuel accounting, no alpha heating
--   CASE = 'ignition' burn + bremsstrahlung, all alpha energy deposited;
--                     T0_KEV below/above the P_alpha = P_br crossing (~4.6 keV)
-- CASE and T0_KEV are set from the inputs file. See check.py for the gates.

-- ======== REFERENCE QUANTITIES ==========
-- chosen so that code T is in keV and code mass is in deuteron masses
ref_length = 1.0e-4                                       -- m
ref_density = 1000.0                                      -- kg/m^3  (n_ref ~ 3e29 m^-3)
ref_mass = 3.34358e-27                                    -- kg, deuteron
ref_temp = 1.0e3 * 1.602176634e-19 / 1.38064852e-23       -- K, 1 keV

-- ======== PROBLEM ==========
m_D, m_T, m_He4 = 1.0, 1.497877, 1.987289 -- code units (deuteron masses)
rho0 = 1.0
T0 = T0_KEV or 2.0

-- 50/50 D/T by number, no ash: mass fractions
alpha_D = m_D / (m_D + m_T)
alpha_T = m_T / (m_D + m_T)
n0 = rho0 / (0.5 * (m_D + m_T)) -- total ion number density
p0 = n0 * T0                    -- p = n T in code units

-- ======== SETTINGS ==========
verbosity = 1
time_integration_scheme = 'strang'
force_dt = FORCE_DT or 1.0e-3

-- ======== STATES ==========
states = {
  plasma = {
    type = 'hydro',
    gas = {
      type = 'thermally_perfect',
      mass = { m_D, m_T, m_He4 },
      charge = { 1.0, 1.0, 2.0 },
      gamma = { 5 / 3, 5 / 3, 5 / 3 },
      names = { 'D', 'T', 'He4' },
    },
    reconstruction = 'constant',
    flux = 'HLLC', -- required key; never evaluated (no CTU action)
    value = {
      rho = rho0,
      x_vel = 0,
      y_vel = 0,
      z_vel = 0,
      p = p0,
      alpha = { alpha_D, alpha_T },
    },
  },
}

-- ======== ACTIONS ==========
brems = {
  type = 'bremsstrahlung',
  states = { 'plasma' },
  species = { 'D', 'T', 'He4' },
}

burn = {
  type = 'fusion',
  states = { 'plasma' },
  species = { 'D', 'T', 'He4' },
  fuel = { 'D', 'T' },
  product = 'He4',
  alpha_fraction = (CASE == 'ignition') and 1.0 or 0.0,
}

if CASE == 'brems' then
  actions = { brems = brems }
elseif CASE == 'fusion' then
  actions = { burn = burn }
elseif CASE == 'ignition' then
  actions = { brems = brems, burn = burn }
else
  error("CASE must be 'brems', 'fusion' or 'ignition'")
end
