-- Checks of the `laser` action (inverse-bremsstrahlung ray march, no
-- refraction) against closed forms. No flux action: the plasma is frozen
-- apart from the deposited energy, so the first source evaluation sees the
-- initial fields exactly (gates read the [laser] log line; see check.py).
--   CASE = 'ramp1d'  plane wave into a linear density ramp up to n_c
--   CASE = 'ramp2d'  single ray into the same ramp at 30 deg incidence
--   CASE = 'beam2d'  top-hat beam through a uniform absorber

-- ======== CONSTANTS / REFERENCE QUANTITIES ==========
eps0, m_e, e, c, kB, m_p = 8.8541878128e-12, 9.1093837015e-31, 1.602176634e-19, 299792458.0, 1.38064852e-23, 1.67262192e-27
lambda = 0.351e-6
n_crit = eps0 * m_e * (2 * math.pi * c / lambda) ^ 2 / e ^ 2 -- m^-3

ref_length = 1.0e-4             -- m
ref_mass = m_p                  -- kg: code mass 1 = proton, Z = 1 -> n_e = n_i
ref_density = n_crit * m_p      -- kg/m^3: code n = 1 is the critical density
ref_temp = 1.0e3 * e / kB       -- K: code T in keV

-- ======== PROBLEM ==========
L = 1.0           -- ramp length, code units (1e-4 m)
n_floor = 1.0e-3
T0 = 1.0          -- keV, uniform
I0 = 1.0e17       -- W/m^2
LNLAMBDA = 7.0

function n_profile(x)
  if CASE == 'beam2d' then
    return 0.5
  end
  if x < L then
    return math.max(x / L, n_floor)
  end
  return 1.2
end

function density(dat)
  return n_profile(dat['x'])
end
function pressure(dat)
  return n_profile(dat['x']) * T0
end

-- ======== SETTINGS ==========
verbosity = 2
time_integration_scheme = 'strang'
force_dt = 1.0e-6

-- ======== STATES ==========
states = {
  plasma = {
    type = 'hydro',
    gas = { type = 'thermally_perfect', mass = 1.0, charge = 1.0, gamma = 5 / 3 },
    reconstruction = 'constant',
    flux = 'HLLC', -- required key; never evaluated (no CTU action)
    value = { rho = density, x_vel = 0, y_vel = 0, z_vel = 0, p = pressure },
  },
}

-- ======== ACTIONS ==========
theta = math.rad(30)
beams = {
  ramp1d = { origin = { -0.1e-4, 0, 0 }, direction = { 1, 0, 0 }, radius = 0 },
  ramp2d = { origin = { -0.1e-4, 0.1e-4, 0 }, direction = { math.cos(theta), math.sin(theta), 0 }, radius = 0 },
  beam2d = { origin = { -0.1e-4, 0.5e-4, 0 }, direction = { 1, 0, 0 }, radius = 0.2e-4 },
}
b = beams[CASE] or error("CASE must be 'ramp1d', 'ramp2d' or 'beam2d'")

actions = {
  drive = {
    type = 'laser',
    state = 'plasma',
    origin = b.origin,
    direction = b.direction,
    radius = b.radius,
    profile = 'tophat',
    rays_per_cell = 2,
    wavelength = lambda,
    intensity = I0,
    lnLambda = LNLAMBDA,
    critical_fraction = 1.0,
  },
}
