-- ======== PROBLEM ==========

-- we force the internal state to use the explicit mu0 by using
-- this reference temperature
kB = 1.38064852 * 10 ^ -23
ref_temp = 1 / kB

G = 1.0
t0 = 10.0
mu0 = 0.08
rho0 = 1.0

nu = mu0 / rho0

-- get the pressure profile
dofile('./pressure.lua')
require('./linear_interp')

function tangential_flow_velocity(r)
  return G / (2 * math.pi * r) * (1 - math.exp(-r ^ 2 / (4 * nu * t0)))
end

function flow_u(dat)
  local r = math.sqrt(dat['x'] ^ 2 + dat['y'] ^ 2)
  local u_t = tangential_flow_velocity(r)
  return u_t * dat['y'] / r
end

function flow_v(dat)
  local r = math.sqrt(dat['x'] ^ 2 + dat['y'] ^ 2)
  local u_t = tangential_flow_velocity(r)
  return -u_t * dat['x'] / r
end

function flow_p(dat)
  local r = math.sqrt(dat['x'] ^ 2 + dat['y'] ^ 2)
  return linear_interp(r, x_data_pressure, y_data_pressure)
end

viscosity = { Pr = 0.71, mu0 = mu0, T0 = ref_temp, n = 1, type = 'PowerLaw' }

-- === SETTINGS ===

verbosity = 1
cfl = 0.5

time_integration_scheme = 'euler'

refine_boxes = {
  { { -8, -8 }, { 8, 8 }, type = 'force_refine' },
}

-- === DEFINE STATES ===

states = {

  air = {
    type = 'hydro',
    gas = {
      type = 'thermally_perfect',
      mass = 1.0,
      charge = 0.0,
      gamma = 5 / 3,
    },
    reconstruction = 'O6',
    flux = 'HLLC',
    viscosity = viscosity,
    value = {
      rho = 1,
      x_vel = flow_u,
      y_vel = flow_v,
      p = flow_p,
      alpha = 0,
    },
    refinement = { name = 'hydro_gradient', rho = 0.1 },
  },
}

-- === ACTIONS ===

actions = {

  hydro_fluxes = {
    type = 'CTU',
    corner_transport = true,
    states = { 'air' },
  },
}
