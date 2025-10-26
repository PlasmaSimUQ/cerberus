--########## Cerberus Numerical Settings ##########--
verbosity = 2
cfl = 0.3
time_integration_scheme = 'strang'

--########## Problem conditions ##########--
L_y = 4 * math.pi
L_x = 8 * math.pi

B_0 = 1.0
zeta_0 = 0.1
n_0 = 1.0
mass_ion = 1.0
mass_ele = 1. / 25.
mass_ratio = mass_ion / mass_ele
n_inf = 0.2 * n_0
kappa = 0.5

gamma_ele = 5. / 3.
gamma_ion = gamma_ele

q_ion = 1.0
q_ele = -1.0

mu_0 = 1.0 -- normalised mu_0
k_B = 1.0 -- normalised Boltzmann contsant

v_A = B_0 / math.sqrt(mu_0 * n_0 * mass_ion)

c_nd = 10 * v_A

eps_0 = 1 / ((c_nd * c_nd * mu_0) ^ 2)

T_ei = 0.2 --T_E/T_i

print('L_x = ', L_x, 'L_y = ', L_y, 'gamma_ele = ', gamma_ele)
print(
  'B_0: ',
  B_0,
  'zeta_0: ',
  zeta_0,
  'n_0: ',
  n_0,
  'mass_ion: ',
  mass_ion,
  'mass_ele: ',
  mass_ele,
  'mass_ratio: ',
  mass_ratio,
  'n_inf: ',
  n_inf,
  'kappa: ',
  kappa,
  'gamma_ele: ',
  gamma_ele,
  'gamma_ion: ',
  gamma_ion,
  'q_ion: ',
  q_ion,
  'q_ele: ',
  q_ele,
  'mu_0: ',
  mu_0,
  'k_B: ',
  k_B,
  'v_A: ',
  v_A,
  'c_nd: ',
  c_nd,
  'eps_0: ',
  eps_0,
  'T_ei: ',
  T_ei
)

--########## Non-dimensionalisation ##########--
lightspeed = c_nd

--[[
ref_lightspeed = 299792458.0
mu_0_dim = 1.25663706e-6
ep_0_dim = 8.85418782e-12
kb_dim = 1.38064852e-23
q_dim = 1.60217662e-19 -- Coulombs

--skin_depth = 1.0e+2
--lightspeed = 2000.0
beta = 1e-3 -- p_ref * 2 * mu0/ B_ref

ref_mass = 1.6726219000e-27
ref_velocity = ref_lightspeed / lightspeed
ref_density = ref_mass * 1e15

ref_T = ref_mass * ref_velocity * ref_velocity / kb_dim

n_ref = ref_density / ref_mass

ref_length = n_ref ^ (1. / 3.)

B_ref = math.sqrt(2 * mu_0_dim * n_ref * ref_mass * ref_velocity * ref_velocity / beta)
ref_omega_c = q_dim * B_ref / ref_mass
ref_omega_p = math.sqrt(n_ref * q_dim * q_dim / ref_mass / ep_0_dim)

tau_i = 12
  * math.pi ^ (3. / 2.)
  * ep_0_dim
  * ep_0_dim
  * math.sqrt(ref_mass)
  * (kb_dim * ref_T) ^ (3. / 2.)
  / (10 * q_dim ^ 4 * ref_density / ref_mass)

tau_e = 6
  * math.sqrt(2)
  * (math.pi ^ (3. / 2.))
  * ep_0_dim
  * ep_0_dim
  * math.sqrt(ref_mass / mass_ratio)
  * ((kb_dim * ref_T) ^ (3. / 2.))
  / (10 * q_dim ^ 4 * n_ref)

ref_nu_p = 1 / tau_i
ref_nu_e = 1 / tau_e
print('tau_i:\t', tau_i, '\ntau_e:\t', tau_e)

ref_larmor_dim = ref_velocity / ref_omega_c
ref_skin_dim = ref_mass / (q_dim * math.sqrt(mu_0_dim * ref_density))
ref_skin_nd = ref_skin_dim / ref_length
print('\nNon dimensional ion skin depth:\t', ref_skin_nd)
print('Dimensional ion skin depth:\t', ref_skin_dim)

Larmor = ref_larmor_dim / ref_length --======================================important
Debye = ref_skin_nd / c_nd --=======================================important

ref_time = ref_length / ref_velocity
print(
  'omega_c_tau\t',
  ref_omega_c * ref_time,
  '\nomega_p_tau\t',
  ref_omega_p * ref_time,
  '\nnu_p_tau\t',
  ref_nu_p * ref_time
)

betaBoi = 2 * (Larmor / ref_skin_nd) ^ 2 -- magnetic interaction parameter
print('\nbeta_0', betaBoi)

print('Knudsen number in this regiem')
Kn_dim = tau_i * ref_velocity / ref_length
print(Kn_dim)

if (Kn_dim > 10e-5) and (Kn_dim < 10e-2) then
  print('Braginskii model appropriate')
end

eta0 = (ref_mass / q_dim / ref_density) * (ref_mass / mass_ratio / q_dim / tau_e) -- background resistivity
v_a = B_ref / math.sqrt(mu_0_dim * n_ref * ref_mass) -- alfven velocity
Re_m = ref_density * v_a * ref_length / eta0

--Sutherland = { Pr = 0.72, mu0 = mu_0, T0 = 273, S = 110.4, type = 'Sutherland' }
--viscosity = { Pr = 1.0, mu0 = mu_0, type = 'UserDefined' }
--]]
-- === DEFINE PROBLEM ===
-- Auxiliary functions --

function tanh(x)
  if x == 0 then
    return 0.0
  end
  local neg = false
  if x < 0 then
    x = -x
    neg = true
  end
  if x < 0.54930614433405 then
    local y = x * x
    x = x
      + x
        * y
        * ((-0.96437492777225469787e0 * y + -0.99225929672236083313e2) * y + -0.16134119023996228053e4)
        / (((0.10000000000000000000e1 * y + 0.11274474380534949335e3) * y + 0.22337720718962312926e4) * y + 0.48402357071988688686e4)
  else
    x = math.exp(x)
    x = 1.0 - 2.0 / (x * x + 1.0)
  end
  if neg then
    x = -x
  end
  return x
end

function sech(x)
  return 2 / (math.exp(x) + math.exp(-x))
end

--########## State Functions ##########--
--number density
function number_density(dat)
  x = dat['x']
  y = dat['y']

  return n_0 * (sech((y - L_y / 2) / kappa)) ^ 2 + n_inf
end

--density
function rho_ion(dat)
  return mass_ion * number_density(dat)
end

function rho_ele(dat)
  return mass_ele * number_density(dat)
end

--pressure
function p_ion(dat)
  return number_density(dat) * (B_0 ^ 2) / (2 * mu_0 * k_B * n_0 * (1 + T_ei))
end

function p_ele(dat)
  return T_ei * p_ion(dat)
end

--velocity
function current_density_z(dat)
  return -B_0 / (mu_0 * kappa) * sech((dat['y'] - L_y / 2) / kappa) ^ 2
end

function w_ele(dat)
  return current_density_z(dat) / (number_density(dat) * q_ele)
end

--fluid energy density
function e_ion(dat)
  return (p_ion(dat) / rho_ion(dat)) * (1 / (gamma_ion - 1))
end

function e_e(dat)
  return w_ele(dat) / 2 + (p_ele(dat) / rho_ele(dat)) * (1 / (gamma_ele - 1))
end

--magnetic field
-- (0, 0, ez) cross (dx, dy, dz)phi = i( dy phi ) - j ( dx phi)  + k ( 0 )
function B1_x(x, y)
  return -zeta_0
    * (math.pi / L_y)
    * math.cos((2 * math.pi * (x - L_x / 2)) / L_x)
    * math.sin((math.pi * (y - L_y / 2)) / L_y)
end

function B1_y(x, y)
  return -zeta_0
    * (2 * math.pi / L_x)
    * math.sin(2 * math.pi * (x - L_x / 2) / L_x)
    * math.cos(math.pi * (y - L_y / 2) / L_y)
end

function B_x(dat)
  return B_0 * tanh((dat['y'] - L_y / 2) / kappa) - B1_x(dat['x'], dat['y'])
end

function B_y(dat)
  return -B1_y(dat['x'], dat['y'])
end

-- test
dat = { ['x'] = 1.0, ['y'] = 1.0 }

-- === DEFINE STATES ===

states = {
  ions = {
    type = 'hydro',
    gas = {
      type = 'thermally_perfect',
      mass = 1.0,
      charge = 1.0,
      gamma = 1.667,
    },
    reconstruction = 'vanLeer',
    flux = 'HLLC',
    --viscosity = viscosity,
    --refinement = { name = 'hydro_gradient', rho = 0.1 },
    value = {
      rho = rho_ion,
      x_vel = 0,
      y_vel = 0,
      z_vel = 0,
      p = p_ion,
    },
  },

  electrons = {
    type = 'hydro',
    gas = {
      type = 'thermally_perfect',
      mass = 0.04,
      charge = -1,
      gamma = 1.667,
    },
    reconstruction = 'vanLeer',
    flux = 'HLLC',
    --viscosity = viscosity,
    refinement = { name = 'hydro_gradient', rho = 0.1 },
    value = {
      rho = rho_ele,
      x_vel = 0,
      y_vel = 0,
      z_vel = w_ele,
      p = p_ele,
    },
  },

  field = {
    type = 'field',
    reconstruction = 'O6',
    flux = 'RankineHugoniot',
    value = {
      x_D = 0,
      y_D = 0,
      z_D = 0,
      x_B = B_x,
      y_B = B_y,
      z_B = 0,
    },
    refinement = { name = 'field_gradient', x_D = 0.2, y_D = 0.2, min_value = 1e-2 },
  },
}

-- === DEFINE ACTIONS ===

actions = {
  hydro_fluxes = {
    type = 'CTU',
    corner_transport = true,
    states = { 'ions', 'electrons' },
  },

  fluxes = {
    type = 'CTU',
    corner_transport = true,
    states = { 'field' },
  },

  plasma = {
    type = 'plasma5',
    --solver = 'implicit',
    solver = 'explicit',
    states = { 'ions', 'electrons', 'field' },
  },
  --[[
  divergence_cleaning = {
    type = 'elliptic',
    projection = 1,
    state = 'field',
  },
  --]]
}

-- === GEOMETRY ===
--Geometry for top and bottom PEC
-- options

refine_cutcells = true
merge_fraction = 0.5

wire_length = 3
wire_thickness = 0.2
wire_separation = 1.6

function make_rectangles(x, y, collection)
  local d, d1, dx, dy, dx_, dy_

  local coords = { x, y }

  for i, v in ipairs(collection) do
    dx = math.max(coords[1] - v[1][2], v[1][1] - coords[1])
    dy = math.max(coords[2] - v[2][2], v[2][1] - coords[2])

    dx_ = math.max(dx, 0.0)
    dy_ = math.max(dy, 0.0)

    d1 = math.sqrt(dx_ * dx_ + dy_ * dy_) + math.min(0.0, math.max(dx, dy))

    if i == 1 then
      d = d1
    else
      d = math.min(d, d1)
    end
  end

  return d
end

function wire_1(x, y)
  local rect = {
    { { 0., 25.2 }, { 12, 12.6 } }, -- {{xlo, xhi}, {ylo, yhi}}
  }
  return make_rectangles(x, y, rect)
end

function wire_2(x, y)
  local rect = {
    { { 0., 25.2 }, { -0.1, 0.5 } },
  }
  return make_rectangles(x, y, rect)
end

--Actually define the embedded boundaries for the code
embedded_boundaries = {
  solenoid_part_1 = {
    geom = wire_1,
    bcs = {
      field = { type = 'conductor' },
      ions = { type = 'slip_wall' },
      electrons = { type = 'slip_wall' },
    },
    boolean_operation = 'or',
    inside = 0,
  },

  solenoid_part_2 = {
    geom = wire_2,
    bcs = {
      field = { type = 'conductor' },
      ions = { type = 'slip_wall' },
      electrons = { type = 'slip_wall' },
    },
    boolean_operation = 'and',
    inside = 0,
  },
}

-- outputs --
plot = {
  variables = { 'all' },
  functions = {},
}
