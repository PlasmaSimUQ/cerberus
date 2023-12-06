-- ========== SOLVER SETTINGS ========
verbosity = 0
cfl = 0.5
cfl_viscous = 0.5
--force_dt = 1e-5

--[[This does nothing in the new code, face sources may not even be in here?
do_CTU = 1 
do_face_sources = 1
--]]

-- === DEFINE PROBLEM ===
print('\n\n==========\nEelctron mass set to 0.001\n=========\n\n')
hydro_mass = { 0.01, 1.0 }
hydro_charge = { -1., 1. }
hydro_gamma = { 5 / 3.0, 5 / 3.0 }

--constants
mu_0_dim = 1.25663706e-6
ep_0_dim = 8.85418782e-12
kb = 1.38064852e-23

ref_lightspeed = 299792458.0 --dimensional speed of light 
q_ref = 1.60217662e-19 -- Coulombs
m_ref = 1.6726219000e-27 -- kg

-- REFERENCE VALUES
lightspeed = 2000
lightspeed_nd = lightspeed
ref_length = 1e-8
ref_mass = m_ref
ref_density = 1e31 * m_ref
beta = 1

ref_velocity = ref_lightspeed / lightspeed
ref_T = ref_mass * ref_velocity * ref_velocity / kb

n_ref = ref_density / ref_mass
B_ref = math.sqrt(2 * mu_0_dim * n_ref * ref_mass * ref_velocity * ref_velocity / beta)
ref_omega_c = q_ref * B_ref / ref_mass
ref_omega_p = math.sqrt(n_ref * q_ref * q_ref / ref_mass / ep_0_dim)
ref_nu_p = 12
  * math.pi ^ (3. / 2.)
  * ep_0_dim
  * ep_0_dim
  * math.sqrt(ref_mass)
  * (kb * ref_T) ^ (3. / 2.)
  / (10 * q_ref ^ 4 * ref_density / ref_mass)
ref_nu_p = 1 / ref_nu_p

ref_larmor_dim = ref_velocity / ref_omega_c
ref_skin_dim = ref_mass / (q_ref * math.sqrt(mu_0_dim * ref_density))
ref_skin_nd = ref_skin_dim / ref_length
print('\nNon dimensional ion skin depth:\t', ref_skin_nd)
print('Dimensional ion skin depth:\t', ref_skin_dim)
Larmor = ref_larmor_dim / ref_length --======================================important
Debye = ref_skin_nd / lightspeed_nd --=======================================important

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

-- NON-DIMENSIONAL CONSTANTS
v_plate = 1e-3
n0_nd = 1.0
Bx0_nd = 1.0
By0_nd = 0.0
Bz0_nd = 0
T0_nd = 1.0
print('\nBx0_nd = ', Bx0_nd)

--  magnetic field strength
--[[
B = {0,0,0}
axis_beta = {5e-1, 0, 0} --{5e4, 0, 0}
p = T0_nd*n0_nd * 2 -- use total pressure 

for i, b in ipairs(axis_beta) do
    if b > 0 then
        B[i] = math.sqrt(p*beta/b)
        print('i = ', i, '\tBi = ', B[i])
    end
end
--]]
-- ======== FLUID FUNCTIONS ==========

-- === DEFINE STATES ===
states = {
  ions = {
    type = 'hydro',
    gas = {
      type = 'thermally_perfect',
      mass = hydro_mass[2],
      charge = hydro_charge[2],
      gamma = hydro_gamma[2],
    },
    reconstruction = 'vanLeer',
    flux = 'HLLE',
    refinement = { name = 'hydro_gradient', z_vel = 0.000001, y_vel = 0.01 },
    value = {
      rho = n0_nd * hydro_mass[2],
      p = T0_nd * n0_nd,
    },
    bc = {
      x = {
        lo = { fill_hydro_bc = 'noslipwall', y_vel = -v_plate },
        hi = { fill_hydro_bc = 'noslipwall', y_vel = v_plate },
      },
    },
  },

  electrons = {
    type = 'hydro',
    gas = {
      type = 'thermally_perfect',
      mass = hydro_mass[1],
      charge = hydro_charge[1],
      gamma = hydro_gamma[1],
    },
    reconstruction = 'vanLeer',
    flux = 'HLLE',
    refinement = { name = 'hydro_gradient', z_vel = 0.0001, y_vel = 0.01 },
    value = {
      rho = n0_nd * hydro_mass[1],
      p = T0_nd * n0_nd,
    },
    bc = {
      x = {
        lo = { fill_hydro_bc = 'noslipwall', y_vel = -v_plate },
        hi = { fill_hydro_bc = 'noslipwall', y_vel = v_plate },
      },
    },
  },

  field = {
    type = 'field',
    reconstruction = 'O6',
    flux = 'HLLE',
    value = { x_B = Bx0_nd },
    dynamic = { 'x_B' },
    bc = {
      x = {
        lo = { fill_B_bc = 'asymmetry', fill_D_bc = 'asymmetry' },
        hi = { fill_B_bc = 'asymmetry', fill_D_bc = 'asymmetry' },
      },
    },
    --project_divergence=1,
  },
}

actions = {
  --[[
    hydro_fluxes = {
        type = 'CTU',
        corner_transport=true,
        states = {'ions', 'electrons'},
    },
    --]]

  em_fluxes = {
    type = 'CTU',
    corner_transport = true,
    states = { 'field' },
  },

  braginskii = {
    -- this handles inviscid and viscous fluxes as well as inter-species collisions
    type = 'BraginskiiCTU',
    corner_transport = true,
    hall_correction = false,
    DebyeReference = Debye,
    LarmorReference = Larmor,
    srin_switch = false,
    anisotropic = false,
    cfl = cfl_viscous,
    force_ion_viscosity = 0.010779311214565553,
    force_electron_viscosity = 0.0005819809305517557,
    do_inter_species = true,
    do_intra_species = true,
    time_refinement_factor = 10,
    max_time_refinement_levels = 100,
    states = { ion = 'ions', electron = 'electrons', field = 'field' },
  },

  plasma = {
    type = 'plasma5',
    solver = 'implicit',
    --solver = 'explicit',
    states = { 'ions', 'electrons', 'field' },
  },

  --[[
    divergence_cleaning = {
      type='elliptic',
      projection=1,
      state = 'field',
    },
    --]]
}

--source_collision_frequency_constraint = 0; //TODO not implemented

-- === PLOTTING ===
--[
plot = {
  variables = {
    'all',
    'x_D-field-dx',
    'y_D-field-dx',
    'z_D-field-dx',
    'x_B-field-dx',
    'y_B-field-dx',
    'z_B-field-dx',
    'p-electrons-dx',
    'p-ions-dx',
    'x_vel-ions-dx',
    'y_vel-ions-dx',
    'z_vel-ions-dx',
    'x_vel-electrons-dx',
    'y_vel-electrons-dx',
    'z_vel-electrons-dx',
    'T-ions-dx',
    'T-electrons-dx',
    'charge-electrons',
    'charge-ions',
    'mass-ions',
    'mass-electrons',
    'rho-electrons',
    'rho-ions',
  },
  functions = {},
}
--]
