time_integration_scheme = 'RK2'

verbosity = 2
cfl = 0.5

-- force_dt = 0.1

-- === DEFINE PROBLEM ===

kB = 1.38064852e-23

ref_temp = 4000.0 -- K
ref_prs = 1.0e5 -- Pa

ref_length = 1 -- m
ref_mass = 2.3258671e-26 -- kg
ref_vel = math.sqrt(kB * ref_temp / ref_mass) -- m/s
ref_density = ref_prs / (ref_vel ^ 2) -- kg/m^3

Larmor = 1.0
Debye = 1.0

-- === DEFINE STATES ===

states = {

  N = {
    type = 'hydro',
    gas = {
      type = 'eilmer',
      gas_model = 'nitrogen-2sp.lua',
      names = 'N',
      charge = 0,
    },
    reconstruction = 'centre',
    flux = 'HLLE',
    value = {
      p = 0.095,
      T = 1.0,
    },
  },

  N2 = {
    type = 'hydro',
    gas = {
      type = 'eilmer',
      gas_model = 'nitrogen-2sp.lua',
      names = 'N2',
      charge = 0,
    },
    reconstruction = 'centre',
    flux = 'HLLE',
    value = {
      p = 0.905,
      T = 1.0,
    },
  },

  combined = {
    type = 'hydro',
    gas = {
      type = 'eilmer',
      gas_model = 'nitrogen-2sp.lua',
      names = { 'N', 'N2' },
      charge = { 0, 0 },
    },
    reconstruction = 'centre',
    flux = 'HLLE',
    value = {
      p = 1.0,
      T = 1.0,
      alpha = { 0.05 },
    },
  },
}

actions = {

  gas = {
    type = 'gas_kinetics',
    states = { 'N', 'N2' },
    gas_model = 'nitrogen-2sp.lua',
    chemistry_update = 'chem.lua',
  },

  gas2 = {
    type = 'gas_kinetics',
    states = { 'combined' },
    gas_model = 'nitrogen-2sp.lua',
    chemistry_update = 'chem.lua',
  },
}
