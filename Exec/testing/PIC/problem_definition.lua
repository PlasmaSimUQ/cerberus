-- === SETTINGS ===

verbosity = 1
cfl = 0.5

time_integration_scheme = 'symplectic'

skin_depth = 1.0
beta = 1.0
lightspeed = 1.0

-- === FUNCTIONS ===

function density(dat)

    local x = dat['x']
    local y = dat['y']
    local z = dat['z']
    
    local n
    local r = (x)^2 + (y)^2
    
    if (r<1^2) then
      n = 1.1e6
    else
      n = 0.0
    end
    
    return n
  end

-- === DEFINE STATES ===

states = {

    electron = {
        type='charged_particle',
        particles={
            type= 'distribution',
            distribution='Maxwell',
            mass = 0.1,
            charge = -0.01,
            nd = density,
            x_vel = 0.1,
            y_vel = 0.0,
            z_vel = 0.0,
            T = 1e-5,
            particles_per_cell = 20
          },
    },

    ion = {
        type='charged_particle',
        particles={
            type= 'distribution',
            distribution='Maxwell',
            mass = 1.0,
            charge = 0.01,
            nd = density,
            x_vel = -0.1,
            y_vel = 0.0,
            z_vel = 0.0,
            T = 1e-5,
            particles_per_cell = 20
          },
    },

    field = {
        type='field',
        reconstruction='O6',
        static = 1,
        boundary_conditions = {
            x={
                lo={
                    fill_D_bc='symmetry', fill_B_bc='symmetry',
                },
                hi={
                    fill_D_bc='inflow', x_D=0.0, y_D=0.0, z_D=0.0, 
                    fill_B_bc='inflow', x_B=0.0, y_B=0.0, z_B=0.0
                },
            },
            y={
                lo={
                    fill_D_bc='symmetry', fill_B_bc='symmetry',
                },
                hi={
                    fill_D_bc='inflow', x_D=0.0, y_D=0.0, z_D=0.0, 
                    fill_B_bc='inflow', x_B=0.0, y_B=0.0, z_B=0.0
                },
            },
          }
    },

}

actions = {
    plasma={
        type='symplectic',
        time_order = 2,
        allow_field_generation = 1,
        states = {'ion', 'electron', 'field'},
     },
}
