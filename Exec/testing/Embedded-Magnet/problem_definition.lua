
-- ======== PROBLEM ==========


ref_temp = 273.0

Larmor = 0.1
Debye = 0.1

-- === SETTINGS ===

verbosity = 1
cfl = 0.25

do_face_sources = 0
do_CTU = 1

-- === DEFINE STATES ===

shock_x = -200

shock_mach = 2.0
density = 1.0

mass_ion = 1.0
mass_electron = 0.01

gam = 5/3
pressure = 0.5
axis_beta = {0,0,0}

-- computations

p1 = pressure
p0 = p1*(1 + ((2*gam)/(gam+1))*(shock_mach^2 - 1))

rho1 = density
rho0 = rho1/(1 - (2/(gam+1))*(1 - 1/shock_mach^2))

a0 = math.sqrt(gam*p0/rho0)
a1 = math.sqrt(gam*p1/rho1)

u1 = 0.0
u0 = shock_mach*a1 - a0*math.sqrt(((gam-1)*shock_mach^2 + 2)/(2*gam*shock_mach^2 - (gam-1)))


-- functions

function shock_interface(x, L, R)
    if x <= shock_x then
	    return L
    else
	    return R
    end
end


-- fluids

function number_density(dat)

    x = dat['x']
    y = dat['y']

    n = shock_interface(x, rho0, rho1)

    return n
end

function ion_density(dat)
    return mass_ion*number_density(dat)
end

function electron_density(dat)
    return mass_electron*number_density(dat)
end

function tracer(dat)
    x = dat['x']
    y = dat['y']
    t = shock_interface(x, 0, 1)
    
    return t
end

function pressure(dat)
    x = dat['x']
    y = dat['y']
    return shock_interface(x, p0, p1)
end

function velocity_x(dat)
    x = dat['x']
    y = dat['y']
    return shock_interface(x, u0, u1)
end





states = {

  electron = {
    type='hydro',
    mass=mass_electron, 
    charge=-1.0, 
    gamma=5/3, 
    reconstruction='minmod',  
    flux='HLLC',
    refine_grad_threshold = {rho=0.1},
    value = {
      rho   = electron_density,
      x_vel = velocity_x,
      p     = pressure,
      alpha = tracer,
    },
  },

  ion = {
    type='hydro',
    mass=1.0,  
    charge= 1.0, 
    gamma=5/3, 
    reconstruction='minmod', 
    flux='HLLC',
    refine_grad_threshold = {rho=0.1},
    value = {
      rho   = ion_density,
      x_vel = velocity_x,
      p     = pressure,
      alpha = tracer,
    },

  },

  field = {
    type='field',
    static=1,
    reconstruction='O6',
  },
}

-- === SOURCES ===

sources = {

  plasma={
      solver = 'implicit',
      options = {order=3},
      sources = {
          plasma={'ion', 'electron', 'field',
              type='Lorentz',
          },
      },
  },    
}


-- === GEOMETRY ===

-- options
refine_cutcells = 1
merge_fraction = 0.5

wire_length = 40
wire_thickness = 5
wire_separation = 5
shield_radius = 25

Az_magnitude = 1

function make_circles(x,y,collection)

  local d, dd

  for i, v in ipairs(collection) do
    dd = v[1]^2 - ((x-v[2])^2 + (y-v[3])^2)
    if (i == 1) then
        d = dd
    end
    d = math.max(d, dd)
  end

  return -d

end

-- collection = {{{x_lo, x_hi},{y_lo,y_hi}}, ... }
function make_rectangles(x,y,collection)

  local d, d1, dx, dy, dx_, dy_

  local coords = {x,y}

  for i, v in ipairs(collection) do
    dx = math.max(coords[1] - v[1][2], v[1][1] - coords[1])
    dy = math.max(coords[2] - v[2][2], v[2][1] - coords[2])

    dx_ = math.max(dx, 0.0)
    dy_ = math.max(dy, 0.0)

    d1 = math.sqrt(dx_*dx_ + dy_*dy_) + math.min(0.0, math.max(dx, dy))

    if (i == 1) then
      d = d1
    else 
      d = math.min(d, d1)
    end
  end

  return d

end

function shield_def(x,y)
  local circles = {
    {shield_radius, 0.0, 0.0},
  }
  return make_circles(x,y,circles)
end


function wire_1(x,y)
  local rect = {
    {{-wire_length/2, wire_length/2},{wire_separation, wire_separation+wire_thickness}},
  }
  return make_rectangles(x,y,rect)
end

function wire_2(x,y)
  local rect = {
    {{-wire_length/2, wire_length/2},{-wire_separation-wire_thickness, -wire_separation}},
  }
  return make_rectangles(x,y,rect)
end

embedded_boundaries = {

  shield = {
    geom=shield_def,
    bcs={
      ion={
        type='slip_wall'
      },
      electron={
        type='slip_wall'
      },
    },
    boolean_operation='or',
    inside=0,
  },

  solenoid_part_1 = {
    geom=wire_1,
    bcs={
      field={
        types={
          'scalar_potential', 'vector_potential'
        }, 
        phi=0.0, 
        A0=0.0,
        A1=0.0,
        A2=-Az_magnitude,
        align_with_boundary=false,
      },
    },
    boolean_operation='or',
    inside=0,
  },

  solenoid_part_2 = {
    geom=wire_2,
    bcs={
      field={
        types={
          'scalar_potential', 'vector_potential'
        }, 
        phi=0.0, 
        A0=0.0,
        A1=0.0,
        A2=Az_magnitude,
        align_with_boundary=false,
      },
    },
    boolean_operation='and',
    inside=0,
  }
}

-- === PLOTTING ===

function outside(dat)
  if ((dat['vfrac-electron'] == 0.0) or (dat['vfrac-ion'] == 0.0)) then
      return true
  else 
      return false
  end
end

function charge_density(dat)

  if (outside(dat)) then
      return 0.0
  end

  local cd_e = dat['charge-electron']*dat['rho-electron']/dat['mass-electron']
  local cd_i = dat['charge-ion']*dat['rho-ion']/dat['mass-ion']
  
  return (cd_e + cd_i)*(dat['Larmor']/dat['Debye']^2)
end

function div_D(dat)

  if (outside(dat)) then
      return 0.0
  end
  
  return dat['x_D-field-dx'] + dat['y_D-field-dy']
end

function err(dat)
  local cd = charge_density(dat)
  local div = div_D(dat)
  return cd - div
end

plot = {
  variables = {
      'all',
      'vfrac-electron',
      'vfrac-ion',
      'charge-electron',
      'charge-ion',
      'rho-electron',
      'rho-ion',
      'mass-electron',
      'mass-ion',
      'x_D-field-dx',
      'y_D-field-dy',
  },
  functions = {
      charge_density=charge_density,
      div_D=div_D,
      err_div_D=err,
  },
}

