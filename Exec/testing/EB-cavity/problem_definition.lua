
-- ======== PROBLEM ==========

ref_density = 1.20 -- kg/m^-3
ref_mass = 4.65*10^-26  -- kg
ref_temp = 273.0 -- K
Sutherland = {Pr=0.72, mu0=1.71e-5, T0=273, S=110.4, type='Sutherland'}
ref_length = 1.0

gravity = -1.0
freestream_density = 1.0
freestream_temperature = 1.0

-- === SETTINGS ===

verbosity = 1
cfl = 0.5
time_integration_scheme = 'strang'

-- === GEOMETRY ===

refine_cutcells = 1

merge_fraction = 0.5

--svg_C = ReadSVG.new('drip_C.svg')
svg_C = ReadSVG.new('cerberus.svg', 0.01)

function svgC(x, y)
  return -svg_C:query(x, y,'/Layer 1')
end


embedded_boundaries = {

  C = {
    geom=svgC,
    bcs={air={type='no_slip_wall', T=10},
    },
    boolean_operation='or',
    inside=1,
  },

}

function inside(funcs, x, y)

  local r
  for i, f in ipairs(funcs) do
    if i == 1 then
      r = f(x,y)
    else
      r = math.max(r, f(x,y))
    end
  end

  if r <= 0 then
    return 1
  else
    return 0
  end
end

function make_particles(N)

  funcs = {}
  for k,v in pairs(embedded_boundaries) do
    table.insert(funcs, v['geom'])
  end

  pos = {}

  n = 0
  while n < N do
    x = math.random()
    y = math.random()

    if inside(funcs, x, y) > 0 then
      table.insert(pos,{x,y})
      n = n+1
    end
  end

  return pos
end

-- === DEFINE STATES ===

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

  local d, dd

  local coords = {x,y}

  for i, v in ipairs(collection) do
    for j=1,2 do
      dd = math.max(coords[j] - v[j][2], v[j][1] - coords[j])
      if (i == 1) then
          d = dd
      end
      d = math.max(d, dd)
    end
  end

  return -d

end

function pressure(dat)
  local h = dat['y']
  return math.exp(gravity*h/freestream_temperature)
end

function density(dat)
  local p = pressure(dat)
  return p/freestream_temperature
end

function blob(dat)

  local circles = {
    {0.1, 0.5, 0.8},
  }

  if make_circles(dat['x'],dat['y'],circles) < 0 then
    return 20.0
  else
    return density(dat)
  end
end

function tracer(dat)
  if blob(dat) == 20.0 then
    return 1.0
  else 
    return 0.0
  end
end



states = {

    air = {
        type='hydro',
        mass=1.0,
        charge= 0.0,
        gamma=1.4,
        reconstruction='constant', 
        flux='HLLC',
        viscosity=Sutherland,
        value = {
            rho =     blob,
            x_vel =   0,
            y_vel =   0,
            z_vel =   0,
            p =       pressure,
            alpha =   tracer,
        },
        refinement={name='hydro_gradient', rho=0.1},
        merge_threshold = 0.5,

        bc={
          x={
            lo={
              fill_hydro_bc='slipwall',
            },
            hi={
              fill_hydro_bc='slipwall',
            }
          },
          y={
            lo={
              fill_hydro_bc='slipwall',
            },
            hi={
              fill_hydro_bc='slipwall',
            }
          }
        }

    },

    tracer = {
      type='tracer',
      particles=make_particles(1000),
  },
}

-- === ACTIONS ===

actions = {

  hydro_fluxes = {
    type = 'CTU',
    corner_transport=true,
    states = {'air'},
},

  gravity = {
    type='acceleration',
    states={'air'},
    vector={0,gravity,0},
  },

  tracer={
    type='hydro_tracer',
    particles='tracer',
    fluid='air'
  },
  
}
