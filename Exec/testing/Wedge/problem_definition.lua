
verbosity = 2
cfl = 0.5

-- === DEFINE PROBLEM ===

ref_density = 4.4e-2 --1.20 -- kg/m^-3
ref_mass = 4.65*10^-26  -- kg
ref_temp = 4415.0 -- K
kB = 1.38064852*10^-23
R = kB/ref_mass
ref_vel = math.sqrt(R*ref_temp) -- m/s
ref_length = 1.0 -- m

flow_velocity = 6360.0 -- m/s
flow_velocity = flow_velocity/ref_vel

Sutherland = {Pr=0.72, mu0=1.71e-5, T0=273, S=110.4, type='Sutherland', cfl=0.5}


p_ref = ref_density*R*ref_temp

p_0 = 5/p_ref
T_0 = 300/ref_temp
rho_0 = p_0/T_0


-- === DEFINE STATES ===

states = {

    fluid = {
        type='hydro',
        mass=1.0,  
        charge= 0.0, 
        gamma=1.4, 
        reconstruction='MC', 
        flux='HLLE/HLLC',
        shock_detector={name='pressure_jump_detector', threshold=0.1},
        refine_grad_threshold = {rho=0.5},
        -- viscosity=Sutherland,
        eb_divergence={
          type='merge',
          merge_threshold=0.5,
        },
        value = {
            rho = rho_0,
            p =   p_0,
        },
        bc = {
            x={
                lo={
                    fill_hydro_bc = 'outflow',
                    rho = 1.0,
                    p = 1.0,
                    x_vel = flow_velocity,
                }
            },
            y={
                lo={
                    fill_hydro_bc = 'symmetry',
                }
            },
        },
    },
}

-- === GEOMETRY ===

refine_cutcells = 1

r = 0.005 -- radius of leading edge
theta = math.rad(7.5) -- half angle of wedge
max_x = 1.0 -- how far back to go

p0 = {0,0}

p1 = {
    p0[1] - r*math.cos(math.pi/2 - theta),
    p0[2] + r*math.sin(math.pi/2 - theta),
}

p11 = {p1[1], -p1[2]}

p2 = {
    p0[1] + max_x*math.cos(theta) - p1[2],
    p0[2] + max_x*math.sin(theta) + p1[2],
}

p22 = {p2[1], -p2[2]}




spline = PolySpline.new()

spline:addLine({
  p2, p1, p11, p22, p2
  }
)


function lines(x, y)
  return spline:query(x, y)
end

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

function leading_edge(x,y)

    local circles = {
      {r, p0[1], p0[2]},
    }

    return make_circles(x,y,circles)
end


embedded_boundaries = {

  nose = {
    geom=leading_edge,
    bcs={fluid={type='slip_wall'},
    },
    boolean_operation='or',
    inside=0,
  },

  wedge = {
    geom=lines,
    bcs={fluid={type='slip_wall'},
    },
    inside=0,
    boolean_operation='and',
  },
}
