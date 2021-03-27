
-- ======== PROBLEM ==========
ref_temp = 300.0

r_outer = 2.5
r_inner = 1
omega_outer = 0.1
omega_inner = 0.5

v_outer = r_outer*omega_outer
v_inner = r_inner*omega_inner

wall_velocity = math.max(v_outer, v_inner)

viscosity = {Pr=1.0, mu0=0.005, type='UserDefined'}


-- === SETTINGS ===

verbosity = 1
cfl = 0.5

do_face_sources = 0
do_CTU = 1

-- === DEFINE STATES ===

states = {

    air = {
        type='hydro',
        mass=1.0,
        charge= 0.0,
        gamma=1.4,
        reconstruction='MC', 
        flux='HLLC',
        viscosity=viscosity,
        value = {
            rho = 1,
            p =   1,
        },
        eb_divergence={
          type='merge',
          merge_threshold=0.5,
        }
    },
}


-- === GEOMETRY ===

refine_cutcells = 1


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

function outer(x,y)

    local circles = {
      {r_outer, 0.0, 0.0, 0.0},
    }

    return make_circles(x,y,circles)
end

function inner(x,y)

    local circles = {
      {r_inner, 0.0, 0.0, 0.0},
    }

    return -make_circles(x,y,circles)
end


embedded_boundaries = {

    inner = {
      geom=inner,
      bcs={air={type='no_slip_wall', v1=v_inner},
      },
      boolean_operation='or',
    },

    source = {
      geom=outer,
      bcs={air={type='no_slip_wall', v1=-v_outer},
      },
      boolean_operation='and',
    },
}
