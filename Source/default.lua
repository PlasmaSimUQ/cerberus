/*

--Non-standard AMReX options:

-- archive level data within a checkpoint folder with the aid of the system tar utility
amr.check_archive = [0,1]

--*/R"(
-- ======== MISC ==========

-- physical definitions
ref_length = 1.0 -- m
ref_density = 1.0 -- kg/m^-3
ref_mass = 1.0 -- kg

ref_temp = -1 -- K
lightspeed = -1 -- non-dim

beta = 1 -- non-dim
skin_depth = 1 -- non-dim
Larmor = -1 -- non-dim
Debye = -1 -- non-dim

verbosity = 0
linear_solver_verbosity = 0

cfl = 0.0

-- refine around cutcells
refine_cutcells = 0

time_integration_scheme = 'strang'

force_dt = 0 -- force specific time step

-- list of boxes with {{{x_lo, y_lo, z_lo}, {x_hi, y_hi, z_hi}, type='?'}, ...}
-- type: 'force_refine' : refine to max level (default)
--       'no_refine' : do not refine
--       'only_refine' : only refine in regions marked as such
refine_boxes = {}

tile_size = {1024, 1024, 1024}

-- list of variables and functions to generate additional plot data
-- all variables used by functions must be listed in variables
-- if 'all' is specified as a variable then all primitive quantities will be plotted
-- e.g. variables = {'rho-electron', 'mass-electron'}
-- e.g. functions = {nd_electron = function f(dat) return dat['rho-electron']/dat['mass-electron'] end}
plot = {
    variables = {'all'},
    functions = {},
}

-- ======== STATES ==========

states = {}

actions = {}

--[[
embedded boundary definitions including boundary condition
embedded_boundaries = {
  name={
    geom=[function f(x,y,z) returning 0 for wall location, -ve for fluid, +ve for wall, or stl geometry file],
    bcs={[state-0]={
      type=[name of bc],
      other=...
      },
      ...
    }, -- boundary conditions for any state
    boolean_operation=[type of insertion], -- 'or', 'and'
    inside=1, -- is this geometry internal or external (allows for inversion of sign)
    -- etc
  }
}

The bcs list holds the definition of the boundary for the state that it is applied to.
The type of the bc must be given and is particular to the state which it is being applied.
All other entries in the list are defined by what is required by that particular type of bc.
NOTE: The order of the entries in the embedded_boundaries list is by sorted order
 of the names

EB bcs:

slip_wall : inviscid wall condition

no_slip_wall : viscous wall condition with specified wall velocity and temperature
  options = v1 : tangential velocity 1 (local y)
            v2 : tangential velocity 2 (local z)
            T  : wall temperature

conductor: perfectly conducting wall, can (optionally) specify surface B = (current density / relative permeability) and D = charge density
  options = B_t1 : magnetic field tangent 1 (local y)
            B_t2 : magnetic field tangent 2 (local z)
            D_n  : electrix flux normal (local x)

scalar_potential : sets electric scalar potential, can be function
  options = phi

vector_potential : sets magnetic vector potential, can be function
  options = A0, A1, A2 : x,y,z components
            align_with_boundary : set tangent with local wall-aligned coordinate system (x normal to wall)

surface_charge : sets the surface charge, can be function
  options = charge : surface charge density

surface_current : sets surface current density, can be function
  options = j1 : tangent 1 (local y)
            j2 : tangent 2 (local z)

Note: EB geometry must be defined over the whole domain, this can include ghost cells for
  boundary domains!
--]]

embedded_boundaries = {}

-- ======== LUA UTILITIES ==========

-- sorted pairs
function spairs(t, order)
    -- collect the keys
    local keys = {}
    for k in pairs(t) do keys[#keys+1] = k end

    -- if order function given, sort by it by passing the table and keys a, b,
    -- otherwise just sort the keys
    if order then
        table.sort(keys, function(a,b) return order(t, a, b) end)
    else
        table.sort(keys)
    end

    -- return the iterator function
    local i = 0
    return function()
        i = i + 1
        if keys[i] then
            return keys[i], t[keys[i]]
        end
    end
end

-- get the length of a table
function len(T)
  local count = 0
  for _ in pairs(T) do
    count = count + 1
  end
  return count
end

-- get a list of all the keys in a table
function get_sorted_keys(T)
  local keys = {}
  if len(T) < 1 then
    return keys
  end
  local n = 0
  for k,v in spairs(T) do
    n = n + 1
    keys[n] = k
  end
  return keys
end

-- given a table or a single value, expand it to be a table of n values
function expand(v, opt)
  local out = {}
  if (type(v) ~= "table") then
    out = {v,v or opt}
  else
    out = v
  end
  return out
end

-- check if something exists
function exists(v)
  if not v then
    return 0
  end
  return 1
end


-- checks if a variable is a function
function is_function(v)
  if (type(v) == "function") then
    return 1
  end
  return 0
end

-- checks if a variable is a table
function is_table(v)
  if (type(v) == "table") then
    return 1
  end
  return 0
end

-- get all integer values keys from a table
function get_ipairs(v)
  local out = {}
  for key, value in ipairs(v) do
    table.insert(out,value)
  end
  return out
end


-- gets the names of the arguments for a function
-- https://stackoverflow.com/a/29246308
function get_args(fun)
  local args = {}
  local hook = debug.gethook()

  local argHook = function( ... )
      local info = debug.getinfo(3)
      if 'pcall' ~= info.name then return end

      for i = 1, math.huge do
          local name, value = debug.getlocal(2, i)
          if '(*temporary)' == name then
              debug.sethook(hook)
              error('')
              return
          end
          table.insert(args,name)
      end
  end

  debug.sethook(argHook, "c")
  pcall(fun)

  return args
end

-- get number of arguments for a function
function get_num_args(fun)
    args = get_args(fun)
    return #args
end

-- currently empty function that is called after user code is defined
function preprocess()
end

--)"
