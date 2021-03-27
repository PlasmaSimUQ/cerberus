/*

--Non-standard AMReX options:

-- archive level data within a checkpoint folder with the aid of the system tar utility
amr.check_archive = [0,1]

--*/R"(
-- ======== MISC ==========

zero_dimensional = 0 -- flag to allow transport

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

-- minimum volume fraction to allow before cell merging is carried out
merge_fraction = 0.5

-- if all density in a block is below effective_zero no calculations will occur for that block
effective_zero = 1e-14

do_CTU = 1 -- use corner-transport-upwinding

do_face_sources = 1 -- apply source terms to reconstructed face values

force_dt = 0 -- force specific time step

plot_shock_detector = 0 -- option to save out the shock detector

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

--[[
Define each state and how it is to be solved.

Note on functions:
  Any Lua function that is to be used internally by the c++ portion of the code
  is expected to accept a single list with named elements
  i.e. 'f(a)' where 'a' is something like [x=1, y=2, z=3, t=0.5, rho=1, ...]


Note on boundary conditions:
For boundary condition functions the input 'a' will carry whatever local
  primitive state the boundary has when the function was called and thus there
  will be elements in 'a' with names corresponding to the primitive state vector.
  This data will always be taken from the internal side of the boundary face and
  is done for both cell and face boundary conditions.

Note on dimensionality:
All fluid definition quantities and boundary conditions are expected to be in non-dimensional form.
All other inputs (viscosity, collisions, etc) are expected to be in dimensional form
  and will be non-dimensionalised on input according to the reference parameteres
  set above i.e. ref_length, ref_density, etc

====================================
Valid options for all states:

- reconstruction : method for computing face centred values, options are:
    * 'constant'
    * 'minmod'
    * 'vanLeer'
    * 'MC'
    * 'centre'
    * 'O6'

- reflux : do refluxing (this is a good idea)

- refine_max_lev : maximum level to refine to, -1 = no limit

- refine_grad_threshold : list of refinement thresholds in range [0,1] 0: refine everything, 1:refine nothing
    can choose from primitive or conserved state vectors, tagging from conserved is recommended as this does
    not require extra computation to convert the conserved variables to primitives
    'min_value' may be specified to give an approximate threshold below which gradients will not be computed
    e.g. refine_grad_threshold = {rho=0.1, p=0.2, min_value=1e-6}

- dynamic : a list naming entries in the 'value' list that should be called at each time-step to
    re-define the associated flow state
    e.g. dynamic={'ep', 'mu'}

- eb_divergence : a list naming the type of divergence calculation to be used to accommodate
    boundaries in the divergence calculation
    e.g. eb_divergence = {type='redistribute', strategy='volume'}
    available types & options:
      * type = 'redistribute'
        - strategy = 'volume', 'density', 'energy' or 'uniform'
        - reredistribution_threshold = real number > 0.0
      * type = 'merge'
        -  merge_threshold = real number in [0,1]

--]]

eb_divergence_default = {type='merge', merge_threshold=0.5}

--[[


====================================
Valid options for a hydro state are:

- type : type of hydrodynamic fluid, options are:
    * 'hydro'
    * 'hydro_2p'
    * 'hydro_BT'
    * 'mhd'

- mass : mass of the species, can be single value or tuple to correspond with tracer value
    e.g. mass = 1
    e.g. mass = {1,2}

- charge : charge of species, single value or tuple
    e.g. charge = 1
    e.g. charge = {1,2}

- gamma : ratio of specific heats, single value or tuple
    e.g. gamma = 5/3
    e.g. gamma = {5/3, 7/5}

- value : list of named values/functions that provide the initial flow state
    a number density may be specified instead of a density i.e. use 'nd' instead of 'rho'
    e.g. value = {rho=1, u=func_u, v=0, w=0, p=func_p, alpha=0}

- flux : method for computing face fluxes, options are:
    * 'HLLE'
    * 'HLLC'
    * 'HLLE/HLLC' (needs the shock_threshold parameter)
    * 'AUSMDV'
    * 'EFM'

- bc : multi-element list defining values for all domain edge boundary conditions:
    note: not everything needs to be defined, defaults to outflow fill_bc
    fill_bc options are : fill_hydro_bc
    variable names are for primitives
    with values of
    * 'interior'
    * 'outflow'  (default)
    * 'symmetry'
    * 'asymmetry'
    * 'slipwall'
    * 'noslipwall'
    e.g. inflow_bc={x={lo={fill_hydro_bc='outflow, rho=1, u=0, v=func_v, w=0, p=1, alpha=0},
                       hi={fill_hydro_bc='noslipwall', rho=2, ...}
                      }
                   }

- viscosity : multi-element list defining the viscosity model used for this fluid
    entries are in dimensional form and thus the reference quantities must be defined for
    reasonable results
    e.g. viscosity={Prandtl=0.72, mu0=1.716e-5, T0=273, n=2/3, type='PowerLaw'}           --> power law
    e.g. viscosity = {Prandtl=0.72, mu0=1.7894e-5, T0=273.11, S=110.563, type='Sutherland'} --> Sutherland

    optional 'cfl' variable scales the reported maximum wave speed due to viscosity

- pressure_relaxation : gives the relaxation rate for two-pressure model

- shock_threshold : level for which a shock is considered present [0 = none, inf = all]

- particles : name of ASCII file defining particles

- enforce_positivity : force Density and Eden > 0.0 when doing source terms [1 = apply, 0 = don't apply, default = 0]
                       currently only works for the CVODE source term solver and only if CVODE_apply_constraints=1.

======================================
Valid options for the mhd state are:


- type : type of fluid, options are:
    * 'mhd'

- mass : mass of the species, can be single value or tuple to correspond with tracer value

- gamma : ratio of specific heats, single value or tuple
    e.g. gamma = 5/3
    e.g. gamma = {5/3, 7/5}

- flux : method for computing face fluxes, options are:
    * 'mhdHlle'

- bc : multi-element list defining values for all domain edge boundary conditions:
    note: not everything needs to be defined, defaults to outflow fill_bc
    fill_bc options are : fill_hydro_bc, fill_psi_bc, and fill_B_bc
    variable names are for primitives
    with values of
    * 'interior'
    * 'outflow'  (default)
    * 'symmetry'
    * 'asymmetry'
    * 'slipwall'
    * 'noslipwall'

  e.g. bc={x={lo={fill_hydro_bc='4'slipwall', rho=1, u=0, v=0, w=0, p=1, alpha=0,
                 fill_B_bc='outflow', Bx=0, By=0, Bz=0,
                 fill_psi_bc='outflow', psi=0},
            }
         }

- shock_threshold : level for which a shock is considered present [0 = none, inf = all]

- div_transport : single value giving the speed of cleaning factor transport as a ratio with
    the maximum wave speed in the domain (0 for off)

- project_divergence : flag to indicate if the approximate projection method should be used for divergence cleaning
    if the field is linked to fluids these will be used to calculate the divergence source for the D field
   div_transport & damping_ratio should be set to zero if this option is enabled
    0 = off, 1 = on

- enforce_positivity : force Density and Eden > 0.0 when doing source terms [1 = apply, 0 = don't apply, default = 0]
                       currently only works for the CVODE source term solver and only if CVODE_apply_constraints=1.
======================================
Valid options for the field state are:

- type : 'field'

- static : 0 or 1

- value : list of named values/functions that provide the initial flow state
    e.g. value = {Bx=1, By=func_By, Bz=0, Dx=0, Dy=func_Dx, Dz=0,
                  mu=1, ep=1}

- flux : method for computing face fluxes, options are:
    * 'HLLE'
    * 'fieldRH'


- bc : multi-element list defining values for all domain edge boundary conditions:
    note: not everything needs to be defined, defaults to outflow fill_bc
    variable names are for primitives
    fill_bc options are fill_D_bc, fill_B_bc, fill_ep_bc, fill_mu_bc, fill_psi_bc,
    and fill_phi_bc, with values of
      * 'interior'
      * 'outflow'  (default)
      * 'symmetry'
      * 'asymmetry'

    e.g. bc={x={lo={fill_D_bc=2, Dx=0, Dy=0, Dz=0,
                   fill_B_bc=2, Bx=0, By=0, Bz=0,
                   fill_psi_bc=2, psi=0,
                   fill_phi_bc=2, phi=0,
                   fill_ep_bc=2, ep=1,
                   fill_mu_bc=2, mu=1,
              }
           }

- div_transport : single value giving the speed of cleaning factor transport as a ratio with
    the maximum wave speed in the domain (0 for off)

- project_divergence : flag to indicate if the approximate projection method should be used for divergence cleaning
    if the field is linked to fluids these will be used to calculate the divergence source for the D field
   div_transport & damping_ratio should be set to zero if this option is enabled
    0 = off, 1 = on

--]]

states = {}

--[[
Define the source terms to apply.
The intent is to define systems of ODEs that are solved separately in order that
different solvers may be used for different components of the source terms. In order
to define these systems we use a nested list defining ODE systems along with the
source terms that are to be solved by that system and the solver to be used.

The structure is as follows:

sources = {
  [system_name] = {
    solver = [type_of_solver],
    options = {[solver_options]=value, ...},
    sources = {
      [name_of_source]={
        [names_of_states],...,
        type=[type_of_source],
        source_options={[source_options], ...},
        value={cons_var=const/func, ...} -- only for UDF
      }
    }
  }
}

where anything in [] is something to be defined and can be so as follows:

[system_name] = any string of your choosing (must be a valid lua variable name), is used as a descriptor for output purposes

[type_of_solver] = how to solve the system of ODEs. Options are:
  - 'explicit'
  - 'implicit

[solver_options] = solver specific options, including (not all implemented):
  - abstol (absolute tolerance)
  - reltol (relative tolerance)
  - nsteps (maximum number of solver steps
  - refine_factor (if the solver fails recursively refine the domain by this factor)
  - refine_limit (how many recursion stesp to use)
  - jacobian (use the provided function for the jacobian, otherwise use the solvers inbuilt approach)
  - order (RK steps)

  - CVODE ONLY:
    - CVODE_method (default = 2):
      - 1-> use the Adams method
      - 2-> use BDF method.
    - CVODE_nls_fixed_point (default = 0):
      - 0-> use the newton non-linear solver
      - 1-> use the fixed-point non-linear solver
    - CVODE_apply_constraints (default = 0, apply positivity constraints to the states where enforce_positvity=1):
      - 1-> apply positivity constraints, should only be used if absolutely needed.
            see CVODE user guide, section 4.5.2, "advice on controlling unphysical negative values"
      - 0-> default value, don't apply constraints

[type_of_source] = the type of source, options are:
  - plasma5 : enables interaction between ions, electrons, and em fields according to the five-moment plasma model
  - damp_divergence : implements damping of the hyperbolic/parabolic divergence cleaning, can/should be used in conjunction with hydro_em
  - collisions : introduces collisions between fluid species (may require hydro_ccs to be defined for any neutral fluid also)
  - two_pressure : two-pressure relaxation
  - UDF : user defined, supply a constant or function for any of the conserved variables in a 'value' list (like for initial conditions)
  - acceleration : bulk acceleration defined by a constant or a function that has access to location and time, i.e. {type='acceleration', x=func(), y=1, z=0}
  - current : current source defined by a constant or a function that has access to location and time, i.e. {type='current', x=func(), y=1, z=0}

[names_of_states] = any of the available named states, defined in the 'states' list

[name_of_source] = any text string of your choosing, is used as a descriptor for output purposes

[source_options] = source specific options
  - collisions:
    - cross_sections = list of collision cross sections between states in dimensional form (m^2).
      matrix filled out when processed so as long as a pair is mentioned then it will be in both
      locations in the full matrix. Example:
      cross_sections = {ion={electron=1e-10, neutral=1e-10}, electron={neutral=1e-10}}
  - damp_divergence:
    - damping_ratio : single value giving the ratio of damping rate to cleaning transport (0 for off)

--]]
sources = {}

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
