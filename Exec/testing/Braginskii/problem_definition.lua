-- ======== PROBLEM ==========

shock_x = -0.2
interface_x = 0.0
interface_amplitude = 0.1
interface_transition = 0.01

shock_mach = 2.0
density_L = 1.0
density_R = 3.0

mass_ion = 1.0
mass_electron = 0.01

gam = 5/3
pressure = 0.5
axis_beta = {0,0,0}

-- computations

p1 = pressure
p2 = pressure
p0 = p1*(1 + ((2*gam)/(gam+1))*(shock_mach^2 - 1))

rho1 = density_L
rho2 = density_R
rho0 = rho1/(1 - (2/(gam+1))*(1 - 1/shock_mach^2))

a0 = math.sqrt(gam*p0/rho0)
a1 = math.sqrt(gam*p1/rho1)
a2 = math.sqrt(gam*p2/rho2)

u1 = 0.0
u2 = 0.0
u0 = shock_mach*a1 - a0*math.sqrt(((gam-1)*shock_mach^2 + 2)/(2*gam*shock_mach^2 - (gam-1)))

--  magnetic field strength
B = {0,0,0}
p = 2*pressure -- use total pressure
for i, b in ipairs(axis_beta) do
    if b > 0 then
        B[i] = math.sqrt(p*beta/b)
    end
end


-- functions

local function islarge (x) return x > 2 ^ 28 end
local function issmall (x) return x < 2 ^ (-28) end

local function log1p (x) -- not very precise, but works well
  local u = 1 + x
  if u == 1 then return x end -- x < eps?
  return math.log(u) * x / (u - 1)
end

function tanh (x)
if x == 0 then return 0.0 end
local neg = false
if x < 0 then x = -x; neg = true end
if x < 0.54930614433405 then
    local y = x * x
    x = x + x * y *
        ((-0.96437492777225469787e0  * y +
        -0.99225929672236083313e2) * y +
        -0.16134119023996228053e4) /
        (((0.10000000000000000000e1  * y +
            0.11274474380534949335e3) * y +
            0.22337720718962312926e4) * y +
            0.48402357071988688686e4)
else
    x = math.exp(x)
    x = 1.0 - 2.0 / (x * x + 1.0)
end
if neg then x = -x end
return x
end


function atanh (x)
  y = math.abs(x)
  if y < .5 then
    if issmall(y) then return x end
    a = 2 * y
    a = .5 * log1p(a + a * y / (1 - y))
  else
    if y < 1 then
      a = .5 * log1p(2 * y / (1 - y))
    elseif y > 1 then
      return (x - x) / (x - x) -- nan
    else -- y == 1
      return x / 0 -- inf with sign
    end
  end
  return x < 0 and -a or a -- transfer sign
end

function RMI_interface_A(x, L, R)
    if x <= shock_x then
            return L
    else
            return R
    end
end

function RMI_interface_B(x, y, L, R)
    if math.abs(L - R)  < 1e-14 then
        return L
    end


    rr = x

    centre = interface_x + interface_amplitude*math.cos(2*math.pi*y)

    if (L < R) then
        at = (10.0*L - 9.0*R)/(10.0*(R-L))
    else
        at = (10.0*R - 9.0*L)/(10.0*(L-R))
    end

    slope = (2.0/interface_transition)*atanh(at)
    out = -((tanh(-slope*(rr - centre))-1.0)/2.0)*(L-R)+R

    return out
end

function RMI_interface_x(y)
    return interface_x + interface_amplitude*math.cos(2*math.pi*y)
end

function make_particles(N)

    p = {}

    local x, y, dy

    dy = 1/N

    for i=0,N-1,1 do
        y = i*dy
        table.insert(p, {RMI_interface_x(y), y})
    end

    return p

end

-- fluids

function number_density(dat)

    x = dat['x']
    y = dat['y']

    if (x <= shock_x) then
        n = RMI_interface_A(x, rho0, rho1)
    else
        n = RMI_interface_B(x, y, rho1, rho2)
    end

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
    if (x <= shock_x) then
        t = RMI_interface_A(x, 0, 1)
    else
        t = RMI_interface_B(x, y, 0, 1)
    end

    return t
end

function tracer_shock(dat)
    x = dat['x']
    y = dat['y']

        t = RMI_interface_A(x, 0, 1)

    return t
end

function pressure(dat)
    x = dat['x']
    y = dat['y']
    return RMI_interface_A(x, p0, p1)
end

function velocity_x(dat)
    x = dat['x']
    y = dat['y']
    return RMI_interface_A(x, u0, u1)
end

-- === SETTINGS ===

verbosity = 1
cfl = 0.3

time_integration_scheme = 'strang'

force_dt = 0

ref_mass = 1.6726219000e-27
ref_density = 1.6726219000e+04
ref_length = 1e-8
skin_depth =  7.2008467405e+00 
beta = 1.0
lightspeed = 2.0000000000e+03

-- debug delete 
--[[
yy = 0.5
for i=0,100,1 do
  xx = 2/100*i + -1
  print(xx, yy)
  dat = {x=xx, y=yy}
  print(ion_density(dat), velocity_x(dat), pressure(dat), tracer(dat),tracer_shock(dat))
  print(electron_density(dat), velocity_x(dat), pressure(dat), tracer(dat), tracer_shock(dat), '\n')
end 
--]]
--


-- === DEFINE STATES ===


states = {

    ion = {
        type='hydro',
        gas={
          type='thermally_perfect', 
          mass=mass_ion,
          charge=1,
          gamma=5/3,
        }, 
        reconstruction = 'minmod',
        refinement={name='hydro_gradient', rho=0.1},
        flux = 'HLLC',
        value = {
            rho   = ion_density,
            x_vel = velocity_x,
            p     = pressure,
            alpha = {tracer,tracer_shock},
        },
    },

    electron = {
        type='hydro',
        gas={
          type='thermally_perfect', 
          mass=mass_electron,
          charge=-1.0,
          gamma=5/3,
        }, 
        reconstruction = 'minmod',
        refinement={name='hydro_gradient', rho=0.1},
        flux = 'HLLC',
        value = {
            rho   = electron_density,
            x_vel = velocity_x,
            p     = pressure,
            alpha = {tracer,tracer_shock},
        },
    },


    field = {
        type='field',
        reconstruction='O6',
        flux = 'RankineHugoniot',
    },
  
    electron_tracer = {
        type='tracer',
        particles=make_particles(100)
    },

    ion_tracer = {
        type='tracer',
        particles=make_particles(100)
    }
}

actions = {
    --[[
    hydro_fluxes = {
        type = 'CTU',
        corner_transport=true,
        states = {'ion', 'electron'},
    },
    --]]

    braginskii = {
        -- this handles inviscid and viscous fluxes as well as inter-species collisions
        type = 'BraginskiiCTU',
        corner_transport=true,
        srin_switch = false,
        anisotropic = true,
        cfl=1.0,
        force_ion_viscosity = 1e-3,
        force_electron_viscosity = 1e-5,
        time_refinement_factor = 10,
        max_time_refinement_levels = 100,
        states = {ion='ion', electron='electron', field='field'},
    },

    plasma={
        type='plasma5',
        solver = 'explicit',
        states = {'ion', 'electron', 'field'},
     },
     divergence_cleaning = {
        type='elliptic',
        projection=1,
        state = 'field',
     },

     tracer_electron={
         type='hydro_tracer',
         particles='electron_tracer',
         fluid='electron'
     },

     tracer_ion={
         type='hydro_tracer',
         particles='ion_tracer',
         fluid='ion'
     }
}
