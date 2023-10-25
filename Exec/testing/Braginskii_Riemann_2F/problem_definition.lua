verbosity = 1
cfl = 0.25

-- === DEFINE PROBLEM ===
mu_0_dim = 1.25663706e-6 
ep_0_dim = 8.85418782e-12
kb = 1.38064852e-23
ref_lightspeed = 299792458.0; --dimensional speed of light 
q_ref = 1.60217662e-19 -- Coulombs
m_ref = 1.6726219000e-27 -- kg 

-- REFERENCE VALUES 
lightspeed = 2000; lightspeed_nd = lightspeed
ref_length = 1e-8 
ref_mass = m_ref
ref_density = 1e31*m_ref
beta = 1
mass_ratio = 100 --1836.0

ref_velocity = ref_lightspeed/lightspeed;
ref_T = ref_mass*ref_velocity*ref_velocity/kb

n_ref = ref_density/ref_mass
B_ref = math.sqrt(2*mu_0_dim*n_ref*ref_mass*ref_velocity*ref_velocity/beta)
ref_omega_c = q_ref*B_ref/ref_mass; 
ref_omega_p = math.sqrt(n_ref*q_ref*q_ref/ref_mass/ep_0_dim); 

t_i =  12*(math.pi)^(3./2.)*ep_0_dim*ep_0_dim*math.sqrt(ref_mass)*(kb*ref_T)^(3./2.)/(10*(q_ref)^(4)*ref_density/ref_mass)
t_e = 6*math.sqrt(2)*((math.pi)^(3./2.))*ep_0_dim*ep_0_dim*math.sqrt(ref_mass/mass_ratio)*((kb*ref_T)^(3./2.))/(10*(q_ref)^(4)*n_ref)

ref_nu_p = 1/t_i
ref_nu_e = 1/t_e
print("t_i:\t", t_i, "\nt_e:\t",  t_e)

ref_larmor_dim = ref_velocity/ref_omega_c;
ref_skin_dim = ref_mass/(q_ref*math.sqrt(mu_0_dim*ref_density)) 
ref_skin_nd = ref_skin_dim/ref_length;
print('\nNon dimensional ion skin depth:\t', ref_skin_nd)
print('Dimensional ion skin depth:\t', ref_skin_dim)
Larmor = ref_larmor_dim/ref_length --======================================important 
Debye = ref_skin_nd/lightspeed_nd --=======================================important 

ref_time = ref_length/ref_velocity
print('omega_c_tau\t', ref_omega_c*ref_time, '\nomega_p_tau\t', ref_omega_p*ref_time, '\nnu_p_tau\t', ref_nu_p*ref_time)

betaBoi = 2*(Larmor/ref_skin_nd)^2 -- magnetic interaction parameter
print('\nbeta_0', betaBoi)


--[[
lightspeed = 2000 -- 50 --100.0
Larmor = 5.09 -- 1/30 --1/300
Debye = 3.6e-3 --0.1 --Larmor/lightspeed --Larmor/100
skinyBoi = Debye*lightspeed --skin depth
betaBoi = 2*(Larmor/skinyBoi)^2 -- magnetic interaction parameter


mass_ratio = 100 --1836.0
mu_0_dim = 1.25663706e-6 ; -- m kg s-2 A-2

ep0_dim = 8.85418782e-12
kb = 1.38064852e-23

n_ref = 1e31-- 1e20 * 1e-6 ; --cm3 to m3 --1e7;
u_ref = 299792458.0/lightspeed;
ref_mass = 1.6726219000e-27-- 1.6726219000e-27
B_ref = math.sqrt(2*mu_0_dim*n_ref*ref_mass*u_ref*u_ref/betaBoi)
q_ref = 1.60217662e-19 -- Coulombs

ref_length =  1e-8 -- 1e11 -- m
ref_mass = B_ref^2*betaBoi/(2*mu_0_dim*n_ref*u_ref^2) -- kg
ref_density = ref_mass*n_ref -- kg/m^3
i_skin = ref_mass/(q_ref*math.sqrt(mu_0_dim*ref_density))

time_ref = ref_length/u_ref

p_ref = ref_density*u_ref^2
T_ref = p_ref/n_ref

print('Ion skin depth:\t', i_skin)

print('Non-dimensional ion skin depth:\t', i_skin/ref_length)

print('Dimensional Collision time scales, t_e, t_i')
t_e = 6*math.sqrt(2)*(math.pi)^(3./2.)*ep0_dim*ep0_dim*math.sqrt(ref_mass/mass_ratio)*(kb*T_ref)^(3./2.)/(10*(q_ref)^(4)*n_ref)
t_i = 12*(math.pi)^(3./2.)*ep0_dim*ep0_dim*math.sqrt(ref_mass)*(kb*T_ref)^(3./2.)/(10*(q_ref)^(4)*n_ref)
print(t_e)
print(t_i)

print('Non-dimensional collision time scales, t_e, t_i')
print(t_e/time_ref, t_i/time_ref)

--]]
print('Knudsen number in this regiem')
Kn_dim = t_i*ref_velocity/ref_length
print(Kn_dim)

if (Kn_dim > 10e-5) and (Kn_dim<10e-2) then
    print('Braginskii model appropriate')
end 

eta0 = (ref_mass/q_ref/ref_density)*(ref_mass/mass_ratio/q_ref/t_e) -- background resistivity
v_a = B_ref/math.sqrt(mu_0_dim*n_ref*ref_mass) -- alfven velocity
Re_m = ref_density* v_a * ref_length / eta0

interface_x = 0.0

rho0 = 1.0
p0 = 0.5
u0 = 0.0
v0 = 0.0
w0 = 0.0

Bx0 = 0.75
By0 = 1.0
Bz0 = 0.0

rho1 = 0.125
p1 = 0.05
u1 = 0.0
v1 = 0.0
w1 = 0.0

Bx1 = 0.75
By1 = -1
Bz1 = 0.0

function step(A, B, x)
    if (x <= interface_x) then
        return A
    else
        return B
    end
end

-- === DEFINE STATES ===

states = {

    ion = {
        type='hydro',
        gas = {
          type='thermally_perfect', 
          mass=1.0,  
          charge= 1.0, 
          gamma=5/3, 
        }, 
        reconstruction='vanLeer',
        flux='HLLE',
        refinement={name='hydro_gradient', rho=0.1},
        value = {
            rho   = function(dat) return step(rho0, rho1, dat['x']) end,
            x_vel = function(dat) return step(u0,   u1,   dat['x']) end,
            y_vel = function(dat) return step(v0,   v1,   dat['x']) end,
            z_vel = function(dat) return step(w0,   w1,   dat['x']) end,
            p =   function(dat) return step(p0,   p1,   dat['x']) end,
        },
    },

    electron = {
        type='hydro',
        gas = {
          type='thermally_perfect', 
          mass=1.0/mass_ratio,  
          charge= -1.0, 
          gamma=5/3, 
        }, 
        reconstruction='vanLeer',
        flux='HLLE',
        refinement={name='hydro_gradient', rho=0.1},
        value = {
            rho   = function(dat) return step(rho0/mass_ratio, rho1/mass_ratio, dat['x']) end,
            x_vel = function(dat) return step(u0,   u1,   dat['x']) end,
            y_vel = function(dat) return step(v0,   v1,   dat['x']) end,
            z_vel = function(dat) return step(w0,   w1,   dat['x']) end,
            p =   function(dat) return step(p0,   p1,   dat['x']) end,
        },
    },

    field = {
        type='field',
        reconstruction='O6',
        flux='RankineHugoniot',
        value = {
            x_B  = function(dat) return step(Bx0, Bx1, dat['x']) end,
            y_B  = function(dat) return step(By0, By1, dat['x']) end,
            z_B  = function(dat) return step(Bz0, Bz1, dat['x']) end,
        },
        project_divergence = 1,
    },
}

actions = {
    em_fluxes = {
        type = 'CTU',
        corner_transport=true,
        states = {'field'},
    },

    braginskii = {
        -- this handles inviscid and viscous fluxes as well as inter-species collisions
        type = 'BraginskiiCTU',
        corner_transport=true,
        DebyeReference=10., 
        LarmorReference=10., 
        srin_switch = false,
        anisotropic = false,
        cfl=1.0,
        --force_ion_viscosity = 1e-3,
        --force_electron_viscosity = 1e-5,
        do_inter_species=false, 
        do_intra_species=true, 
        time_refinement_factor = 10,
        max_time_refinement_levels = 100,
        states = {ion='ion', electron='electron', field='field'},
    },

    plasma={
        type='plasma5',
        solver = 'explicit',
        states = {'ion', 'electron', 'field',},
     },

    divergence_cleaning = {
      type='elliptic',
      projection=1,
      state = 'field',
    },
}

