import numpy as np
import pdb 
from scipy.optimize import fsolve
import scipy.integrate as integrate
from sympy import Symbol, solve, nsolve
import matplotlib.pyplot as plt


#===================================================================================================#
def referenceParametersTemperature(x_ref, rho_ref, m_ref, T_ref, skin_nd, beta):
  ref_length = x_ref  
  ref_mass = m_ref
  ref_n = rho_ref/ref_mass
  ref_density = rho_ref
  ref_temp = T_ref
  ref_velocity = np.sqrt(kB*T_ref/m_ref)
  c_nd = c_dim/ref_velocity;

  ref_B = (2*mu0_dim*ref_n*ref_mass*ref_velocity*ref_velocity/beta)**0.5
  ref_E = c_dim*ref_B;
  ref_pressure = ref_n*ref_mass*ref_velocity*ref_velocity # reference pressure (dnyamic)
  ref_time = ref_length/ref_velocity

  n0_ref = ref_n*ref_length**3; print(f"n0_ref:\t{n0_ref}")
  # check the Larmor and Debye and consistent with reference parameters 
  ref_omega_c = ref_q*ref_B/ref_mass;
  ref_omega_p = np.sqrt(ref_n*ref_q*ref_q/ref_mass/ep0_dim);
  ref_skin = ref_mass/(ref_q*np.sqrt(mu0_dim*ref_density)) # dimensional ref skin depth
  ref_skin_nd = ref_skin/ref_length; # nondiensional ref skin depth 
  ref_larmor_dim = ref_velocity/ref_omega_c;
  print(f'Dimensional ref skin depth:\t{ref_skin}\nNon dimensional ref skin depth:\t{ref_skin_nd}')  
  ref_Larmor_nd = ref_larmor_dim/ref_length
  ref_Debye_nd  = ref_skin_nd/c_nd
  print(f"\nref_length = {ref_length:.3E}\nref_density = {ref_density:.3E}\nref_mass = {ref_mass:.3E}\nref_temp = {ref_temp:.3E}\nref_n = {ref_n:.3E}\nref_nhat = {n0_ref:.3E}\nref_t = {ref_time:.3E}\nref_B = {ref_B:.3E}\n")

  return ref_length, ref_density, c_dim, c_nd, ref_velocity, \
         ref_pressure, ref_B, ref_time, ref_temp, ref_Debye_nd, ref_Larmor_nd, ref_n, n0_ref, beta ;

def getIonCoeffs(alpha_e, charge_e, mass_e, T_e, nd_e, alpha_i, charge_i, mass_i, T_i, nd_i, \
                 B_xyz, Debye, Larmor, n0_ref, p_lambda):
  #Magnetic field
  Bx = B_xyz[0]; By = B_xyz[1]; Bz = B_xyz[2];

  # absence of the boltzmann constant due to usage of nondimensional variables.
  # Note that here charge_i = e^4*Z^4
  t_collision_ion = (Debye)**(4)*n0_ref\
                    *(12*np.sqrt(mass_i)*(3.14159*T_i)**(3./2.))\
                    /(p_lambda * (charge_i)**(4) * nd_i);
  print("\ngetIon function")
  #print(f"\tcharge_i\t{charge_i*ref_q:.3E}\n\tB_mag\t{ref_B*np.sqrt(B_xyz[0]**2 + B_xyz[1]**2 + B_xyz[2]**2):.3E}\n\tmass\t{mass_i*ref_mass}\n\tLarmor_nd\t{Larmor}\n\tLarmor\t{Larmor*ref_length}")
  omega_ci = charge_i * np.sqrt(B_xyz[0]**2 + B_xyz[1]**2 + B_xyz[2]**2)/mass_i/Larmor;
  omega_p  = np.sqrt(nd_i*charge_i*charge_i/mass_i/Debye/Debye) ;

  x_coef = omega_ci*t_collision_ion;

  #print(f"\tNondimensional omega_ci:\t{omega_ci:.3E}\n\tDimensional omega_ci:\t{omega_ci/ref_time:.3E}\n\tNondimensional x_coef_i:\t{x_coef:.3E}")
  #print(f"\tDimensional (via nonondimensional values) omega_ci:\t{omega_ci/ref_time:.3E}")
  #TODO fix up coefficients here also with tabled depending atomic number

  delta_kappa = x_coef*x_coef*x_coef*x_coef + 2.700*x_coef*x_coef + 0.677;
  delta_eta   = x_coef*x_coef*x_coef*x_coef + 4.030*x_coef*x_coef + 2.330;
  delta_eta2  = 16*x_coef*x_coef*x_coef*x_coef + 4*4.030*x_coef*x_coef + 2.330;
  #print(f"\tdelta_eta\t{delta_eta}\n\tdelta_eta2\t{delta_eta}")
  print(f"\tni_nd:\t{nd_i}\n\tTi_nd:\t{T_i}\n\tt_i:\t{t_collision_ion}")
  print(f"\tni_dim:\t{nd_i*ref_density/ref_mass}\n\tTi_dim:\t{T_i*ref_temp:.3E}\n\tt_i:\t{t_collision_ion*ref_time:.3E}")
  eta0 = 0.96*nd_i*T_i*t_collision_ion ;#* n0_ref;
  eta2 = nd_i*T_i*t_collision_ion*(6./5.*x_coef*x_coef+2.23)/delta_eta;
  eta1 = nd_i*T_i*t_collision_ion*(6./5.*(2*x_coef)*(2*x_coef)+2.23)/delta_eta2;

  #print(f"\talternate eta2 = {6./5. *nd_i*T_i*t_collision_ion/x_coef**2*ref_density*ref_length*ref_velocity:.3E}");

  eta4 = nd_i*T_i*t_collision_ion*x_coef*(x_coef*x_coef + 2.38)/delta_eta;
  eta3 = nd_i*T_i*t_collision_ion*(2*x_coef)*((2*x_coef)*(2*x_coef) + 2.38)/delta_eta2;
  kappa1 = 3.906*nd_i*T_i*t_collision_ion/mass_i;
  kappa2 = (2.*x_coef*x_coef + 2.645)/delta_kappa*nd_i*T_i*t_collision_ion/mass_i;
  kappa3 = (5./2.*x_coef*x_coef + 4.65)*x_coef*nd_i*T_i*t_collision_ion/mass_i/delta_kappa;
  #print("alternate kappa2 = ", 2/x_coef/x_coef*nd_i*T_i*t_collision_ion/mass_i*n0_ref*ref_thermal_conductivity);
  return (eta0, eta1, eta2, eta3, eta4, kappa1, kappa2, kappa3, t_collision_ion) ; 

def getElectronCoeffs(alpha_e, charge_e, mass_e, T_e, nd_e, alpha_i, charge_i, mass_i, T_i, nd_i, \
                 B_xyz, Debye, Larmor, n0_ref, p_lambda): 
    print("\ngetElectronCoeffs function")
    Bx=B_xyz[0]; By=B_xyz[1]; Bz=B_xyz[2];

    t_collision_ele = (Debye)**(4)*n0_ref\
                      *(6*np.sqrt(2*mass_e)*(3.14159*T_e)**(3./2.)) /\
                      (p_lambda*(charge_e)**(4)*(charge_i/-charge_e)*nd_e);

    #print(f"charge_e\t{charge_e*ref_q:.3E}\nB_mag\t{ref_B*np.sqrt(B_xyz[0]**2 + B_xyz[1]**2 + B_xyz[2]**2):.3E}\nmass\t{mass_e*ref_mass}\nLarmor\t{Larmor*ref_length}")
    omega_ce = -charge_e * np.sqrt(B_xyz[0]**2 + B_xyz[1]**2 + B_xyz[2]**2)/mass_e/Larmor;
    omega_p  = np.sqrt(nd_e*charge_e*charge_e/mass_e/Debye/Debye) ;

    x_coef = omega_ce*t_collision_ele;

    #print(f"Nondimensional omega_ce:\t{omega_ce:.3E}\nNondimensional x_coef_e:\t{x_coef:.3E}")
    #print(f"Dimensional (via nonondimensional values) omega_ce:\t{omega_ce/ref_time:.3E}")
    # TODO fix up these table 2 page 251 BT
    delta_0=3.7703; delta_1=14.79;
    delta_kappa= x_coef*x_coef*x_coef*x_coef+delta_1*x_coef*x_coef + delta_0;
    delta_eta  = x_coef*x_coef*x_coef*x_coef+13.8*x_coef*x_coef + 11.6;
    delta_eta2 = 16*x_coef*x_coef*x_coef*x_coef+4*13.8*x_coef*x_coef + 11.6;

    print(f"\tne_nd:\t{nd_e}\n\tTe_nd:\t{T_e}\n\tt_e:\t{t_collision_ele}")
    print(f"\tne_dim:\t{nd_e*ref_density/ref_mass}\n\tTe_dim:\t{T_e*ref_temp:.3E}\n\tt_e:\t{t_collision_ele*ref_time:.3E}")

    eta0 = 0.733*nd_e *T_e * t_collision_ele;
    eta2 = nd_e *T_e*t_collision_ele*(2.05*x_coef*x_coef+8.5)/delta_eta;
    eta1 = nd_e *T_e*t_collision_ele*(2.05*(2*x_coef)*(2*x_coef)+8.5)/delta_eta2;
    eta4 = -nd_e*T_e*t_collision_ele*x_coef*(x_coef*x_coef+7.91)/delta_eta;
    eta3 = -nd_e*T_e*t_collision_ele*(2*x_coef)*((2*x_coef)*(2*x_coef)+7.91)/delta_eta2;

    BT_gamma_0=11.92/3.7703;BT_gamma_0_p=11.92;BT_gamma_1_p=4.664;BT_gamma_1_pp=5./2.;
    BT_gamma_0_pp=21.67;

    kappa1 = nd_e*T_e*t_collision_ele/mass_e*BT_gamma_0;

    print(BT_gamma_1_p, x_coef, BT_gamma_0_p, delta_kappa, nd_e, T_e, t_collision_ele, mass_e)
    kappa2 = (BT_gamma_1_p*x_coef*x_coef+BT_gamma_0_p)/\
              delta_kappa*nd_e*T_e*t_collision_ele/mass_e;
    kappa3=(BT_gamma_1_pp*x_coef*x_coef+BT_gamma_0_pp)*x_coef*nd_e*T_e*t_collision_ele\
           /mass_e/delta_kappa;
    #--- terms for the resitivity
    delta_coef = x_coef*x_coef*x_coef*x_coef+delta_1*x_coef*x_coef + delta_0; # TODO delta0 tables 

    a_0 = 0.5129; a_0_p =1.837; a_1_p =6.416; a_0_pp =0.7796; a_1_pp =1.704;
  
    alpha_0 = mass_e*nd_e/t_collision_ele*a_0;
    alpha_1 = mass_e*nd_e/t_collision_ele*( 1 - (a_1_p*x_coef*x_coef + a_0_p)/delta_coef);
    alpha_2 = mass_e*nd_e/t_collision_ele*x_coef*(a_1_pp*x_coef*x_coef+a_0_pp)/delta_coef;
    #--- beta terms for the thermal component of thermal heat flux of the electrons.
    b_0 = 0.711; b_0_pp = 3.053; b_0_p=2.681; b_1_p=5.101; b_1_pp=3./2.;
    beta1 = nd_e*b_0*T_e;
    beta2 = nd_e*(b_1_p*x_coef*x_coef+b_0_p)/delta_kappa*T_e;
    beta3 = nd_e*x_coef*(b_1_pp*x_coef*x_coef+b_0_pp)/delta_kappa*T_e;

    rhor_para = 1/(nd_e*nd_e*charge_e*charge_e/alpha_0)
    rhor_perp = 1/(nd_e*nd_e*charge_e*charge_e/alpha_1)
    rhor_chev = 1/(nd_e*nd_e*charge_e*charge_e/alpha_2)
    return (eta0, eta1, eta2, eta3, eta4, kappa1, kappa2, kappa3, beta1, beta2, beta3, rhor_para, rhor_perp, rhor_chev, t_collision_ele);


def coeffMatrix():
  A = np.zeros((4,4), dtype=complex)

  (eta0, eta1, eta2, eta3, eta4, kappa1, kappa2, kappa3, tau) =\
  getIonCoeffs(0., charge_e_nd, mass_e_nd, Te_nd, ne_nd, 
               0., charge_i_nd, mass_i_nd, Ti_nd, ni_nd,
               B_field, Debye, Larmor, n0_ref, p_lambda_i_cerb)

  eta_i_0 = eta0; tau_i = tau
  
  (eta0, eta1, eta2, eta3, eta4, kappa1, kappa2, kappa3, 
  beta1, beta2, beta3, rhor_para, rhor_perp, rhor_chev, tau) = \
  getElectronCoeffs(0., charge_e_nd, mass_e_nd, Te_nd, ne_nd, 
                    0., charge_i_nd, mass_i_nd, Ti_nd, ni_nd, 
                    B_field, Debye, Larmor, n0_ref, p_lambda_e_cerb)

  omega_ci = charge_i_nd * np.sqrt(B_field[0]**2 + B_field[1]**2 + B_field[2]**2)/mass_i_nd/Larmor;#nondim
  omega_ce =-charge_e_nd * np.sqrt(B_field[0]**2 + B_field[1]**2 + B_field[2]**2)/mass_e_nd/Larmor;#nondim
  omega_pi  = np.sqrt(ni_nd*charge_i_nd*charge_i_nd/mass_i_nd/Debye/Debye) ;#nondim
  
  if not 1/tau_i > 10*omega_ci: xxx
  if not 1/tau_i>= 1e2: xxx
  if not 1/tau_i<= 1e5: xxx
  if not omega_pi > 10*omega_ci: xxx
  if not 10*v_plate < min(np.sqrt(Te_nd/mass_e_nd), np.sqrt(Ti_nd/mass_i_nd)): 
    print(f"v_plate = {v_plate}")
    print(f"v_e = {np.sqrt(Te_nd/mass_e_nd)}\tv_i = {np.sqrt(Ti_nd/mass_i_nd)}") 
    xxx

  eta_e_0 = eta0; tau_e = tau
  a_para = mass_e_nd*ne_nd*0.5129 / tau_e 
  Bx =  B_field[0]
  gamma_ei = a_para / eta_e_0
  gamma_ie = a_para / eta_i_0
  lambda_e = np.sqrt(2/beta)*ne_nd*charge_e_nd/skin_nd /  eta_e_0 * Bx
  lambda_i = np.sqrt(2/beta)*ne_nd*charge_e_nd/skin_nd /  eta_i_0 * Bx
  
  A[0,0] = gamma_ei; A[0,1] = -lambda_e; A[0,2] = -gamma_ei; A[0,3] = 0     
  A[1,0] = lambda_e; A[1,1] =  gamma_ei; A[1,2] =  0       ; A[1,3] = -gamma_ei
  A[2,0] = gamma_ie; A[2,1] =  0       ; A[2,2] = -gamma_ie; A[3,3] = -lambda_i
  A[3,0] = 0       ; A[3,1] =  gamma_ie; A[3,2] =  lambda_i; A[3,2] = -gamma_ei
  return A

def printArray(AA):
  try:
    (nRows, nCols) = AA.shape
  except:
    print("\t(vector)")
    nRows = AA.shape[0]
    nCols = 1
  
  outputString = ""
  for i in range(nRows):
    rowString = "|\t"
    for j in range(nCols):
      if nCols == 1: rowString += f"{AA[i]:10.4E}\t"
      else: rowString += f"{AA[i,j]:10.4E}\t"
    outputString += rowString + "\t|\n"
  print(outputString)

def genSol(inputs, Sil, Sir, Di, Xl, Xr):
  #global Sil, Sir, Di, Xl, Xr;  
  #print(Sil, Sir, Di, Xl, Xr)
  C0 = inputs[0]+inputs[1]*1j
  C1 = inputs[2]+inputs[3]*1j

  #return np.array([Sil - C0*np.exp(np.sqrt(Di)*Xl) - C1*np.exp(-np.sqrt(Di)*Xl),
  #                 Sir - C0*np.exp(np.sqrt(Di)*Xl) - C1*np.exp(-np.sqrt(Di)*Xl)])

  val1 = Sil - C0*np.exp(np.sqrt(Di)*Xl) - C1*np.exp(-np.sqrt(Di)*Xl)
  val2 = Sir - C0*np.exp(np.sqrt(Di)*Xr) - C1*np.exp(-np.sqrt(Di)*Xr)

  return np.array([np.real(val1), np.imag(val1), np.real(val2), np.imag(val2)])

def sGenSol(inputs):
  C0 = inputs[0]
  C1 = inputs[1]
  x  = inputs[2]
  Di = inputs[3]
  return C0*np.exp(np.sqrt(Di)*x) + C1*np.exp(-np.sqrt(Di)*x)

def sFundmentalSol(Di, x):
  return ( np.exp(np.sqrt(Di)*x), np.exp(-np.sqrt(Di)*x) )

def dBzdx(ne, qe, uey, ni, qi, uiy, ds, beta0):
  return -1/ds*np.sqrt(beta0/2)*(ne*qe*uey + ni*qi*uiy)

def dBydx(ne, qe, uez, ni, qi, uiz, ds, beta0):
  return 1/ds*np.sqrt(beta0/2)*(ne*qe*uez + ni*qi*uiz)

#===================================================================================================#
if __name__ == '__main__':
  # Set parameters not varying 
  c_dim = 299792458.0 #m/s
  mu0_dim = 1.25663706e-6 
  ep0_dim = 8.85418782e-12
  kB = 1.38064852e-23
  eV = 1.602176634e-19 # Joules
  ref_mass = 1.6726219e-27 # kg 
  mass_i_dim = ref_mass; mass_e_dim = mass_i_dim/1836
  ref_q = 1.602176634e-19 # Coulombs

  #Option 16
  x_ref = 1.0e-8;  m_ref = 1.6726219000e-27;  rho_ref = 1e31 * m_ref # 1.6726219000e+03
  T_ref = 2.7220497703e+06;  beta = 1;  skin_nd =7.2008467405e+00; c_nd = 2.0e+03
  p_lambda_i_cerb = 10;  p_lambda_e_cerb = 10;
  #standard :
  n0_nd = 1.0; Bx0_nd = 0.; By0_nd = 0.0; Bz0_nd = 0; T0_nd = 0.5;

  #Hartman flow parameters 
  v_plate = 1e-3 ; # dimensionless
  Bx0_nd = 0.0044721359549996; # dimensionless 5e-3 #
  Xdomain = 1;
  Ydomain = 1;
  #ref parameters for SRMI sim 
  ref_length, ref_density, ref_lightspeed, lightspeed, ref_velocity, ref_pressure, ref_B, ref_time, ref_temp, Debye, Larmor, ref_n, n0_ref, beta = \
  referenceParametersTemperature(x_ref, rho_ref, m_ref, T_ref, skin_nd, beta)

  #=====Non-dimensional quantities 
  mass_e_nd = 1/100; mass_i_nd = 1.; #masses 
  Te_nd = T0_nd; Ti_nd = Te_nd; # temepratures
  charge_i_nd = 1; charge_e_nd = -1; #charges
  ne_nd = n0_nd; ni_nd = n0_nd; #number densities 
  B_field = np.array([Bx0_nd, By0_nd, Bz0_nd])
  B_mag_nd = (Bx0_nd**2 + By0_nd**2 + Bz0_nd**2)**0.5  

  #=====Dimensional quantities 
  Te_dim = ref_temp * Te_nd; Ti_dim = ref_temp*Ti_nd; 
  T_grad_dim = -5 *eV / kB
  ne_dim = ref_n*ne_nd; ni_dim = ref_n*ni_nd;

  A = coeffMatrix()
  print("Playing with different eigvalue functions")  
  
  Eigenvalues, Q = np.linalg.eig(A)
  Lambda = np.diag(Eigenvalues)
  A_back = Q @ Lambda @ np.linalg.inv(Q)

  #===============================Solve boundary value problem======================================#
  
  normaliseValues = False

  #boundary value problem 
  U_left = np.array([-v_plate, 0, -v_plate, 0])
  U_right= np.array([ v_plate, 0,  v_plate, 0])

  s_left = np.linalg.inv(Q)@U_left
  s_right= np.linalg.inv(Q)@U_right

  cArray = np.zeros((s_left.shape[0], 2), dtype=complex)
  global Sil, Sir, Di, Xl, Xr;
  Xl = -Xdomain/2.; Xr = Xdomain/2.
  #C0 = Symbol('C0')
  #C1 = Symbol('C1')

  guesses = []
  for i in range(1,2):
    guesses.append(np.array([10000/i*(-1)**i]*4))

  for g in guesses:
    print(f"\nGuess = {g}")
    for i in range(Eigenvalues.shape[0]):
      Sil = s_left[i]; Sir = s_right[i]; Di = Eigenvalues[i];
      #eq1 = Sil - C0*np.exp(np.sqrt(Di)*Xl) - C1*np.exp(-np.sqrt(Di)*Xl)
      #eq2 = Sir - C0*np.exp(np.sqrt(Di)*Xl) - C1*np.exp(-np.sqrt(Di)*Xl)
  
      #sSol = nsolve( (eq1, eq2), (C0, C1), (0, 0) )
      sSol = fsolve( genSol, g, (Sil, Sir, Di, Xl, Xr))#, full_output=True)

      #print(f"\nSolution:\t{sSol[0]}+{sSol[1]}j\t {sSol[2]}+{sSol[3]}j")
      cArray[i,0] = sSol[0]+sSol[1]*1j; cArray[i,1] = sSol[2]+sSol[3]*1j;
      #print(f"Stored value:\t{cArray[i,0]}\t{cArray[i,1]}")
    xs = np.linspace(Xl, Xr, 1000)
    
    sx = np.zeros((s_left.shape[0], xs.shape[0]), dtype=complex)
    sFundamental = np.zeros((s_left.shape[0], xs.shape[0], 4))

    for i in range(s_left.shape[0]):
      Di = Eigenvalues[i]
      sx[i,:] = sGenSol([cArray[i,0], cArray[i, 1], xs, Di])
      (val1, val2) = sFundmentalSol(Di, xs)
      sFundamental[i,:,0] = np.real(val1)
      sFundamental[i,:,1] = np.real(val2)
      sFundamental[i,:,2] = np.imag(val1) 
      sFundamental[i,:,3] = np.imag(val2)

    # make a figure to plot as we get the data 
    fig_u = plt.figure(figsize=(16,20), dpi=300); ax_u = [] # physical solution 
    fig_f = plt.figure(figsize=(16,20), dpi=300); ax_f = [] # fundamental solutions

    # titles for each plot 
    titles = [r"$u^e_y/V$", r"$u^e_z/V$", r"$u^i_y/V$", r"$u^i_z/V$"]

    ### plot of physical values 
    ux = np.zeros(sx.shape, dtype=complex)
    for j in range(xs.shape[0]):
      ux[:,j] = Q @ sx[:,j] #/abs(v_plate)
      if normaliseValues: ux[:,j] /= abs(v_plate)
    
    # polyfits     # Polynomial fits for simulation
    dataFits = []
    dataFitList = [ux[0,:], ux[1,:], ux[2,:], ux[3,:]]
    xList = [xs, xs, xs, xs]
    degList = [4]*4
    for i in range(len(dataFitList)):
      dataFits.append(np.polyfit(xList[i], np.real(dataFitList[i]), 3))
      print(dataFits[-1])

    for i in range(s_left.shape[0]):
      ax_u.append(fig_u.add_subplot(2,2,i+1))
      ax_u[i].plot(xs, np.real(ux[i,:]))
      ax_u[i].plot(xs, np.poly1d(dataFits[i])(xs), label='fit')
      ax_u[i].set_ylabel(titles[i])
      ax_u[i].legend()

      if i < 2: ax_u[i].set_xticks([])
      else: ax_u[i].set_xlabel('x')
    
    fig_u.tight_layout(pad=10)

    ### plot of physical values 
    titles = [r"$s^e_y$", r"$s^e_z$", r"$s^i_y$", r"$s^i_z$"]
    for i in range(s_left.shape[0]):
      ax_f.append(fig_f.add_subplot(2,2,i+1))
      ax_f[i].plot(xs, sFundamental[i,:,0], label='g1 Real')
      ax_f[i].plot(xs, sFundamental[i,:,1], label='g2 Real')
      ax_f[i].plot(xs, sFundamental[i,:,2], label='g1 Imag')
      ax_f[i].plot(xs, sFundamental[i,:,3], label='g2 Imag')
      ax_f[i].legend()
      ax_f[i].set_ylabel(titles[i])
      if i < 2: ax_f[i].set_xticks([])
      else: ax_f[i].set_xlabel('x')

    fig_f.tight_layout(pad=10)

    #================================================= Derived quantities ===============================================#
    Jy_vals = (ne_nd*charge_e_nd*ux[0,:] + ni_nd*charge_i_nd*ux[2,:])#/ni_nd/abs(v_plate)
    Jz_vals = (ne_nd*charge_e_nd*ux[1,:] + ni_nd*charge_i_nd*ux[3,:])#/ni_nd/abs(v_plate)
    if normaliseValues:
      Jy_vals /= (ni_nd*abs(v_plate))  ; Jz_vals /= (ni_nd*abs(v_plate));

    Jy_vals = np.real(Jz_vals) ; Jz_vals = np.real(Jz_vals)

    dBydx_vals = dBydx(ne_nd, charge_e_nd, np.real(ux[1,:]), 
                       ni_nd, charge_i_nd, np.real(ux[3,:]), skin_nd, beta)
    dBzdx_vals = dBzdx(ne_nd, charge_e_nd, np.real(ux[0,:]), 
                       ni_nd, charge_i_nd, np.real(ux[2,:]), skin_nd, beta)

    By_vals = np.zeros(len(xs)-1); Bz_vals = np.zeros(len(xs)-1)
    xs_B = xs[:-1] + (xs[1]-xs[0])/2
    for i in range(len(xs)-1):
      By_vals[i] = integrate.simps(dBydx_vals[:i+2], xs[:i+2])*1e5 # numerically integrte using simposons rule from smaple
      Bz_vals[i] = integrate.simps(dBzdx_vals[:i+2], xs[:i+2])*1e5 # numerically integrte using simposons rule from smaples

    #By_vals = integrate.simps(dBydx_vals, xs) # numerically integrte using simposons rule from smaples 
    #Bz_vals = integrate.simps(dBzdx_vals, xs) # numerically integrte using simposons rule from smaples 

    fig_d = plt.figure(figsize=(16,20), dpi=300); ax_d = [] # derived solutions

    # titles for each plot 
    titles = [r"$B_y\cdot 10^5$", r"$B_z\cdot 10^5$", r"$J_y/nV$", r"$J_z/nV$"]; 
    propList = [By_vals, Bz_vals, Jy_vals, Jz_vals]
    xList = [xs_B, xs_B, xs, xs]

    # polyfits     # Polynomial fits for simulation
    dataFits = []
    dataFitList = [By_vals, Bz_vals]
    degList = [3]*2
    for i in range(len(dataFitList)):
      dataFits.append(np.polyfit(xList[i], dataFitList[i], 2))
      print(dataFits[-1])

    for i in range(len(titles)):
      ax_d.append(fig_d.add_subplot(2,2,i+1))
      print(f"OK {i}\t{titles[i]}")
      ax_d[i].plot(xList[i], propList[i])
      
      if i < 2:
        ax_d[i].plot(xs_B, np.poly1d(dataFits[i])(xs_B), label='fit')
        ax_d[i].legend()
      ax_d[i].set_ylabel(titles[i])
      if i < 2: ax_u[i].set_xticks([])
      else: ax_u[i].set_xlabel('x')

    fig_u.tight_layout(pad=10)

    #=================================================== Plot all graphs ================================================#
    plt.show()




    # print checkcs   
    if False:
      print("\nThe coefficient matrix is")
      printArray(A)
  
      print("The eigenvalues and eigenvectors of A are:")
      stringOutput = "\t" 
      for i in range(A.shape[0]):
        stringOutput += f"{Eigenvalues[i]:10.4E}\t"
      print(stringOutput, "\n")
      printArray(Q)
  
      print("\nThe diagonal eigenmatrix is")
      printArray(Lambda)
  
      print("\nThe product Q Lambda Q^{-1} is:\n")
      printArray(A_back)
  
      print("Relativ error in the Eigendecomposition is:\n")
      printArray( (A-A_back)/A_back )
      pdb.set_trace()
      print("Left boundary conditions")
      print("\tU")
      printArray(U_left)
      print("\ts")
      printArray(s_left)
      print("Right boundary conditions")
      print("\tU")
      printArray(U_right)
      print("\ts")
      printArray(s_right)
  

