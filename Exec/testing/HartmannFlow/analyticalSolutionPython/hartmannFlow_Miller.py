import numpy as np
import pdb 
from scipy.optimize import fsolve
from sympy import Symbol, solve, nsolve
import matplotlib.pyplot as plt


#============================================================================================#
def coeffMatrix():
  A = np.zeros((4,4), dtype=complex)

  global nu_p_tau, omega_c_tau, B_field;
  global Z_e, nu_e, n_e, P_e, A_e;
  global Z_i, nu_i, n_i, P_i, A_i;

  Bx =  B_field[0]

  lambda_e = nu_p_tau * omega_c_tau * Z_e*nu_e*n_e*Bx/ P_e
  lambda_i = nu_p_tau * omega_c_tau * Z_i*nu_i*n_i*Bx/ P_i

  gamma_ei = nu_p_tau**2 * A_e * n_e * nu_ei * nu_e/P_e
  gamma_ie = nu_p_tau**2 * A_i * n_i * nu_ie * nu_i/P_i
  
  A[0,0] = gamma_ei; A[0,1] =  lambda_e; A[0,2] = -gamma_ei; A[0,3] = 0     
  A[1,0] =-lambda_e; A[1,1] =  gamma_ei; A[1,2] =  0       ; A[1,3] = -gamma_ei
  A[2,0] =-gamma_ie; A[2,1] =  0       ; A[2,2] =  gamma_ie; A[3,3] = lambda_i
  A[3,0] = 0       ; A[3,1] = -gamma_ie; A[3,2] = -lambda_i; A[3,2] = gamma_ie
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

def genSol(inputs):
  global Sil, Sir, Di, Xl, Xr;  
  #print(Sil, Sir, Di, Xl, Xr)
  C0 = inputs[0]+inputs[1]*1j
  C1 = inputs[2]+inputs[3]*1j

  #return np.array([Sil - C0*np.exp(np.sqrt(Di)*Xl) - C1*np.exp(-np.sqrt(Di)*Xl),
  #                 Sir - C0*np.exp(np.sqrt(Di)*Xl) - C1*np.exp(-np.sqrt(Di)*Xl)])

  val1 = Sil - C0*np.exp(np.sqrt(Di)*Xl) - C1*np.exp(-np.sqrt(Di)*Xl)
  val2 = Sir - C0*np.exp(np.sqrt(Di)*Xr) - C1*np.exp(-np.sqrt(Di)*Xr)

  return np.array([np.real(val1), np.imag(val1), np.real(val2), np.imag(val2)])

def sGenSol(inputs):
  global Sil, Sir, Di, Xl, Xr;
  C0 = inputs[0]
  C1 = inputs[1]
  x  = inputs[2]
  return C0*np.exp(np.sqrt(Di)*x) + C1*np.exp(-np.sqrt(Di)*x)

#===================================================================================================#
if __name__ == '__main__':
  # Set Miller parameters 

  global nu_p_tau, omega_c_tau, B_field;
  global Z_e, nu_e, n_e, P_e, A_e;
  global Z_i, nu_i, n_i, P_i, A_i;

  nu_p_tau = 1e2 ; omega_c_tau = 1; omega_p_tau = 1e2; 
  Z_e = -1; n_e = 1;P_e = 1; A_e = 1e-2;
  Z_i =  1; n_i = 1;P_i = 1; A_i = 1   ;

  nu_e = 1e1;
  nu_i = 1e1;
  nu_ei = 1e1;
  nu_ie = 1e1;
  B_field = np.array([1,0,0])

  #standard :
  #n0_nd = 1.0; Bx0_nd = 0.; By0_nd = 0.0; Bz0_nd = 0; T0_nd = 0.5;

  #Hartman flow parameters 
  v_plate = 0.01 ; # dimensionless
  Bx0_nd = 0.1 ; # dimensionless 
  Xdomain = 1;
  Ydomain = 1;

  #checks according to model constrains by Milelr i.e. region of applicabilty
  T_e = P_e / n_e ; T_i = P_i / n_i;
  if not nu_p_tau > 10*omega_c_tau: xxx
  if not nu_p_tau >= 100: xxx
  if not nu_p_tau <= 1e5: xxx
  if not omega_p_tau > 10*omega_c_tau: xxx
  if not 10*v_plate < min(np.sqrt(T_e/A_e), np.sqrt(T_i/A_i)): 
    print(f"v_plate = {v_plate}")
    print(f"v_e = {np.sqrt(T_e/A_e)}\tv_i = {np.sqrt(T_i/A_i)}") 
    xxx
  #===============================================Solve======================================#
  A = coeffMatrix()
  Eigenvalues, Q = np.linalg.eig(A)
  Lambda = np.diag(Eigenvalues)
  A_back = Q @ Lambda @ np.linalg.inv(Q)

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
  for i in range(1,5):
    guesses.append(np.array([1/i*(-1)**i]*4))

  for g in guesses:
    print(f"\nGuess = {g}")
    for i in range(Eigenvalues.shape[0]):
      Sil = s_left[i]; Sir = s_right[i]; Di = Eigenvalues[i];
      #eq1 = Sil - C0*np.exp(np.sqrt(Di)*Xl) - C1*np.exp(-np.sqrt(Di)*Xl)
      #eq2 = Sir - C0*np.exp(np.sqrt(Di)*Xl) - C1*np.exp(-np.sqrt(Di)*Xl)
  
      #sSol = nsolve( (eq1, eq2), (C0, C1), (0, 0) )
      sSol = fsolve( genSol, g)#, full_output=True)
      #print(f"\nSolution:\t{sSol[0]}+{sSol[1]}j\t {sSol[2]}+{sSol[3]}j")
      cArray[i,0] = sSol[0]+sSol[1]*1j; cArray[i,1] = sSol[2]+sSol[3]*1j;
      #print(f"Stored value:\t{cArray[i,0]}\t{cArray[i,1]}")
    xs = np.linspace(Xl, Xr, 1000)
    
    sx = np.zeros((s_left.shape[0], xs.shape[0]), dtype=complex)

    for i in range(s_left.shape[0]):
      sx[i,:] = sGenSol([cArray[i,0], cArray[i, 1], xs])
    
    # make a figure to plot as we get the data 
    fig = plt.figure(figsize=(16,20), dpi=300); ax = []
    # titles for each plot 
    titles = [r"$u^e_y$", r"$u^e_z$", r"$u^i_y$", r"$u^i_z$"]
    ux = np.zeros(sx.shape, dtype=complex)
    for j in range(xs.shape[0]):
      ux[:,j] = Q @ sx[:,j]
  
    for i in range(s_left.shape[0]):
      ax.append(fig.add_subplot(2,2,i+1))
      ax[i].plot(xs, np.real(ux[i,:]))
      ax[i].set_ylabel(titles[i])
      if i < 2: ax[i].set_xticks([])
      else: ax[i].set_xlabel('x')

    fig.tight_layout(pad=10)
    plt.show()
  
    if True:
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

