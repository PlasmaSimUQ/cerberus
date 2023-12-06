#!/usr/bin/env/ python3
# -*- coding: utf-8 -*-
"""
Created Tuesday 20200421 14:00

Adapted from numericalPostProcessing_v18. 

Purpose of this script is to serve as a centralised store of all functions 
commonly used across all my post processing functions for deriving properties
fron the PHM MFP solver writen by Dr. Daryl Bond.

Updated on 20200731 to be compatible with the June 11 push of cerberus by Daryl. 

Updated 20211111 to use the Boxlib file format. 
"""
#===============================================================================
import os, sys, gc, numpy as np # standard modules

visulaisation_code_folder ="/home/kyriakos/Documents/Code/000_cerberus_dev/githubRelease-cerberus/cerberus/vis"

if visulaisation_code_folder not in sys.path:
  sys.path.insert(0, visulaisation_code_folder)

#from get_hdf5_data_kyri_edit import ReadBoxLib
#from get_hdf5_data import ReadBoxLib
from get_boxlib import ReadBoxLib, get_files

import h5py, pylab as plt, matplotlib.gridspec as gridspec, matplotlib as mpl
#mpl.use('agg')
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from multiprocessing import Pool
import pdb, math
import copy 
#===============================================================================
def interp_val(input_alpha, t_low, t_high): # interperate value 
  if input_alpha < 0.5:
    return t_low
  else:   
    return t_high

def get_interface(ch, name):
  #x, alpha = ch.get_flat("alpha-%s"%name) # problem when reading for the first time without intiialiing get of something first 
  x, alpha = ch.get("alpha-%s"%name)

  if alpha.max() > 1. and  alpha.min() > 0.5: pdb.set_trace()
  y= x[1]; x = x[0]; dx = x[1] - x[0]; dy = y[1] - y[0]
  n_x = np.shape(alpha)[0]; n_y = np.shape(alpha)[1]
  #get interface stats for search
  interface_x = 0.0; interface_amp = 0.1; interface_transition = 0.01;
  #averaged interface stats 
  alpha_y_collapsed = np.sum(alpha, 1)  # find more accurate integration, inbuilt funcitons maybe 
  alpha_y_collapsed = alpha_y_collapsed/n_y #p.shape(alpha)[1]  #max_alpha_y_collapsed
  tracer_tol = 0.45
  tracer_low = 0.5 - tracer_tol
  tracer_high = 0.5 + tracer_tol
  # ------------------------Volume of fluid interface
  grab3 = (alpha_y_collapsed > tracer_low) & (alpha_y_collapsed < tracer_high)
  
  avg_int_coords = []; start_switch = 0;
  for i in range( np.shape(grab3)[0]):
    if grab3[i] == True and start_switch ==0: #positive transition
      start_switch = 1.
      int_start = i
      interp_val_1 = interp_val(alpha_y_collapsed[i], tracer_low, tracer_high)

    if grab3[i] == False and start_switch == 1.:
      int_end = i
      start_switch = 0.
      interp_val_2 = interp_val(alpha_y_collapsed[i], tracer_low, tracer_high)
      x1 = np.interp([interp_val_1], alpha_y_collapsed[int_start-1:int_start+1], x[int_start-1:int_start+1])
      x2 = np.interp([interp_val_2], alpha_y_collapsed[int_end-1:int_end+1], x[int_end-1:int_end+1])

      avg_int_coords.append([x1,x2])
      del int_start, int_end, x1, x2 ; gc.collect() ;
      if len(avg_int_coords[-1])%2 != 0:
        print("Interface tracking is defective") 
        print("time\n", ch.time )

  # ------------------------ eah row
  interface_tracking = {}; #np.zeros((n_y, 1)) # np.zeros((n_y, 2)) 
  for i in range(n_y):
#    interface_tracking[i] = {}
    interface_tracking[i] = []
    
  global_interface_start = x[-1]
  global_interface_end   = x[0]
  for j in range(n_y):
    grab3 = (alpha[:,j] > tracer_low) & ( alpha[:,j] < tracer_high )
    start_switch = 0;     
    for i in range( np.shape(grab3)[0]):
      if grab3[i] == True and start_switch ==0: #positive transition
        start_switch = 1.
        int_start = i
        interp_val_1 = interp_val(alpha[i,j], tracer_low, tracer_high)
      if grab3[i] == False and start_switch == 1.:
        int_end = i
        start_switch = 0.  
        interp_val_2 = interp_val(alpha[i,j], tracer_low, tracer_high) #TODO check this change
        x1 = np.interp([interp_val_1], alpha[int_start-1:int_start+1,j], x[int_start-1:int_start+1])
        x2 = np.interp([interp_val_2], alpha[int_end-1:int_end+1,j], x[int_end-1:int_end+1])

        interface_tracking[j].append(x1)
        interface_tracking[j].append(x2)
        del int_start, int_end ; gc.collect()
    if len(interface_tracking[j])%2 != 0:
      print("Interface_tracking is broken")
      print("time\n", ch.time)
      pdb.set_trace()      
    if len(interface_tracking[j]) == 0:
      print("No interface on row: ", j, "height: ", y[j] )
      print("time\n", ch.time)
      pdb.set_trace()
    try:
      global_interface_start = min(min(interface_tracking[j]), global_interface_start)
      global_interface_end   = max(max(interface_tracking[j]), global_interface_end)
    except:
      continue

  return x, y, avg_int_coords, global_interface_start, global_interface_end, interface_tracking 
  #avg_interface_start, avg_interface_end, global_interface_start, global_interface_end, interface_tracking

def get_growth(int_width_n, int_width_nm1, t_n, t_nm1):
  """
  calculate the growth rate from the nth and n+1 data file --- each file represents
  a time step. Uses backward difference, chosen for ease of use.
  """
  return (int_width_n - int_width_nm1)/(t_n-t_nm1)

def get_baroclinic_torque(ch, name):
    
    x, rho = ch.get_flat("rho-%s"%name)
    try:
        x, prs = ch.get_flat("p-%s"%name)
    except:
        x, y, mx = ch.get_flat("x_mom-%s"%name)
        x, y, my = ch.get_flat("y_mom-%s"%name)
        x, y, mz = ch.get_flat("z_mom-%s"%name)
        x, y, nrg = ch.get_flat("rho_E-%s"%name)
        x, y, gam = ch.get_flat("gamma-%s"%name)
        prs = (nrg - (mx**2 + my**2 + mz**2)/(2*rho))*(gam - 1.0)
    
    y = x[1]; x = x[0]; dx = x[1] - x[0]; dy = y[1] - y[0]

    rho_dx = np.gradient(rho, dx, axis=0)
    rho_dy = np.gradient(rho, dy, axis=1)
    
    prs_dx = np.gradient(prs, dx, axis=0)
    prs_dy = np.gradient(prs, dy, axis=1)
    
    baro = (1/rho**2)*(rho_dx*prs_dy - rho_dy*prs_dx)
    
    return x, y, baro

#### Old standard get_Lorentz_torque, there is another strictly for get_Lorentz
def get_Lorentz_torque(ch, name):
    
    x, rho = ch.get_flat("rho-%s"%name)
    try:
        u = ch.get_flat("x_vel-%s"%name)[1]
        v = ch.get_flat("y_vel-%s"%name)[1]
        w = ch.get_flat("z_vel-%s"%name)[1]
    except:
        x, mx = ch.expression("{x-mom-%s}"%name)
        x, my = ch.expression("{y-mom-%s}"%name)
        x, mz = ch.expression("{z-mom-%s}"%name)
        
        u = mx/rho; v = my/rho; w = mz/rho
    
    Ex = ch.get_flat("x_D-field")[1]; Ey = ch.get_flat("y_D-field")[1]; 
    Ez = ch.get_flat("z_D-field")[1]
    
    Bx = ch.get_flat("x_B-field")[1]; By = ch.get_flat("y_B-field")[1]
    Bz = ch.get_flat("z_B-field")[1]
     
    mass = ch.get("mass-"+name)[1]; charge = ch.get("charge-"+name)[1]
    
    beta = ch.data["beta"]; dS = ch.data["skin_depth"]; c = ch.data["lightspeed"]

    # find Lorentz torque component from electric and magnetic field contributions
    Lx_E = np.sqrt(2/beta)*charge/(mass*dS)*c*Ex
    Ly_E = np.sqrt(2/beta)*charge/(mass*dS)*c*Ey 
    Lx_B = np.sqrt(2/beta)*charge/(mass*dS)*(v*Bz - w*By)
    Ly_B = np.sqrt(2/beta)*charge/(mass*dS)*(w*Bx - u*Bz)

    Lx = np.sqrt(2/beta)*charge/(mass*dS)*(c*Ex + v*Bz - w*By)
    Ly = np.sqrt(2/beta)*charge/(mass*dS)*(c*Ey + w*Bx - u*Bz)
    y= x[1]; x = x[0]; dx = x[1] - x[0]; dy = y[1] - y[0]

    dLx_E_dy = np.gradient(Lx_E, dy, axis=1)
    dLy_E_dx = np.gradient(Ly_E, dx, axis=0)
    dLx_B_dy = np.gradient(Lx_B, dy, axis=1)
    dLy_B_dx = np.gradient(Ly_B, dx, axis=0)
  
    dLx_dy = np.gradient(Lx, dy, axis=1)
    dLy_dx = np.gradient(Ly, dx, axis=0)

    tau_E = dLy_E_dx - dLx_E_dy
    tau_B = dLy_B_dx - dLx_B_dy
    tau = dLy_dx - dLx_dy
    
    if abs(np.sum(tau) - (np.sum(tau_E) + np.sum(tau_B))) > 0.1:
      print('EM torque do not match')
      print('np.sum(tau)',np.sum(tau))
      print('np.sum(tau_E)+...',np.sum(tau_E) + np.sum(tau_B))
      xxxx

    if True: #non constnat q and m values 
      dyq = np.gradient(charge, y, axis=1)
      dyminv = np.gradient(mass**-1, y, axis=1)
      dxq = np.gradient(charge, x, axis=0)
      dxminv = np.gradient(mass**-1, x, axis=0)
      tauE_qm = np.sqrt(2/beta)/dS*\
        ((mass**-1*dxq + charge*dxminv)*Ey - \
          (mass**-1*dyq - charge*dyminv)*Ex )
      tauB_qm = np.sqrt(2/beta)/dS*\
        ((mass**-1*dxq + charge*dxminv)*(w*Bx - u*Bz) -\
          (mass**-1*dyq - charge*dyminv)*(v*Bz - w*By) )
      return x, y, tau_E + tauE_qm, tau_B+tauB_qm, tau
    else: print("\tWarning Lorentz torque m, q assumed constant")
    return x, y, tau_E, tau_B, tau

def get_vorticity(ch, input_options):#name):
    """Calculates the vorticity (omega) and change in vorticity due to the 
    omega*(grad dot u) {flow1} and (u dot grad)omega {flow2} terms  
    added flow3 for (omega cdot nabla) u
    """
    if 'dim' not in input_options.keys(): dim = 2 # z direction
    else: 
      dim = input_options['dim']
      #TODO
      #print("\t\tNote the 2D slab assumption for x- and y-vorticity - zero z gradient")

    name = input_options['name']
    quantity = input_options['quantity']
    try:
      x, rho = ch.get("rho-%s"%name)
    except:
      print(f"Component not available {'rho-%s'%name}")
      print(str(np.ravel(ch.data["names"])))
      pdb.set_trace()

    u = ch.get_flat("x_vel-%s"%name)[1]
    v = ch.get_flat("y_vel-%s"%name)[1]
    w = ch.get_flat("z_vel-%s"%name)[1]

    y = x[1]; x = x[0]; dx = x[1] - x[0]; dy = y[1] - y[0]

    if dim == 0: # vorticity-x
      try: # check if gradients given in output data
        dw_dy = ch.get_flat("z_vel-%s-dy"%name)[1]
        #dv_dz = ch.get_flat("y_vel-%s-dz"%name)[1]
      except:
        dw_dy = np.gradient(w, dy, axis=1)
        #dv_dz = np.gradient(v, dz, axis=2)
      omega = dw_dy #- dv_dz
    elif dim == 1: # vorticity-y
      try: # check if gradients given in output data
        dw_dx = ch.get_flat("z_vel-%s-dx"%name)[1]
        #du_dz = ch.get_flat("x_vel-%s-dz"%name)[1]
      except:
        dw_dx = np.gradient(w, dx, axis=0)
        #du_dz = np.gradient(u, dz, axis=2)
      omega = -(dw_dx )  #- du_dz)
    elif dim == 2: # vorticity-z
      try: # check if gradients given in output data
        du_dy = ch.get_flat("x_vel-%s-dy"%name)[1]
        dv_dx =  ch.get_flat("y_vel-%s-dx"%name)[1]
      
        du_dx =  ch.get_flat("x_vel-%s-dx"%name)[1]
        dv_dy =  ch.get_flat("y_vel-%s-dy"%name)[1]
      except:
        du_dy = np.gradient(u, dy, axis=1)
        dv_dx = np.gradient(v, dx, axis=0)
      
        du_dx = np.gradient(u, dx, axis=0)
        dv_dy = np.gradient(v, dy, axis=1)
  
      omega = dv_dx - du_dy
    if quantity == 'omega':
      return x, y, omega

    # if return has not been reached then calculate flow1 and flow2 assuming z dimension
    if dim != 2: xxxx # don't support non z dimension flow1 and flow 2 
    # velocity derivatives for vorticity 
    try:
      du_dx =  ch.get_flat("x_vel-%s-dx"%name)[1]
      dv_dy =  ch.get_flat("y_vel-%s-dy"%name)[1]
    except:
      du_dx = np.gradient(u, dx, axis=0)
      dv_dy = np.gradient(v, dy, axis=1)
    d_omega_dx = np.gradient(omega, dx, axis=0)
    d_omega_dy = np.gradient(omega, dy, axis=1)
    flow1 = -omega*(du_dx + dv_dy)
    flow2 = -u*d_omega_dx - v*d_omega_dy # not used usually because this is just transport 

    # velocity derivatives for the x and y components of vorticity

    try: # check if gradients given in output data
      dw_dy = ch.get_flat("z_vel-%s-dy"%name)[1]
      #dv_dz = ch.get_flat("y_vel-%s-dz"%name)[1]
      dw_dx = ch.get_flat("z_vel-%s-dx"%name)[1]
      #du_dz = ch.get_flat("x_vel-%s-dz"%name)[1]
    except:
      dw_dy = np.gradient(w, dy, axis=1)
      #dv_dz = np.gradient(v, dz, axis=2)
      dw_dx = np.gradient(w, dx, axis=0)
      #du_dz = np.gradient(u, dz, axis=2)


    omega_x = dw_dy #- dv_dz
    omega_y = -(dw_dx )  #- du_dz)

    flow3 = omega_x*dw_dx + omega_y*dw_dy
    return x, y, omega, flow1, flow2, flow3

def get_charge_number_density(rc, name, returnMany = True):
    """
    Inputs --- data h5 data construct holding simulation output data *see 
               get_hdf5_data.py for syntax and h5data construct access* and name, 
               'ion', 'electron', 'neutral'
    Outputs --- depending on spatial dimensions --- x, [y, z], charge-rho, x-current rho, 
                y-current rho...
    """
    em_prop = {}
    fluid_prop = {}
    
    fluid_prop['mass-%s'%name] = rc.get('mass-%s'%name)[-1]
    em_prop['charge-%s'%name] = rc.get('charge-%s'%name)[-1]
    specialBoi = rc.get("rho-%s"%name)
    fluid_prop['rho-%s'%name] = specialBoi[-1]

    try:
      fluid_prop['x_vel-%s'%name] = rc.get_flat("x_vel-%s"%name)[-1]
      fluid_prop['y_vel-%s'%name] = rc.get_flat("y_vel-%s"%name)[-1]
      fluid_prop['z_vel-%s'%name] = rc.get_flat("z_vel-%s"%name)[-1]
    except:
      mx = rc.get_flat("x-mom-%s"%name)[-1]
      my = rc.get_flat("y-mom-%s"%name)[-1]
      mz = rc.get_flat("z-mom-%s"%name)[-1]
      
      fluid_prop['x_vel-%s'%name] = mx/fluid_prop['rho-%s'%name]
      fluid_prop['y_vel-%s'%name] = my/fluid_prop['rho-%s'%name]
      fluid_prop['z_vel-%s'%name] = mz/fluid_prop['rho-%s'%name]
    
    fluid_prop['n-%s'%name] = np.divide(fluid_prop['rho-%s'%name], fluid_prop['mass-%s'%name])

    charge_rho = np.multiply(em_prop['charge-%s'%name], fluid_prop['n-%s'%name])
    current_rho_x = np.multiply( fluid_prop['x_vel-%s'%name], charge_rho) 
    current_rho_y = np.multiply( fluid_prop['y_vel-%s'%name], charge_rho) 
    current_rho_z = np.multiply( fluid_prop['z_vel-%s'%name], charge_rho) 
    if returnMany:
      return specialBoi[0][0], charge_rho, fluid_prop['x_vel-%s'%name], fluid_prop['y_vel-%s'%name], fluid_prop['z_vel-%s'%name], current_rho_x, current_rho_y, current_rho_z
    else:
      return specialBoi[0][0], charge_rho, current_rho_x, current_rho_y, current_rho_z

def get_pressure(ch, inputs):
    name = inputs['name']
    x, rho = ch.get("rho-%s"%name)
    try:
        x, prs = ch.get("p-%s"%name)
    except:
        x, mx = ch.get("x_mom-%s"%name); x, my = ch.get("y_mom-%s"%name); 
        x, mz = ch.get("z_mom-%s"%name)
        x, nrg = ch.get("rho_E-%s"%name); x, gam = ch.get("gamma-%s"%name)
        prs = (nrg - (mx**2 + my**2 + mz**2)/(2*rho))*(gam - 1.0)
    return x[0], x[1], prs

def get_single_data(din):
    save_name = get_save_name(din['key'], din['folder'], din['level'], os.path.split(din["dataName"])[1] )
    braginskiiVorticity = din['braginskiiVorticity']
    bragVortRhoVar = din['bragVortRhoVar']
    calcIntra = din['calcIntra']; calcInter = ['calcInter']; isoSwitch = din['isoSwitch']
    #print( save_name )
    if os.path.isfile(din['dir_name']+"/" +save_name):
        print(os.path.split(save_name)[1]," already exists")
        return

    rc = ReadBoxLib(din["dataName"], max_level=din["level"], limits=din["window"]) 
    #interface info     
    constOC = (2/rc.data['beta']/rc.data['skin_depth']**2)**0.5
     
    data = {}
    for name in rc.names:
        if 'field' not in name: data[name] = {} # ignore field state 

    print("Processing: ", os.path.split(din["dataName"])[1]," @ ",rc.time) 

    #### Extract collisional terms
    if braginskiiVorticity: 
      #print("Viscous vorticity contribution")
      #TODO Change the axis since the web page i was on seems to say axis=0     takes gradient acrosos rows???
      if isoSwitch: print("\t\tIso overide ON! Isotropic modelling ")
      else: print("\t\tIso overide OFF! Anisotropic modelling ")

      #TODO reduce the search area and offset indexes
      x, y, dudt_fluxes, srcDst = get_transportProperties(rc, 
        ["ions", "electrons"], din['level'], calcIntra, calcInter, 
        isoOveride=isoSwitch, useNPROC=1)

      rc = ReadBoxLib(din["dataName"], max_level=din["level"], limits=din["window"]) #reopen after operatons ar done for the remaineder
      fluxDict = {}; srcDict = {}; OC = {}# vorticity vontribution 
      for name in ['ions', 'electrons']: #delete useless data
        fluxDict[name] = {}
        srcDict[name] = {}
        # divide by density 
        rhoName = rc.get("rho-%s"%name)[-1]

        if 'ion' in name: nameStr = '_i_'
        elif 'electron' in name: nameStr = "_e_"

        dxrhoinv = np.gradient(rhoName**-1, x, axis=0)
        dyrhoinv = np.gradient(rhoName**-1, y, axis=1)

        if bragVortRhoVar: 
          if calcInter: OC['R'+nameStr+'rho'] = constOC*\
            (dxrhoinv*(srcDst[name][:,:,1]) - dyrhoinv*(srcDst[name][:,:,0]))
          if calcIntra: OC['PI'+nameStr+'rho'] = constOC*\
            (dxrhoinv*(dudt_fluxes[name][:,:,1]) -\
            dyrhoinv*(dudt_fluxes[name][:,:,0]))
        else:
          if calcInter: OC['R'+nameStr+'rho'] = np.zeros(srcDst[name][:,:,0].shape)  
          if calcIntra: OC['PI'+nameStr+'rho'] = np.zeros(dudt_fluxes[name][:,:,1].shape)

        if calcIntra:
          for key in range(dudt_fluxes[name].shape[-1]):
            if key in [0, 1]: # 0:Xmom, 1:Ymom
              fluxDict[name][key] = np.copy(dudt_fluxes[name][:,:,key])/rhoName*constOC
        if calcInter:
          for key in range(srcDst[name].shape[-1]):
            if key in [0, 1]: # 0:Xmom, 1:Ymom
              srcDict[name][key] = np.copy(srcDst[name][:,:,key])/rhoName*constOC

        del rhoName;
      if calcIntra: del dudt_fluxes, constOC; gc.collect()
      if calcInter: del srcDst; gc.collect()

      if calcInter:
        OC['R_i'] = np.gradient(srcDict['ions'][1], x, axis=0) - \
              np.gradient(srcDict['ions'][0], y, axis=1)
        OC['R_e'] = np.gradient(srcDict['electrons'][1], x, axis=0) - \
              np.gradient(srcDict['electrons'][0], y, axis=1)
        del srcDict; gc.collect()
      if calcIntra:
        OC['PI_i'] = np.gradient(fluxDict['ions'][1], x, axis=0) - \
              np.gradient(fluxDict['ions'][0], y, axis=1)
        OC['PI_e'] = \
          np.gradient(fluxDict['electrons'][1], x, axis=0) - \
          np.gradient(fluxDict['electrons'][0], y, axis=1)
        del fluxDict; gc.collect();

      #if rc.time > 0. and True: pdb.set_trace()
      
    ####### Loop for each state and interace 
    for name, d in data.items():
      tracerDefined = False
      
      try: # check for tracer and find interface if present
        tracer = rc.get("alpha-%s"%name)[-1]
        tracerDefined = True
      except:
        print(f"\ttracer not available for {name}")
      
      get_vorticity_inputs = {}
      get_vorticity_inputs['name'] = name
      get_vorticity_inputs['quantity'] ='all'

      x, y, omega, flow1, flow2, flow3 = get_vorticity(rc, get_vorticity_inputs)
      get_vorticity_inputs['quantity'] = 'omega'; get_vorticity_inputs['dim'] = 0; 
      deleteMe, deleteMe1, omega_x = get_vorticity(rc, get_vorticity_inputs)
      get_vorticity_inputs['dim'] = 1; 
      deleteMe, deleteMe1, omega_y = get_vorticity(rc, get_vorticity_inputs)
      del deleteMe, deleteMe1; gc.collect()

      x, y, baro = get_baroclinic_torque(rc, name)
      dx = x[1] - x[0]; dy = y[1] - y[0] # note each level is constant in grid 

      #print("\n\n#==========Note hardcoded ignore electrons in interface anaylsis ln 450 PHM ")
      if tracerDefined:# and name != 'electrons':
        # int tracking 
        x, y, avg_int_coords, global_interface_start, global_interface_end, \
        interface_tracking = get_interface(rc, name)    
        if din['interfaceHeuristic']: # When brag pushes back interface use this 
          xx, rho = rc.get("rho-%s"%name)
          yy = xx[1]; xx = xx[0]
          drhodx = np.gradient(rho,xx,axis=0)
          drhody = np.gradient(rho,yy,axis=1)
          drho_mag = np.sqrt(drhodx**2 + drhody**2)
          if 'ion' in name or 'electron' in name:
            chargeDensity = get_charge_number_density(rc, "electrons", False)[1] +\
                            get_charge_number_density(rc, "ions", False)[1]
          else:
            print('Neutral fluid - mock cd search with ')
            chargeDensity = drho_mag.copy()
           
          del drhodx, drhody, rho; gc.collect();
          #print("\tInterface Heutristic active")
          for j in range(len(interface_tracking)):
            if int(len(interface_tracking[j])/2) > 1: print("\t###More than one transition")
            for k in range(int(len(interface_tracking[j])/2)):
              #Fist pass at the interface heuristic 
              iStart = int((interface_tracking[j][2*k]-xx[0])/dx)
              iEnd = int((interface_tracking[j][2*k+1]-xx[0])/dx)
              #print(f"\n\t##j:\t{j}\nDefault interval:\t{iStart} {iEnd}")
              peak_drho = 0
              drhoAvgWin = (drho_mag[iStart-1:iEnd,j] + drho_mag[iStart:iEnd+1,j] + drho_mag[iStart+1:iEnd+2,j])/3. # rolling avg of three, centred 
              cdWin = (chargeDensity[iStart-1:iEnd,j] + \
                chargeDensity[iStart:iEnd+1,j] + chargeDensity[iStart+1:iEnd+2,j])/3.
              iSymmetry = np.argmax(drhoAvgWin)  # find max 
              iSymmetryCD = np.argmax(np.abs(cdWin)) #charge density average 
              #print(f"dRho criteria interval:\t{iStart} {iStart + 2*iSymmetry+1}")
              #print("\n\n#================Interface heurist drho and cd search off for HRMI")
              if True: #TODO add switch 
                triggerVal = 0.025
                iBuffer = 10 # offset to check region outside in case not acpturing all 
                iEnd = (iStart) + 2*iSymmetry + 1  + iBuffer
                iEndCD = (iStart) + 2*iSymmetryCD + 1 + iBuffer
                iStart -= iBuffer
                #update the search range
                drhoAvgWin = (drho_mag[iStart-1:iEnd,j] + \
                  drho_mag[iStart:iEnd+1,j] + drho_mag[iStart+1:iEnd+2,j])/3.
                cdWin = (chargeDensity[iStart-1:iEndCD,j] + \
                  chargeDensity[iStart:iEndCD+1,j] + \
                  chargeDensity[iStart+1:iEndCD+2,j])/3.
                iSymmetry = np.argmax(drhoAvgWin) #density average 
                iSymmetryCD = np.argmax(np.abs(cdWin)) #density average 
                iN_cd = np.argmin(cdWin); iP_cd = np.argmax(cdWin); #peaks for cd
                peak_i = iSymmetry; peakCD_i = iSymmetryCD
                drhoAvg_threshold = drhoAvgWin[peak_i]*triggerVal
                cdAvgThreshold_N = np.abs(cdWin[iN_cd]*triggerVal) 
                cdAvgThreshold_P = np.abs(cdWin[iP_cd]*triggerVal) 
                #print(f"rhoAvg: {drhoAvg_threshold} , cdN: {cdAvgThreshold_N}, cdP: {cdAvgThreshold_P}")
                #update the start and finish 
                avgList = interfacePeak(drhoAvgWin, drhoAvg_threshold, iStart , \
                  iEnd, False)
                if cdAvgThreshold_N < 1e-10 or cdAvgThreshold_P < 1e-10: cdAvgList = []
                else: cdAvgList = \
                  interfacePeakSigned(cdWin, [cdAvgThreshold_N,cdAvgThreshold_P], iStart, \
                  iEndCD, iN_cd, iP_cd, False)

                #print(f"rho:\t{avgList}, cd:\t{cdAvgList}")
                if len(avgList) < 2 or len(cdAvgList) < 2:
                  finalStart = iStart; finalEnd = iStart + 2*iSymmetry+1
                else:
                  finalStart = min(min(avgList),min(cdAvgList)) 
                  finalEnd = max(max(avgList),max(cdAvgList))

                #Assign new interface 
                #if iEnd-iStart+1 > finalEnd-finalStart+1: 
                #  print(f"{j}\tinterface decreased")
                interface_tracking[j][2*k][0] = xx[finalStart]
                interface_tracking[j][2*k+1][0] = xx[finalEnd]
                #print(f"dRho and CD criteria interval:\t{finalStart} {finalEnd}\n")

              else: # otherwise just use the interface density gradeint only 
                #Assign new interface 
                if iEnd-iStart+1 >2*iSymmetry+1: print(f"{j}\tinterface decreased")
                else: interface_tracking[j][2*k+1][0]=xx[iStart + 2*iSymmetry+1]
        del drho_mag; gc.collect();

        #d["y_avg_int_start"] = avg_interface_start
        #d["y_avg_int_end"] = avg_interface_end 

        avg_all_points = []
        for i in range(len(avg_int_coords)):
          avg_all_points += avg_int_coords[i]

        if len(avg_all_points)==0:
          print(rc.time, "no interface detected")
        else:
          d["y_avg_int_width"] = max(avg_all_points) - min(avg_all_points) #avg_interface_end-avg_interface_start
#        d["global_int_start"]= global_interface_start
#        d["global_int_end"]= global_interface_end
        d["global_int_width"]= global_interface_end-global_interface_start
        # create attributes within 

        d["interface_location"]=interface_tracking #now an array (j,1) with each entry a list of lists that descrive each transition zone either 0.05 to 0.95 or 0.95 to 0.05
        #Hydrodynamic interface stats 
        omega_interface_sum = 0.; omega_interface_sum_half = 0.
        omega_interface_sum_x = 0.; omega_interface_sum_half_x = 0.
        omega_interface_sum_y = 0.; omega_interface_sum_half_y = 0.

        tb_interface_sum = 0.; tb_interface_sum_half = 0.;
        flow1_interface_sum = 0.; flow1_interface_sum_half = 0.; 
        flow2_interface_sum = 0.; flow2_interface_sum_half = 0.;
        flow3_interface_sum = 0.; flow3_interface_sum_half = 0.;

        interface_area = 0.
     
        if ("ion" in name) or ("electron" in name): #Lorent contributions
          # sum over the area of the interface
          x, y, tau_E, tau_B, tau = get_Lorentz_torque(rc, name)     
          tl_E_interface_sum = 0.; tl_E_interface_sum_half = 0.;
          tl_B_interface_sum = 0.; tl_B_interface_sum_half = 0.;
          tl_interface_sum = 0.; tl_interface_sum_half = 0.;
          #set Collisional sums to zero at the start of each fluid state 
          ##note the values for intra inter etc are calced already together 
          t_brag_intra_sum = 0; t_brag_intra_half = 0; #Storage 
          t_brag_inter_sum = 0; t_brag_inter_half = 0;
          if 'ion' in name: ieKey = '_i'; 
          else: ieKey = '_e'; 

        # changed to accomodate multilpe interface transition regions i.e. 
        # interface_tracking[j,0]= [[x1,x2], [x1,x2]]

        interface_area = 0.
          
        for j in range(len(interface_tracking)):
          for k in range(int(len(interface_tracking[j])/2)):
            for i in range(int((interface_tracking[j][2*k]-x[0])/dx), int((interface_tracking[j][2*k+1]-x[0])/dx)):
              interface_area += dx*dy
              omega_interface_sum += omega[i,j]
              omega_interface_sum_x += omega_x[i,j]
              omega_interface_sum_y += omega_y[i,j]
              tb_interface_sum += baro[i,j]
              flow1_interface_sum += flow1[i,j]
              flow2_interface_sum += flow2[i,j]
              flow3_interface_sum += flow3[i,j]
              if ("ion" in name) or ("electron" in name):
                tl_E_interface_sum += tau_E[i,j]
                tl_B_interface_sum += tau_B[i,j]
                tl_interface_sum += tau[i,j]
                if braginskiiVorticity: 
                  if calcIntra: t_brag_intra_sum += OC['PI'+ieKey][i,j] + \
                    OC['PI'+ieKey+'_rho'][i,j]
                  if calcInter: t_brag_inter_sum += OC['R'+ieKey][i,j] + \
                    OC['R'+ieKey+'_rho'][i,j]; 

              if j <= int(len(interface_tracking)/2):
                omega_interface_sum_half += omega[i,j]
                omega_interface_sum_half_x += omega_x[i,j]
                omega_interface_sum_half_y += omega_y[i,j]
                tb_interface_sum_half += baro[i,j]
                flow1_interface_sum_half += flow1[i,j]
                flow2_interface_sum_half += flow2[i,j]
                flow3_interface_sum_half += flow3[i,j]
                if ("ion" in name) or ("electron" in name):
                  tl_E_interface_sum_half += tau_E[i,j]
                  tl_B_interface_sum_half += tau_B[i,j]
                  tl_interface_sum_half += tau[i,j]
                  if braginskiiVorticity: 
                    if calcIntra: t_brag_intra_half += OC['PI'+ ieKey][i,j]+ \
                      OC['PI'+ieKey+'_rho'][i,j]; 
                    if calcInter: t_brag_inter_half += OC['R' + ieKey][i,j]+ \
                      OC['R'+ieKey+'_rho'][i,j];

          d["circulation_interface_sum"] = omega_interface_sum*dx*dy 
          d["circulation_interface_sum_x"] = omega_interface_sum_x*dx*dy 
          d["circulation_interface_sum_y"] = omega_interface_sum_y*dx*dy 
          d["baroclinic_interface_sum"] = tb_interface_sum*dx*dy
          d["tau_sc_interface_sum"] = flow1_interface_sum*dx*dy
          d["tau_conv_interface_sum"] = flow2_interface_sum*dx*dy
          d["tau_tilt_interface_sum"] = flow3_interface_sum*dx*dy
          d["circulation_interface_sum_half"] = omega_interface_sum_half*dx*dy 
          d["circulation_interface_sum_half_x"] = omega_interface_sum_half_x*dx*dy 
          d["circulation_interface_sum_half_y"] = omega_interface_sum_half_y*dx*dy 
          d["baroclinic_interface_sum_half"] = tb_interface_sum_half*dx*dy
          d["tau_sc_interface_sum_half"] = flow1_interface_sum_half*dx*dy
          d["tau_conv_interface_sum_half"] = flow2_interface_sum_half*dx*dy
          d["tau_tilt_interface_sum_half"] = flow3_interface_sum_half*dx*dy
          d["interface_area"] = interface_area 
          if ("ion" in name) or ("electron" in name):
            d["curl_Lorentz_E_interface_sum_half"] = tl_E_interface_sum_half*dx*dy
            d["curl_Lorentz_B_interface_sum_half"] = tl_B_interface_sum_half*dx*dy
            d["curl_Lorentz_interfaces_sum_half"] = tl_interface_sum_half*dx*dy
            d["curl_brag_inter_sum_half"] = t_brag_inter_half*dx*dy
            d["curl_brag_intra_sum_half"] = t_brag_intra_half*dx*dy

            d["curl_Lorentz_E_interface_sum"] = tl_E_interface_sum*dx*dy
            d["curl_Lorentz_B_interface_sum"] = tl_B_interface_sum*dx*dy
            d["curl_Lorentz_interfaces_sum"] = tl_interface_sum*dx*dy
            d["curl_brag_inter_sum"] =   t_brag_inter_sum*dx*dy 
            d["curl_brag_intra_sum"] =   t_brag_intra_sum*dx*dy

        if din["record_contour"] == True:
          #print("Record contour data:", din["dataName"] )
          d["vorticity"] = omega
          d["flow1"] = flow1
          d["flow2"] = flow2
          d["baroclinic"] = baro
          d["flow1_sum"] = np.sum(flow1*dx*dy)
          d["flow2_sum"] = np.sum(flow2*dx*dy)
          d["baroclinic_sum"] = np.sum(baro*dx*dy)

          if ("neutral" not in name):
            d["curl_Lorentz_E"] = tau_E
            d["curl_Lorentz_B"] = tau_B
            d["curl_Lorentz_E_sum"] = np.sum(tau_E*dx*dy)
            d["curl_Lorentz_B_sum"] = np.sum(tau_B*dx*dy)
            #x, qrho, cd_x, cd_y, cd_z = get_charge_number_density(rc, name, False)  # charge density, current_density_x, currrent_density_y
            #d['qrho_%s'%name] = qrho # chrage density
            #d['cd_%s_x'%name] = cd_x # current density
            #d['cd_%s_y'%name] = cd_y #current density
            #d['cd_%s_z'%name] = cd_y #current density

      if ("ion" in name) or ("electron" in name):
        del tau_E, tau_B, tau, tl_interface_sum, tl_B_interface_sum, tl_E_interface_sum, tl_interface_sum_half, tl_B_interface_sum_half, tl_E_interface_sum_half; #, cd_y, cd_x, qrho, 

    data["x"] = x
    data["y"] = y
    data["t"] = rc.time
    save = h5py.File(din['dir_name']+"/" + save_name, "w")# save data for this time 
    print("\tProcessed: ", os.path.split(din["dataName"])[1]," @ ",rc.time) 
    rc.close() ; del rc ; gc.collect() ;
    save_dict(data, save)
    save.close()

    del x, y, omega, omega_x, omega_y, flow1, flow2, baro, avg_int_coords, global_interface_start, global_interface_end, interface_tracking, omega_interface_sum, flow1_interface_sum, flow2_interface_sum, tb_interface_sum, save, data, interface_area
    gc.collect()

    return []

def save_dict(d, h5, include=[], exclude=[]):
      for key, value in d.items():
          if include:
              if key not in include:
                  continue
              elif type(value) == dict:
                  include += list(value.keys())
          elif key in exclude:
              continue
         
          if type(value) != dict:
              h5.create_dataset(str(key), data=value)
          else:
              grp = h5.create_group(str(key))
              save_dict(value, grp, include, exclude)
      return

def load_dict(d, h5, include=[], exclude=[]):
   
    for key, value in h5.items():
        if include:
            if key not in include:
                continue
            elif type(value) == h5py.Group:
                include += list(value.keys())
        elif key in exclude:
            continue
       
        if type(value) == h5py.Dataset:
            d[key] = value[()]
        else:
            d[key] = {}
            load_dict(d[key], value, include, exclude)
    return

"""
def get_save_name_old(folder, level, output_file = 0):
if output_file == 0:
  path, save_name = os.path.split(folder)
  if 'PHM' in folder:
    return folder + "_level=%i"%level + '_post_processed'
  else: 
    save_name = path.split("SRMI_")[-1].split("-DC")[0] + '-T-' + path.split("T-")[-1].split("\Res")[0]
  save_name += "_level=%i.h5"%level
  return save_name

else:
  if 'chk.' in output_file or 'plt.' in output_file:
    if 'chk' in output_file:
      step_id = output_file.split('chk.')[1].split('.hdf5')[0]
    elif 'plt' in output_file:
      step_id = output_file.split('plt.')[1].split('.hdf5')[0]
    else:   
      print("save_name  error in get_save_name: no chk or plt in name")
  else: 
    if 'chk' in output_file:
      step_id = output_file.split('chk')[1].split('.hdf5')[0]
    elif 'plt' in output_file:
      step_id = output_file.split('plt')[1].split('.hdf5')[0]
    else:   
      print("save_name  error in get_save_name: no chk or plt in name")


  path, save_name = os.path.split(folder)
  if 'PHM' in path or 'PHM' in save_name:
    save_name = save_name.split("PHM_")[-1] + '-T-1-Step'  
  else:
    save_name = path.split("SRMI_")[-1].split("-DC")[0] + '-T-' + path.split("T-")[-1].split("\Res")[0]

  save_name += "_" + step_id + "_level=%i.h5"%level
  return save_name
"""

def get_save_name(key, folder, level, output_file = 0):
    """
    Get the save name of the directoryor specific file belinging to the directory. This uses
    the new format were the ouput file directory nickname passed in is used as the prefix. 
    """
    if output_file == 0:
      return key + "_level=%i.h5"%level

    else:
      if 'chk' in output_file:
        step_id = output_file.split('chk')[1].split('.hdf5')[0]
      elif 'plt' in output_file:
        step_id = output_file.split('plt')[1].split('.hdf5')[0]
      else:
        print("save_name  error in get_save_name: no chk or plt in name")

      save_name = key + "_step" + step_id + "_level=%i.h5"%level
      return save_name

def get_batch_data(key, folder, level, max_res, window, n_increments, 
  nproc=1, outputType=["plt"], braginskiiVorticity=True, bragVortRhoVar=True, 
  interfaceHeuristic=True, calcIntra=True, calcInter=True, isoSwitch=True): 

    print("get_batch_data --- note the function relieson ion, electron, neutral, keywords in state names for interface statistics")
    """
    master function for processing data. Organises the individual processing of 
    time steps and creates the folder to store data. 
    """
    print(f"#=========Parameters==========#\n\tBraginskii vorticity:\t{braginskiiVorticity}", 
      f"\n\tRho variance: {bragVortRhoVar}\n\tInterface heuristic: {interfaceHeuristic}", 
      f"\n\tIntra collisions: {calcIntra}\n\tInter collisions: {calcInter}", 
      f"\n\tForce iso: {isoSwitch}") 

    # create directory for data 
    dir_name = get_save_name(key, folder, level)
    #print( dir_name )
    if os.path.isdir(dir_name):
      print("\nDirectory: ", os.path.split(dir_name)[1], " already exists...\n\t...Checking all files enclosed")
    else:
      os.mkdir( get_save_name(key, folder, level) )
    
    outputFiles = get_files(folder, include=outputType, get_all=False)

    if n_increments == False:
      contour_save = []
    else:
      if n_increments ==  1:
        contour_save_increment = 1
      else:
       contour_save_increment=math.floor(len(outputFiles)/(n_increments-1))

      for i in range(n_increments-1):
        contour_save.append( int(i*contour_save_increment) )
      contour_save.append(int(len(outputFiles)-1))
      print("Save contour information every", contour_save_increment)
      contour_save = []

    din = []
    counter = 0
    for f in outputFiles:
      if counter in contour_save and n_increments != False:
        saveContour = True
      else: saveContour = False
      din.append({"dataName":f, "level":level, "max_res":max_res, 
        "window":window, "record_contour":saveContour, "key":key, 
        "folder":folder, "dir_name":dir_name, 'braginskiiVorticity':braginskiiVorticity, 
        'bragVortRhoVar':bragVortRhoVar, "interfaceHeuristic":interfaceHeuristic, 
        'calcIntra':calcIntra, calcInter:'calcInter', 'isoSwitch':isoSwitch})
      counter += 1 
    print(f"Begin reading:\t {key}")
    data = []
    if nproc == 1:
        data = []
        for d in din: 
          get_single_data(d)
    else:
        p = Pool(nproc)
        p.map(get_single_data, din)
        p.close()
        p.join()
    
    # calculate growth rate data ## moved into the time series art of the post cde
    """ fix this 
    data[0]["growth_rate"] = 0.
    data[0]["global_growth_rate"] = 0.
    for i in range(1, len(data)):
      data[i]["growth_rate"] = get_growth(data[i]["ion"]["y_avg_int_width"], 
                                 data[i-1]["ion"]["y_avg_int_width"], data[i]["t"], 
                                 data[i-1]["t"])
      data[i]["global_growth_rate"] = get_growth(data[i]["ion"]["global_int_width"], 
                                 data[i-1]["ion"]["global_int_width"], data[i]["t"], 
                                 data[i-1]["t"])
     """
#    save.create_dataset("n_data", data=len(data))
#    for i, d in enumerate(data):
#        h5 = save.create_group("data_%i"%i)
#        save_dict(d, h5)
#    
#    save.close()
    
    print("\nData:\t ", folder," : done")
    del data, din, outputFiles ; gc.collect() 
    return

#def get_1D_data(h5, species="ion", quantity="baroclinic_sum", nproc = 1, cumsum=False):
def get_1D_time_series_data(file_list, species="ion", quantity="baroclinic_sum", nproc = 8, cumsum=False):
    print(f"Accumulating {species} {quantity} data for time series...")
    t_output   = []
    data_output = []
    din = []

    for i in range(len(file_list)):
      din.append({'file':file_list[i], 'species':species, 'quantity':quantity, 'cumsum':cumsum})    
    if nproc == 1:# if single processor
        data = []
        for d in din: 
          data.append(get_1D_time_series_data_slave(d))
    else:# if multiprocessing 
        p = Pool(nproc)
        data = p.map(get_1D_time_series_data_slave, din)
        #data = p.imap(get_single_data, din, )
        p.close()
        p.join()
  
    if len(data) < 2: 
      print("insufficient data for time series")
      xxxx
    if "growth_rate" in quantity:
      if False: #backward difference 
        data_output = []
        t_output = []
        for i in range(len(data)):
          t_output.append(data[i][0])
          if i == 0:
            data_output.append(0.)
          else:
            data_output.append(get_growth(data[i][1], data[i-1][1], t_output[i], t_output[i-1] ))

      else: # central difference 
        t_output = []
        for i in range(len(data)):
          t_output.append(data[i][0])
          data_output.append(data[i][1][0])
        data_output = np.gradient( np.array(data_output), np.array(t_output))

    else: 
      for i in range(len(data)):
        t_output.append(data[i][0])
        data_output.append(data[i][1])

    del data; gc.collect()

    print(f"\t...finished getting {species} {quantity} data.")
    return t_output, data_output

#def get_1D_data_slave(input_data):
def get_1D_time_series_data_slave(input_data):
  
  h5 = input_data['file']
  try:
    h5 = h5py.File( h5, "r")
  except:
    print("Read failed\n")
    print("File: \t", input_data['file'])
  species = input_data['species']
  quantity = input_data['quantity']
  cumsum = input_data['cumsum']     
  t = h5["t"][()]
  path = "%s/%s"%(species, quantity)
  if quantity == 'interface_location':
    int_loc = {}
    for j in range(len(h5[path])):
      int_loc[j]= h5[path+'/%i'%j][()]
    data = int_loc
    del int_loc
    gc.collect()
  elif quantity == "growth_rate":
    data = h5[species+'/'+"y_avg_int_width"][()]
  elif quantity == "global_growth_rate":
    data = h5[species + '/' + "global_int_width"][()]
  elif path in h5:
    data = h5[path][()]
  else:
    pdb.set_trace()
  h5.close(); del h5; gc.collect()
  return [t, data]

def get_2D_data(h5, species="ion", quantity="vorticity"):
    #print( t_integer)
    #n = h5["n_data"][()]
    #if t_integer >= n: 
      #print("you ducked up", t_integer, n)
      #xxxx
    #i = t_integer # n-1
    #t = h5["data_%i/t"%i][()]
    #x = h5["data_%i/x"%i][()]
    #y = h5["data_%i/y"%i][()]
    #path = "data_%i/%s/%s"%(i, species, quantity)

    t = h5["t"][()]
    x = h5["x"][()]
    y = h5["y"][()]
    path = "%s/%s"%(species, quantity)

    if path in h5:
        dat = h5[path][()]
    else:
        raise RuntimeError(path + " does not exist") 
    
    return t, x, y, dat

def get_figure_size(width, nR, nC, rSP, wspacing, hspacing, fraction=1):
  #width = 6.7 #inches
  """ Set aesthetic figure dimensions to avoid scaling in latex.
  Parameters
  ----------
  width: float - Width of the page in inches
  nR - number of rows
  nC - number of columns 
  rSP - ratio of the subplots assumed homogeneous, height to width
  spacing - spacing between plots as for gridspec
  fraction: float - Fraction of the width which you wish the figure to occupy
  Returns
  ------
  fig_dim: tuple
          Dimensions of figure in inches
  """
  ratio_figure = rSP*(nR + (nR-1)*hspacing*2)/(nC + (nC-1)*wspacing)
  fig_width = width * fraction # Width of figure
  golden_ratio = (5**.5 - 1) / 2 # Golden ratio to set aesthetic figure height
  fig_height = fig_width * ratio_figure # Figure height in inches
  return (fig_width, fig_height)

def get_grad_components(input_array, input_dx, input_dy=False):
  if input_dy == False: input_dy = input_dx
  """
  Find the pressure gradient components. Return the x and y gradients. """
  input_array_dx = np.gradient(input_array, input_dx, axis=0)
  input_array_dy = np.gradient(input_array, input_dy, axis=1)

  return input_array_dx, input_array_dy 

def get_Lorentz(ch, options):
  name = options['name']
  quantity = options['quantity']
  level = options['level']
  #print(quantity)
  """Calculate the Lorentz force and Lorenz force torque contributions to 
  vorticity and individual contributions and return a specfific contribution 
  or total  
  Forces and contributions:
    L_x_E; L_y_E;  
    L_x_B; L_y_B;
    L_x_total; L_y_total
  Torques and contributions:
    tau_E; tau_B, tau_L - total Lorentz force torque
  """
  x, rho = ch.get("rho-%s"%name)
  try:
      u = ch.get("x_vel-%s"%name)[1]
      v = ch.get("y_vel-%s"%name)[1]
      w = ch.get("z_vel-%s"%name)[1]
  except:
      mx = ch.get_flat("x_mom-%s"%name)[1]
      my = ch.get_flat("y_mom-%s"%name)[1]
      mz = ch.get_flat("z_mom-%s"%name)[1]
      u = mx/rho; v = my/rho; w = mz/rho;
  y = [1]; x = x[0]; dx = ch.data["levels"][level]['dx']; dy = dx[1]; dx = dx[0];

  Ex = ch.get("x_D-field")[1]; Ey = ch.get("y_D-field")[1]; 
  Ez = ch.get("z_D-field")[1]
  
  Bx = ch.get("x_B-field")[1]; By = ch.get("y_B-field")[1]
  Bz = ch.get("z_B-field")[1]
  mass = ch.get("mass-"+name)[1]; charge = ch.get("charge-"+name)[1]
  beta = ch.data["beta"]; dS = ch.data["skin_depth"]; c = ch.data["lightspeed"]

  # find Lorentz torque component from electric and magnetic field contributions
  if quantity == 'L_x_E':
    return x, y, np.sqrt(2/beta)*charge/(mass*dS)*c*Ex
  elif quantity == 'L_y_E':
    return x, y, np.sqrt(2/beta)*charge/(mass*dS)*c*Ey 
  elif quantity == 'L_x_B':
    return x, y, np.sqrt(2/beta)*charge/(mass*dS)*(v*Bz - w*By)
  elif quantity == 'L_y_B': 
    return x, y, np.sqrt(2/beta)*charge/(mass*dS)*(w*Bx - u*Bz)
  elif quantity == 'L_x_total':
   return x, y, np.sqrt(2/beta)*charge/(mass*dS)*c*Ex +\
          np.sqrt(2/beta)*charge/(mass*dS)*(v*Bz - w*By)
  elif quantity == 'L_y_total':
   return x, y, np.sqrt(2/beta)*charge/(mass*dS)*c*Ey + \
          np.sqrt(2/beta)*charge/(mass*dS)*(w*Bx - u*Bz)
 
  else:
    Lx_E = np.sqrt(2/beta)*charge/(mass*dS)*c*Ex
    Ly_E = np.sqrt(2/beta)*charge/(mass*dS)*c*Ey 
    Lx_B = np.sqrt(2/beta)*charge/(mass*dS)*(v*Bz - w*By)
    Ly_B = np.sqrt(2/beta)*charge/(mass*dS)*(w*Bx - u*Bz)

  dx = x[1] - x[0]; dy = y[1] - y[0];
  
  dLx_E_dy = np.gradient(Lx_E, dy, axis=1)
  dLy_E_dx = np.gradient(Ly_E, dx, axis=0)
  dLx_B_dy = np.gradient(Lx_B, dy, axis=1)
  dLy_B_dx = np.gradient(Ly_B, dx, axis=0)
  Lx = Lx_E + Lx_B
  Ly = Ly_E + Ly_B

  dLx_dy = np.gradient(Lx, dy, axis=1)
  dLy_dx = np.gradient(Ly, dx, axis=0)
  
  if quantity == 'tau_E':
    return x, y, dLy_E_dx - dLx_E_dy
  elif quantity == 'tau_B':
    return x, y, dLy_B_dx - dLx_B_dy
  elif quantity == 'tau_L':
    return x, y, dLy_dx - dLx_dy
  else: 
    print("error in get_lorentz oooof")
    pdb.set_trace()
  #tau_E = dLy_E_dx - dLx_E_dy
  #tau_B = dLy_B_dx - dLx_B_dy
  #tau = dLy_dx - dLx_dy
  
  #return x, y, tau_E, tau_B, tau   

def get_vorticity_dot(ch, options):
  name = options['name']; quantity = options['quantity']; level = options['level'];
  """Calculate the terms of the vorticity equations and return  
    vort_dot_total
    vort_dot_comp
    vort_dot_baro
    vort_dot_LorE
    vort_dot_LorB
  """
  x, rho = ch.get("rho-%s"%name)
  try:
      u = ch.get("x_vel-%s"%name)[1]
      v = ch.get("y_vel-%s"%name)[1]; w = ch.get("z_vel-%s"%name)[1]
  except:
      x, mx = ch.get("x-mom-%s"%name)
      x, my = ch.get("y-mom-%s"%name)
      x, mz = ch.get("z-mom-%s"%name)
      u = mx/rho; v = my/rho; w = mz/rho;
  dx = ch.data['levels'][level]['dx'][0]; dy = ch.data['levels'][level]['dx'][1]
 
  ### Comp/Conv
  if quantity == 'vort_dot_total' or quantity == 'vort_dot_comp':
    du_dy = np.gradient(u, dy, axis=1); dv_dx = np.gradient(v, dx, axis=0)
    du_dx = np.gradient(u, dx, axis=0); dv_dy = np.gradient(v, dy, axis=1)
  
    omega = dv_dx - du_dy
  
    d_omega_dx = np.gradient(omega, dx, axis=0); 
    d_omega_dy = np.gradient(omega, dy, axis=1)

    vort_dot_comp = -omega*(du_dx + dv_dy); 
    #vort_dot_conv = -u*d_omega_dx - v*d_omega_dy 
    if quantity == 'vort_dot_comp':
      return x[0], x[1], vort_dot_comp
  ### Baro
  if quantity == 'vort_dot_total' or quantity == 'vort_dot_baro':
    try:
      x, prs = ch.get("p-%s"%name)
    except:
      x, nrg = ch.get("rho_E-%s"%name)
      x, gam = ch.get("gamma-%s"%name)
      prs = (nrg - (mx**2 + my**2 + mz**2)/(2*rho))*(gam - 1.0)

    rho_dx = np.gradient(rho, dx, axis=0); rho_dy = np.gradient(rho, dy, axis=1)
    prs_dx = np.gradient(prs, dx, axis=0); prs_dy = np.gradient(prs, dy, axis=1)
  
    vort_dot_baro = (1/rho**2)*(rho_dx*prs_dy - rho_dy*prs_dx)
    if quantity == 'vort_dot_baro':
      return x[0], x[1], vort_dot_baro

  ### LorE and LorB
  if quantity == 'vort_dot_total' or quantity == 'vort_dot_LorE' or \
     quantity == 'Lorentz' or quantity == 'vort_dot_LorB':
    Ex = ch.get("x_D-field")[-1]; Ey = ch.get("y_D-field")[-1]; Ez = ch.get("z_D-field")[-1]
    Bx = ch.get("x_B-field")[-1]; By = ch.get("y_B-field")[-1]; Bz = ch.get("z_B-field")[-1]
    mass = ch.get("mass-"+name)[-1]; charge = ch.get("charge-"+name)[-1]
    beta = ch.data["beta"]; dS = ch.data["skin_depth"]; c = ch.data["lightspeed"]
  
    # find Lorentz torque component from electric and magnetic field contributions
    Lx_E = np.sqrt(2/beta)*charge/(mass*dS)*c*Ex
    Ly_E = np.sqrt(2/beta)*charge/(mass*dS)*c*Ey 
    Lx_B = np.sqrt(2/beta)*charge/(mass*dS)*(v*Bz - w*By)
    Ly_B = np.sqrt(2/beta)*charge/(mass*dS)*(w*Bx - u*Bz)
  
    Lx = np.sqrt(2/beta)*charge/(mass*dS)*(c*Ex + v*Bz - w*By)
    Ly = np.sqrt(2/beta)*charge/(mass*dS)*(c*Ey + w*Bx - u*Bz)
    if quantity == 'Lorentz':
      return x[0], x[1], Lx, Ly

    dLx_E_dy = np.gradient(Lx_E, dy, axis=1); dLy_E_dx = np.gradient(Ly_E, dx, axis=0)
    dLx_B_dy = np.gradient(Lx_B, dy, axis=1); dLy_B_dx = np.gradient(Ly_B, dx, axis=0)

    dLx_dy = np.gradient(Lx, dy, axis=1); dLy_dx = np.gradient(Ly, dx, axis=0)
  
    tau_E = dLy_E_dx - dLx_E_dy; tau_B = dLy_B_dx - dLx_B_dy
    #tau = dLy_dx - dLx_dy
  
    #if abs(np.sum(tau) - (np.sum(tau_E) + np.sum(tau_B))) > 0.1:
      #print('EM torque do not match')
      #print('np.sum(tau)',np.sum(tau))
      #print('np.sum(tau_E)+...',np.sum(tau_E) + np.sum(tau_B))
      #xxxx
    if quantity == 'vort_dot_LorE':
      return x, y, tau_E
    elif quantity == 'vort_dot_LorB':
      return x, y, tau_B
  #####
  if quantity == 'vort_dot_total':
    return x, y, vort_dot_comp + vort_dot_baro + tau_E + tau_B

def get_rho_E_EM(ch, options):
  level = options['level']

  #"rho_E-EM" in attr_name_access: # special case of a derived property 
  x, Ex = ch.get_flat("x_D-field"); y = x[1]; x = x[0];
  Ey = ch.get_flat("y_D-field")[1]; Ez = ch.get_flat("z_D-field")[1]
  
  Bx = ch.get_flat("x_B-field")[1]; By = ch.get_flat("y_B-field")[1]
  Bz = ch.get_flat("z_B-field")[1]

  beta = ch.data['beta'][1]
  dx = ch.data['levels'][level]['dx'][0]; dy = ch.data['levels'][level]['dx'][1]
  rho_E_EM = (Bx**2 + By**2 + Bz**2 + Ex**2 + Ey**2 + Ez**2)/beta
  energySum = np.sum(rho_E_EM)*dx*dy

  return rho_E_EM, energySum

def get_T_cons(rc, options):
  name = options['name'] 
  x, density = rc.get_flat("rho-%s"%name)
  try:
    x, mx = rc.get_flat("x-mom-%s"%name)
    x, my = rc.get_flat("y-mom-%s"%name)
    x, mz = rc.get_flat("z-mom-%s"%name)
  except:
    mx = rc.get_flat("x_vel-%s"%name)[1]*density
    my = rc.get_flat("y_vel-%s"%name)[1]*density
    mz = rc.get_flat("z_vel-%s"%name)[1]*density
    
  x, nrg = rc.get_flat("rho_E-%s"%name)
  x, gam = rc.get_flat("gamma-%s"%name)    
  prs = (nrg - (mx**2 + my**2 + mz**2)/(2*density))*(gam - 1.0)
  mass = rc.get("mass-%s"%name)[2]
  return prs/(density/mass)

def get_EM(rc, options):
  quantity = options['quantity']
  if quantity not in ['x-electric', 'y-electric', 'z-electric', 'x-magnetic', 'y-magnetic', 'z-magnetic', 'density-charge', 'x-density-current', 'y-density-current']:
    print('incorrect input quantity:\t', quantity)
    return 

  if ('electric' in quantity) or ('magnetic' in quantity):
    return rc.expression('{'+quantity+'}') #rc.get_flat(quantity) #
  else:
    # problem constants
    mass = rc.mass; charge = rc.charge
    density_charge = 0.; x_density_current = 0.; y_density_current = 0.
    for (name,index) in [('electron', 0), ('ion',1)]: 
      m0 = mass[index,0]; m1 = mass[index,1]
      q0 =charge[index,0]; q1 = charge[index,1]
  
      x, y, rho = rc.expression("{rho-%s}"%name)

      if 'current' in quantity:
        try:
          x_vel = rc.get_flat("x_vel-%s}"%name)[2]
          y_vel = rc.get_flat("y_vel-%s}"%name)[2]
          z_vel = rc.get_flat("z_vel-%s}"%name)[2]
        except:
          x, y, mx = rc.expression("{x-mom-%s}"%name)
          x, y, my = rc.expression("{y-mom-%s}"%name)
          x, y, mz = rc.expression("{z-mom-%s}"%name)
        
          x_vel = mx/rho
          y_vel = my/rho
          z_vel = mz/rho

      x, y, tracer  = rc.expression("{tracer-%s}"%name)
      alpha      = np.divide(tracer, rho)
      mass_mix   = (m0*m1)/(m0*alpha + m1*(1.-alpha)) 
      charge_mix = (alpha*m0*q1 + (1.-alpha)*m1*q0)/(m0*alpha + m1*(1.-alpha))
      density_n  = np.divide(rho, mass_mix)
  
      density_charge += np.multiply(charge_mix, density_n)
      if 'x-density' in quantity:
        x_density_current += np.multiply( x_vel, density_charge )
      elif 'y-density' in quantity:
        y_density_current += np.multiply( y_vel, density_charge )

      #del m0, m1, q0, q1, rho, x_vel, y_vel, z_vel, tracer, alpha, mass_mix, charge_mix, density_n; gc.collect()

    if 'charge' in quantity:
      return x, y, density_charge
    elif 'x-' in quantity:
      return x, y, x_density_current
    else:
      return x, y, y_density_current
#3###############################################################333

def get_transportProperties(ch, names, level, calcIntra=True, calcInter=True, isoOveride=False, useNPROC=1):
  """
  function takes as input the data file opened with daryl's post function
  and extract outs primitives, gradients, then calculates the intra and 
  inter-species collisional controbutions.

  inputs:
    ch -  ReadBoxlib object having read data file 
    names - list of strings e.g. ['ions', 'electrons'] for the states
    level - gid depth 
    isoOveride - true or false to foce isotropic 
    useNPROC - integer number of processor to use
    
  outputs:
    x - numoy array of x points (nodes)
    y - numpy array of y points (nodes)
    dudt_fluxes - flux contributions to the cell values (time rate)
    srcDst - source term contributions to the cell evoltions (time rate)
  """

  print(f"\t\tNote isoOveride set to:{isoOveride}")
  Q = {}
  #print(f"Time is:\t{ch.time}")
  Density=0; Xvel = 1; Yvel=2; Zvel=3; Prs=4; Temp=5; Alpha=6;
  x_D = 0; y_D = 1; z_D = 2; x_B = 3; y_B = 4; z_B = 5; muIdx = 6; epIdx = 7;
  #print("\tExtracting primitives...")
  for name in names: # get the ion and electron prims
    Q[name] = {}
    try:
      x, Q[name][Density] = ch.get("rho-%s"%name)
    except:
      x, Q[name][Density] = ch.get_flat("rho-%s"%name)

    try:
        Q[name][Xvel] = ch.get("x_vel-%s"%name)[-1]
        Q[name][Yvel] = ch.get("y_vel-%s"%name)[-1]
        Q[name][Zvel] = ch.get("z_vel-%s"%name)[-1]
    except:
        mx = ch.get_flat("x_mom-%s"%name)[1]
        my = ch.get_flat("y_mom-%s"%name)[1]
        mz = ch.get_flat("z_mom-%s"%name)[1]
        Q[name][Xvel] = mx/rho;
        Q[name][Yvel] = my/rho;
        Q[name][Zvel] = mz/rho;

    Q[name][Prs] = ch.get("p-%s"%name)[-1]
    Q[name][Temp] = ch.get("T-%s"%name)[-1]
    Q[name][Alpha] = ch.get("alpha-%s"%name)[-1]

    Q[name]['mass'] = ch.get("mass-"+name)[-1]; Q[name]['charge'] = ch.get("charge-"+name)[-1]
    Q[name]['gamma'] = ch.get("gamma-"+name)[-1];
    
  y = x[1]; x = x[0]; dx = ch.data["levels"][level]['dx']; dy = dx[1]; dx = dx[0];

  Q["field"] = {}
  Q["field"][x_D] = ch.get("x_D-field")[-1]; Q["field"][y_D] = ch.get("y_D-field")[-1]; 
  Q["field"][z_D] = ch.get("z_D-field")[-1];
  Q["field"][epIdx] = ch.get("ep-field")[-1];
  
  Q["field"][x_B] = ch.get("x_B-field")[-1]; Q["field"][y_B] = ch.get("y_B-field")[-1]
  Q["field"][z_B] = ch.get("z_B-field")[-1]

  beta = ch.data["beta"]; dS = ch.data["skin_depth"]; c = ch.data["lightspeed"]
  n_ref = ch.data["n0"]; x_ref = ch.data["x_ref"]; u_ref = ch.data["u_ref"];
  rho_ref = ch.data["rho_ref"]; T_ref = ch.data["T_ref"]; m_ref = ch.data["m_ref"]; 

  # check if primitives that we have access to the cerberus gradients - super important 
  QD = {} #gradients of primitives 
  properties = [#'x_D-field-dx', 'x_D-field-dy', 
                #'y_D-field-dx', 'y_D-field-dy', 
                #'x_B-field-dx', 'y_B-field-dx', 'y_B-field-dy', 'z_B-field-dx', 'z_B-field-dy',
                #'p-electrons-dx', 'p-electrons-dy', 'p-ions-dx', 'p-ions-dy', 
                'x_vel-ions-dx', 'y_vel-ions-dx', 'x_vel-ions-dy', 'y_vel-ions-dy', 
                'x_vel-electrons-dx', 'y_vel-electrons-dx', 
                'x_vel-electrons-dy', 'y_vel-electrons-dy', 
                'T-electrons-dx', 'T-electrons-dy', 'T-ions-dx', 'T-ions-dy', 
               ]
  propertiesHandle = {'x_D-field-dx':"xDdx", 'x_D-field-dy':"xDdy", 
                      'y_D-field-dx':"yDdx", 'y_D-field-dy':"yDdy", 
                      'x_B-field-dx':"xBdx", 'y_B-field-dx':"yBdx", 
                      'y_B-field-dy':"yBdy", 'z_B-field-dx':"zBdx", 'z_B-field-dy':"zBdy",
                      'p-electrons-dx':"Pedx", 'p-electrons-dy':"Pedy", 
                      'p-ions-dx':"Pidx", 'p-ions-dy':"Pidy", 
                      'x_vel-ions-dx':"Uidx", 'y_vel-ions-dx':"Vidx", 
                      'x_vel-ions-dy':"Uidy", 'y_vel-ions-dy':"Vidy", 
                      'x_vel-electrons-dx':"Uedx", 'y_vel-electrons-dx':"Vedx", 
                      'x_vel-electrons-dy':"Uedy", 'y_vel-electrons-dy':"Vedy", 
                      'T-electrons-dx':"Tedx", 'T-electrons-dy':"Tedy", 
                      'T-ions-dx':"Tidx", 'T-ions-dy':"Tidy"}

  for prop in properties:
    try:
      QD[propertiesHandle[prop]] = ch.get(prop)[-1]
      #print(f"\t\tGradient {prop} read")
    except:
      print(f"\t\t-->>Gradient {prop} unavailable")
  #print("\t...extracted\n")
  ch.close() # clsoe this bith 
  del ch; gc.collect()

  # prepare indivudal component data --- pulled in from 1D code and adapted for 2D
      # needs to be done for every interface - multiprocesses for each row ?
  #print("\tViscousTensor and heat flux calc")
  X = 0; Y=1; Z=2;
  Xmom = 0; Ymom = 1; Zmom = 2; EdenPi = 3; EdenQ = 4

  if useNPROC == 1: print("Warning --- get_transport is written for parallel, 1 cpu is v slow")

  ####-----------------------------------Flux Routine 
  if calcIntra:
    fluxX = {}; fluxY = {}
    #print(x.shape, y.shape)
    for name in names:
      fluxX[name] = np.zeros((x.shape[0] + 1, y.shape[0] + 1, EdenQ+1)); 
      fluxY[name] = np.zeros((x.shape[0] + 1, y.shape[0] + 1, EdenQ+1)) 
      # cell i,j will handle the interface between
      #   fluxX: i-1, j and i, j e.g. the lo flux for corresponding volume cell i,j is fluxX[i,j], 
      #   fluxY: i, j-1 and i, j
      # Ghost cell values are neuman = 0 from interior i.e. zero gradients on boundaries 
  
      #===========================  multiprocessing  ============================
      din = []
      ### decide the partitions in data 
      if useNPROC == 24:
        ny = 4; nx = 6
      else:
        ny = 2; nx = 4; 
       
      dy_cell = int((y.shape[0]+1)/ny) # number of interfaces in x  #TODO 
      dx_cell = int((x.shape[0]+1)/nx) # number of interfaces in y  #TODO
    
      xcell =[i*dx_cell for i in range(nx)]; xcell.append(x.shape[0] + 1) # up to last domain bc we cannot get the right hand state of the last cells boundary
      ycell =[i*dy_cell for i in range(ny)]; ycell.append(y.shape[0] + 1) #TODO
      
      #print(f"Full domain nX: {x.shape[0]}\tnY: {y.shape[0]}")   
      din = []
      for i in range(len(xcell)-1):
        ih = xcell[i+1]
        for j in range(len(ycell)-1):
          jh = ycell[j+1]
          #print(f"Block bounds: {xcell[i]}, {ih}\t{ycell[j]}, {jh}")
          il = xcell[i]
          jl = ycell[j] # il and ih refer to the cells (in prim array that are hanbfdled 
          iRange = [il, ih-1]; jRange = [jl, jh-1];
          #Note python auto makes pointer rferences this isnt copyying the primitives  ---- moight ot be the case for multiprocess  
          din.append({"Q":Q, "QD":QD, 
            "name":name, "ionName":"ions", "eleName":"electrons", "fieldName":"field", 
            "iRange":iRange, "jRange":jRange, "Debye":dS/c, "Larmor":math.sqrt(beta/2)*dS, 
            "lightspeed":c, "xref":x_ref, "n0ref":n_ref, "mref":m_ref, "rhoref":rho_ref, 
            "Tref":T_ref, "uref":u_ref, "dx":dx, "isoOveride":isoOveride, "verbosity":1, 
            'xhiDomain':x.shape[0],'yhiDomain':y.shape[0]}) #xhiDomain etc passed as size bc flux index at ends is last index of Q+1
  
      data = []
    
      nproc = useNPROC
      if nproc == 1:
          for d in din: 
            data.append(braginskiiViscousWrapper(d))
      else:
            p = Pool(nproc)
            data = p.map(braginskiiViscousWrapper, din)
            p.close()
            p.join()
            # stitch back together
    
      counter = 0
      for i in range(len(xcell)-1):
        ih = xcell[i+1]
        for j in range(len(ycell)-1):
          jh = ycell[j+1]
          fluxX[name][xcell[i]:ih, ycell[j]:jh, :] = data[counter][0]
          fluxY[name][xcell[i]:ih, ycell[j]:jh, :] = data[counter][1]
          counter += 1
    # find max in each region and time 
      # exact flux values 
    dudt_flux = {}
    dt = 1 #Nontrivial time step for contrbution instead of the rate of change 
    #print("\tCalculating viscous flux contribution...")
    for name in names:
      # old for loop method
      dudt_flux[name] = np.zeros((x.shape[0], y.shape[0], EdenQ+1)); 
      """
      for prop in range(EdenQ+1):
        for j in range(1, y.shape[0]-1):#TODO replace with multiprocessing 
          for i in range(1, x.shape[0]-1):#TODO replace with multiprocessing 
            #dudt_flux[name][i,j, prop] = dt/dx * ( fluxX[i,j,prop] - fluxX[i-1,j,prop] + \
            #                           fluxY[i,j,prop] - fluxY[i,j-1,prop])
            
            dudt_flux[name][i,j, prop] = dt/dx * ( fluxX[i-1,j,prop] - fluxX[i,j,prop] + \
                                       fluxY[i,j-1,prop] - fluxY[i,j,prop])
      """
      #dUfluxX = dt/dx*(fluxX[name][0:-1, 1:  , :] - fluxX[name][1:, 1:, :])
      dUfluxX = dt/dx*(fluxX[name][0:-1, :-1  , :] - fluxX[name][1:, :-1, :])
      #dUfluxY = dt/dx*(fluxX[name][1:,   0:-1, :] - fluxX[name][1:, 1:, :])
      dUfluxY = dt/dx*(fluxY[name][:-1,   0:-1, :] - fluxY[name][:-1, 1:, :])
      dudt_flux[name][:,:, :] = dUfluxX + dUfluxY

    if calcInter==False: return x, y, dudt_flux, None
  #print("\t\t..Calc done")
  
  # source term contributions 
  if calcInter:
    #print("\tCalculating src term contribution...neglect the boardering cells")
    # here we exlcude the boarder as to match the flux registers controbution to all the interior
    # cells but not the boarder. 
    srcDst = {"electrons":np.zeros((x.shape[0], y.shape[0], EdenQ+1)), 
              "ions":np.zeros((x.shape[0], y.shape[0], EdenQ+1))} #TODO
  
    #break domain into chunks 
    #print(f"i and j components of domain are: {x.shape[0]} {y.shape[0]}")
    if useNPROC == 24:
      ny = 4; nx = 6
    else:
      ny = 2; nx = 4; 

    dy_cell = int((y.shape[0])/ny) #TODO
    dx_cell = int((x.shape[0])/nx) #TODO
  
    xcell = [i*dx_cell for i in range(nx)]; xcell.append(x.shape[0])#TODO
    ycell = [i*dy_cell for i in range(ny)]; ycell.append(y.shape[0])#TODO
    
    din = []
    for i in range(len(xcell)-1):
      ih = xcell[i+1]
      for j in range(len(ycell)-1):
        jh = ycell[j+1]
        Qin = {}; QDin = {}
        #print(f"Block bounds: {xcell[i]}, {ih}\t{ycell[j]}, {jh}")
        for key in Q.keys():
          Qin[key] = {}
          for key2 in Q[key].keys():
            Qin[key][key2] = Q[key][key2][xcell[i]:ih, ycell[j]:jh]
        # for inputs 
        xcells = Qin[key][key2][xcell[i]:ih, ycell[j]:jh].shape[0]
        ycells = Qin[key][key2][xcell[i]:ih, ycell[j]:jh].shape[1]
  
        if len(QD.keys()) != 0: 
          for key in QD.keys():
            QDin[key] = QD[key][xcell[i]:ih, ycell[j]:jh]
  
        #print(f"\txcells: {xcells}\tycells: {ycells}")
        
        din.append({"xcells":xcells, "ycells":ycells, "Q":Qin, "QD":QDin, 
          "ionName":"ions", "eleName":"electrons", "fieldName":"field", 
          "Debye":dS/c, "Larmor":math.sqrt(beta/2)*dS, "lightspeed":c, 
          "xref":x_ref, "n0ref":n_ref, "mref":m_ref, "rhoref":rho_ref, 
          "Tref":T_ref, "uref":u_ref, "dx":dx, "isoOveride":isoOveride, 
          "verbosity":1})
  
    data = []
  
    nproc = useNPROC
    if nproc == 1:
        for d in din: 
          data.append(braginskiiSourceWrapper(d))
    else:
          p = Pool(nproc)
          data = p.map(braginskiiSourceWrapper, din)
          p.close()
          p.join()
          # stitch back together
  
    counter = 0
    for i in range(len(xcell)-1):
      ih = xcell[i+1]
      for j in range(len(ycell)-1):
        jh = ycell[j+1]
        #print(f"\tBlock bounds: {xcell[i]}, {ih}\t{ycell[j]}, {jh}")
  
        srcDst["ions"][xcell[i]:ih, ycell[j]:jh, :] = data[counter]["ions"]
        srcDst["electrons"][xcell[i]:ih, ycell[j]:jh, :] = data[counter]["electrons"]
        counter += 1

    if calcIntra == False: return x, y, None, srcDst

  return x, y, dudt_flux, srcDst

def braginskiiViscousWrapper(din):#name, ionName, eleName, emName, iRange, jRange,
                              #Q, QD, Debye, Larmor, n0_ref, x_ref, u_ref, dX, verbosity = 1):

  name = din["name"]; iRange = din["iRange"]; jRange = din["jRange"];
  Q = din["Q"]; QD = din["QD"]
  isoOveride = din["isoOveride"];
  xhiDomain = din['xhiDomain']; yhiDomain = din['yhiDomain']

  #print("\n\n####Viscous fluxes on flud:\t", name)

  Xmom = 0; Ymom = 1; Zmom = 2; EdenPi = 3; EdenQ = 4
  Density=0; Xvel = 1; Yvel=2; Zvel=3; Prs=4; Temp=5; Alpha=6;
  x_D = 0; y_D = 1; z_D = 2; x_B = 3; y_B = 4; z_B = 5; muIdx = 6; epIdx = 7;
  nFluxesX = iRange[1] - iRange[0] + 1; nFluxesY = jRange[1] - jRange[0] + 1; 
  fluxX = np.zeros((nFluxesX, nFluxesY, EdenQ+1)); 
  fluxY = np.zeros((nFluxesX, nFluxesY, EdenQ+1))
  # note flux register loops below, not cell
  for i in range(iRange[0], iRange[1]+1): # swapped i and j orders -  shouldn't make a difference but #TODO check  
    for j in range(jRange[0], jRange[1]+1):# 
      #mapping from global domain to local block indexing 
      if i == xhiDomain and j == yhiDomain:
        #print("Corner skipped ", i,j)
        continue

      iloc = i - iRange[0]; jloc = j - jRange[0]

      # fluxes in the x-direcion
      if (j != yhiDomain):
        ViscTens, q_flux = \
          braginskiiViscousTensorHeatFlux(name, "ions", "electrons", "field", i, j, 0, 1,
          xhiDomain, yhiDomain,
          din["Q"], din["QD"], din["Debye"], din["Larmor"], din["n0ref"], din["xref"], 
          din["uref"], din["dx"], isoOveride, din["verbosity"])
          # Q, QD, dS/c, math.sqrt(beta/2)*dS, n_ref, x_ref, u_ref, dx, verbosity = 1)
  
        ##  if dimFlux == 0: # in the x-dimension
        ixl = i - 1; ixh = i ;
        iyl = j; iyh = j;
        #special index for ghost cell /interior matching 
        if i == 0: #xlo interface of ghost and first cell of domain or xhi interface of last interior and ghost 
          ixl = i; ixh = i; 
        elif  i == xhiDomain:
          ixl = i-1; ixh=i-1
  
        if j == 0:
          iyl = j; iyh = j; 
        elif j == yhiDomain: 
          iyl = j-1; iyh = j-1; 
  
        fluxX[iloc,jloc,Xmom] += ViscTens[0];
        fluxX[iloc,jloc,Ymom] += ViscTens[3];
        fluxX[iloc,jloc,Zmom] += ViscTens[5];
        fluxX[iloc,jloc,EdenPi] += 0.5*((Q[name][Xvel][ixh,iyh] + Q[name][Xvel][ixl,iyl])*ViscTens[0]+
                           (Q[name][Yvel][ixh,iyh] + Q[name][Yvel][ixl,iyl])*ViscTens[3]+
                           (Q[name][Zvel][ixh,iyh] + Q[name][Zvel][ixl,iyl])*ViscTens[5])
        fluxX[iloc,jloc,EdenQ] += q_flux[0];
  
      # fluxes in the y-direcion
      if (i != xhiDomain):
        ViscTens, q_flux = \
          braginskiiViscousTensorHeatFlux(name, "ions", "electrons", "field", i, j, 1, 1, 
          xhiDomain, yhiDomain,
          din["Q"], din["QD"], din["Debye"], din["Larmor"], din["n0ref"], din["xref"], 
          din["uref"], din["dx"], din["verbosity"])
            #Q, QD, dS/c, math.sqrt(beta/2)*dS, n_ref, x_ref, u_ref, dx, verbosity = 1)  
  
        #else: # in the y dimension 
        ixl = i; ixh = i; 
        iyl = j - 1; iyh = j;
        if i == 0: #xlo interface of ghost and first cell of domain or xhi interface of last interior and ghost 
          ixl = i; ixh = i; 
        elif i == xhiDomain:
          ixl = i-1; ixh = i-1; 
  
        if j == 0: 
          iyl = j; iyh = j; 
        elif j == yhiDomain:
          iyl = j-1; iyh = j-1; 
  
        fluxY[iloc,jloc,Xmom] += ViscTens[3];
        fluxY[iloc,jloc,Ymom] += ViscTens[1];
        fluxY[iloc,jloc,Zmom] += ViscTens[4];
        fluxY[iloc,jloc,EdenPi] += +0.5*((Q[name][Xvel][ixh,iyh]+Q[name][Xvel][ixl,iyl])*ViscTens[3]+
                (Q[name][Yvel][ixh,iyh]+Q[name][Yvel][ixl,iyl])*ViscTens[1]+
                (Q[name][Zvel][ixh,iyh]+Q[name][Zvel][ixl,iyl])*ViscTens[4]) 
  
        fluxY[iloc,jloc,EdenQ] += q_flux[1];
  return fluxX, fluxY

def braginskiiViscousTensorHeatFlux(name, ionName, eleName, emName, i, j, dimFlux, dim,
                              xhiDomain, yhiDomain, 
                              Q, QD, Debye, Larmor, n0_ref, x_ref, u_ref, dX, 
                              isoOveride = False, verbosity = 1):
  Density=0; Xvel = 1; Yvel=2; Zvel=3; Prs=4; Temp=5; Alpha=6;
  x_D = 0; y_D = 1; z_D = 2; x_B = 3; y_B = 4; z_B = 5; muIdx = 6; epIdx = 7;
  #print("Viscous fluxes on flud:\t", name)
  #print("\t", i, j)
  ##mak the interface index play nice with the volume cell index 
  if dimFlux == 0: # in the x-dimension
    xl = i - 1; xh = i;
    yl = j; yh = j
  elif dimFlux ==1 : # in the y dimension 
    xl = i; xh = i; 
    yl = j - 1; yh = j

  if i == 0 : #xlo interface of ghost and first cell of domain or xhi interface of last interior and ghost 
    xl = i; xh = i; # set the ghost cell value to the interior cell value 
  elif i == xhiDomain:
    xl = i-1; xh = i-1; # set the ghost cell value to the interior cell value 

  if j == 0: 
    yl = j; yh = j; 
  elif j == yhiDomain:
    yl = j-1; yh = j-1; 

  ###TODO if we are on the boundary we are using neuman for the fluxes 
  ViscTens = np.zeros((7)); 

  EffectiveZero = 1e-14;
  #collect relevant coefficent data 
    #mag fields on left and right 
  BxL = Q["field"][x_B][xl, yl]; ByL = Q["field"][y_B][xl, yl]; 
  BzL = Q["field"][z_B][xl, yl]
  BL = BxL*BxL + ByL*ByL + BzL*BzL;
  BL_xyz = np.array([BxL, ByL, BzL]);

  BxR = Q["field"][x_B][xh, yh]; ByR = Q["field"][y_B][xh, yh]; BzR = Q["field"][z_B][xh, yh]
  BR = BxR*BxR + ByR*ByR + BzR*BzR;
  BR_xyz = np.array([BxR, ByR, BzR]);

  Bx = 0.5*(BxL + BxR); By = 0.5*(ByL + ByR); Bz = 0.5*(BzL + BzR);
  B = Bx*Bx + By*By + Bz*Bz;
  B_xyz = np.array([Bx, By, Bz]);

  """ just change the actual coefficients applied and the gradients used 
  if isoOveride:
    BL = 0; BL_xyz = np.zeros(3)
    BR = 0; BR_xyz = np.zeros(3)
    Bx = 0; By = 0; Bz = 0.;
    B = 0.0; B_xyz = np.zeros(3)
  """
  if (B < 0.0):
    print("MFP_braginskii.cpp ln 1103 - Negative magnetic field error");
    sys.ext("ln 204 MFP_braginskii.cpp");
  else:
      if (B < EffectiveZero):
        #print("\nzero magnetic field\n");
        B_pp = 0.;
        B_p  = 0.;
      elif ((abs(Bx) < EffectiveZero) and (abs(By) < EffectiveZero) and 
            (abs(Bz) > EffectiveZero)):
        #print("\nzero x and y magnetic field\n");
        B_pp = 1/np.sqrt(B); # B prime prime 
        B_p  = 0.;
      else :  
        B_pp = 1/np.sqrt(B); # B prime prime 
        B_p  = 1/np.sqrt(Bx*Bx + By*By); # B prime 

  bx_pp = Bx*B_pp; bx_p = Bx*B_p;
  by_pp = By*B_pp; by_p = By*B_p;
  bz_pp = Bz*B_pp;
  
  B_unit = np.array([bx_pp, by_pp, bz_pp]); 
    # relative velocity
  u_para = np.zeros((3))  # Velocity parallel to B_unit
  u_perp = np.zeros((3))  # Velocity perpendicular to B_unit
  u_chev = np.zeros((3))  # unit vector perp to u and B_unit
  TG_para = np.zeros((3)) #Temp grad parallel to B_unit 
  TG_perp = np.zeros((3)) #Temp grad perpendicular to B_unit
  TG_chev = np.zeros((3)) #unit vector perp to gradT and B_unit

  # electron primitives 
  u_ele = 0.5*(Q[eleName][Xvel][xl,yl] + Q[eleName][Xvel][xh,yh]);
  v_ele = 0.5*(Q[eleName][Yvel][xl,yl] + Q[eleName][Yvel][xh,yh]);
  w_ele = 0.5*(Q[eleName][Zvel][xl,yl] + Q[eleName][Zvel][xh,yh]);

  # ion primitives
  u_ion = 0.5*(Q[ionName][Xvel][xl,yl] + Q[ionName][Xvel][xh,yh]);
  v_ion = 0.5*(Q[ionName][Yvel][xl,yl] + Q[ionName][Yvel][xh,yh]);
  w_ion = 0.5*(Q[ionName][Zvel][xl,yl] + Q[ionName][Zvel][xh,yh]);

  #pull out the gradient values from the slopes stuff
  if (name == ionName):
    T_i_L, eta0_L, eta1_L, eta2_L, eta3_L, eta4_L, kappa1_L, kappa2_L, kappa3_L, = \
      getIonCoeffs(ionName, eleName, emName, Q, xl, yl, BL_xyz, Debye, Larmor, n0_ref,
                   x_ref, u_ref, verbosity);

    T_i_R, eta0_R, eta1_R, eta2_R, eta3_R, eta4_R, kappa1_R, kappa2_R, kappa3_R, = \
      getIonCoeffs(ionName, eleName, emName, Q, xh, yh, BL_xyz, Debye, Larmor, n0_ref, 
                   x_ref, u_ref, verbosity);

    T_i = 0.5*(T_i_R + T_i_L);
    eta0 = 0.5*(eta0_R + eta0_L);
    eta1 = 0.5*(eta1_R + eta1_L);
    eta2 = 0.5*(eta2_R + eta2_L);
    eta3 = 0.5*(eta3_R + eta3_L);
    eta4 = 0.5*(eta4_R + eta4_L);
    kappa1 = 0.5*(kappa1_R + kappa1_L);
    kappa2 = 0.5*(kappa2_R + kappa2_L);
    kappa3 = 0.5*(kappa3_R + kappa3_L);
    if (verbosity > 4):
      print("Ions") 
    """
      print(f"\n\nIon\neta0\t{eta0}\neta1\t{eta1}\neta2\t{eta2}\neta3\t{eta3}\neta4\t{eta4}")
      print(f"kappa1\t{kappa1}\nkappa2\t{kappa2}\nkappa3\t{kappa3}")
      print(f"beta1\t-\nbeta2\t-\nbeta3\t-")
    """
  elif (name==eleName):
    T_e_L, eta0_L, eta1_L, eta2_L, eta3_L, eta4_L, kappa1_L, kappa2_L, kappa3_L, \
    beta1_L, beta2_L, beta3_L = \
      getElectronCoeffs(ionName, eleName, emName, Q, xl, yl, BL_xyz, Debye, Larmor, n0_ref,
                        x_ref, u_ref, verbosity);
    T_e_R, eta0_R, eta1_R, eta2_R, eta3_R, eta4_R, kappa1_R, kappa2_R, kappa3_R, \
    beta1_R, beta2_R, beta3_R = \
      getElectronCoeffs(ionName, eleName, emName, Q, xl, yl, BL_xyz, Debye, Larmor, n0_ref, 
                        x_ref, u_ref, verbosity);

    T_e = 0.5*(T_e_R + T_e_L);
    eta0 = 0.5*(eta0_R + eta0_L);
    eta1 = 0.5*(eta1_R + eta1_L);
    eta2 = 0.5*(eta2_R + eta2_L);
    eta3 = 0.5*(eta3_R + eta3_L);
    eta4 = 0.5*(eta4_R + eta4_L);
    kappa1 = 0.5*(kappa1_R + kappa1_L);
    kappa2 = 0.5*(kappa2_R + kappa2_L);
    kappa3 = 0.5*(kappa3_R + kappa3_L);
    beta1 = 0.5*(beta1_R + beta1_L);
    beta2 = 0.5*(beta2_R + beta2_L);
    beta3 = 0.5*(beta3_R + beta3_L);

  else:
    sys.exit("non ion/ele rogue state.");

  if name == eleName: species = "e"
  else: species = "i"

  if dimFlux == 0:
    try: # if we have accees to the gradient info 
      dT_dx = 0.5*(QD["T" + species + "dx"][xl,yl]  + QD["T" + species + "dx"][xh,yh])
    except: 
      dT_dx = (Q[name][Temp][xh,yh] - Q[name][Temp][xl,yl])/dX;

    try: # if we have accees to the gradient info 
      du_dx = 0.5*(QD["U" + species + "dx"][xl,yl]  + QD["U" + species + "dx"][xh,yh])
      dv_dx = 0.5*(QD["V" + species + "dx"][xl,yl]  + QD["V" + species + "dx"][xh,yh])
    except:
      du_dx = (Q[name][Xvel][xh,yh] - Q[name][Xvel][xl,yl])/dX; 
      dv_dx = (Q[name][Yvel][xh,yh] - Q[name][Yvel][xl,yl])/dX; 
    dw_dx = (Q[name][Zvel][xh,yh] - Q[name][Zvel][xl,yl])/dX; 

    drho_dz=0; du_dz=0; dv_dz=0; dw_dz=0; dp_dz=0; dT_dz=0; # two d code results only

    if dim == 0: # 1D
      drho_dy=0; du_dy=0; dv_dy=0; dw_dy=0; dp_dy=0; dT_dy=0; 
    elif dim == 1: # 2D # find the derivatives in the dimensions off the flux dimension -  hard coded for x direction flux 
      # if gradeints aren't avaialble from inputs then we need a rectangle stencile on the itnerface
      if j == 0:
        yyl = 0; yyh = j+1; 
      elif j == yhiDomain-1: 
        yyl = j-1; yyh = j; 
      else:
        yyl = j-1; yyh = j+1; 

      try:
        dT_dy = 0.5*(QD["T" + species + "dy"][xl,yl]  + QD["T" + species + "dy"][xh,yh]) 
      except:
        dT_dy = (Q[name][Temp][xh,yyh]+Q[name][Temp][xl,yyh]-Q[name][Temp][xh,yyl]-Q[name][Temp][xl,yyl])/4/dX;

      try:
        du_dy = 0.5*(QD["U" + species + "dy"][xl,yl]  + QD["U" + species + "dy"][xh,yh]) 
        dv_dy = 0.5*(QD["V" + species + "dy"][xl,yl]  + QD["V" + species + "dy"][xh,yh]) 
      except:
        du_dy = (Q[name][Xvel][xh,yyh]+Q[name][Xvel][xl,yyh]-Q[name][Xvel][xh,yyl]-\
                Q[name][Xvel][xl,yyl])/4/dX;
        dv_dy = (Q[name][Yvel][xh,yyh]+Q[name][Yvel][xl,yyh]-Q[name][Yvel][xh,yyl]-\
                Q[name][Yvel][xl,yyl])/4/dX;
      #print(f"yyh: {yyh}\tyyl: {yyl}")
      dw_dy = (Q[name][Zvel][xh,yyl]+Q[name][Zvel][xl,yyh]-Q[name][Zvel][xh,yyl]-\
              Q[name][Zvel][xl,yyl])/4/dX;

  if dimFlux ==1:
    # if gradeints aren't avaialble from inputs then we need a rectangle stencile on the itnerface
    if i == 0:
      xxl = 0; xxh = i+1; 
    elif i == xhiDomain-1: 
      xxl = i-1; xxh = i; 
    else:
      xxl = i-1; xxh = i+1; 

    try:
      dT_dx = 0.5*(QD["T" + species + "dx"][xl,yl]  + QD["T" + species + "dx"][xh,yh]) 
    except:
      dT_dx = (Q[name][Temp][xxh,yh] + Q[name][Temp][xxh,yl] - Q[name][Temp][xxl,yh] - Q[name][Temp][xxl,yl])/4/dX;

    try:
      du_dx = 0.5*(QD["U" + species + "dx"][xl,yl]  + QD["U" + species + "dx"][xh,yh]) 
      dv_dx = 0.5*(QD["V" + species + "dx"][xl,yl]  + QD["V" + species + "dx"][xh,yh])
    except:
      du_dx = (Q[name][Xvel][xxh,yh] + Q[name][Xvel][xxh,yl] -\
              Q[name][Xvel][xxl,yh] - Q[name][Xvel][xxl,yl])/4/dX;
      dv_dx = (Q[name][Yvel][xxh,yh] + Q[name][Yvel][xxh,yl] -\
              Q[name][Yvel][xxl,yh] - Q[name][Yvel][xxl,yl])/4/dX;
    dw_dx = (Q[name][Zvel][xxh,yh] + Q[name][Zvel][xxh,yl] -\
            Q[name][Zvel][xxl,yh] - Q[name][Zvel][xxl,yl])/4/dX;
    try:
      dT_dy = 0.5*(QD["T" + species + "dy"][xl,yl]  + QD["T" + species + "dy"][xh,yh]) 
    except:
      dT_dy = (Q[name][Temp][xh,yh] - Q[name][Temp][xl,yl])/dX;

    try:
      du_dy = 0.5*(QD["U" + species + "dy"][xl,yl]  + QD["U" + species + "dy"][xh,yh]) 
      dv_dy = 0.5*(QD["V" + species + "dy"][xl,yl]  + QD["V" + species + "dy"][xh,yh]) 
    except:
      du_dy = (Q[name][Xvel][xh,yh] - Q[name][Xvel][xl,yl])/dX; 
      dv_dy = (Q[name][Yvel][xh,yh] - Q[name][Yvel][xl,yl])/dX; 
    dw_dy = (Q[name][Zvel][xh,yh] - Q[name][Zvel][xl,yl])/dX; 

    drho_dz=0; du_dz=0; dv_dz=0; dw_dz=0; dp_dz=0; dT_dz=0;

  divu  = du_dx + dv_dy + dw_dz;

  if (verbosity > 4):
    print("divu\t", divu)
  duVec = np.array([u_ele - u_ion, v_ele - v_ion, w_ele - w_ion]); #dU2 = du*du + dv*dv + dw*dw;
  dT_dVec = np.array([dT_dx, dT_dy, dT_dz]);
  #Braginskii directionality stuff. Real dot_B_unit_TG, dot_B_unit_U ;#temp variables
  dot_B_unit_U = bx_pp*duVec[0]+ by_pp*duVec[1]+ bz_pp*duVec[2];
  dot_B_unit_TG = bx_pp*dT_dVec[0] + by_pp*dT_dVec[1] + bz_pp*dT_dVec[2];
  for i_disp in range(3):#i_disp - i disposable 
    u_para[i_disp] = B_unit[i_disp]*dot_B_unit_U ;
    TG_para[i_disp]= B_unit[i_disp]*dot_B_unit_TG ;
    u_perp[i_disp] = duVec[i_disp] - u_para[i_disp];
    TG_perp[i_disp]= dT_dVec[i_disp] - TG_para[i_disp];  
    if (i_disp==0):
      u_chev[i_disp] = B_unit[1]*duVec[2]-B_unit[2]*duVec[1];
      TG_chev[i_disp]= B_unit[1]*dT_dVec[2]-B_unit[2]*dT_dVec[1];
    elif (i_disp==1):
      u_chev[i_disp] = -(B_unit[0]*duVec[2]-B_unit[2]*duVec[0]);
      TG_chev[i_disp]= -(B_unit[0]*dT_dVec[2]-B_unit[2]*dT_dVec[0]);
    else:
      u_chev[i_disp] = B_unit[0]*duVec[1]-B_unit[1]*duVec[0];
      TG_chev[i_disp]= B_unit[0]*dT_dVec[1]-B_unit[1]*dT_dVec[0];

  # Viscous stress tensor 
  Trans = np.zeros((3,3));
  Strain = np.zeros((3,3)); 
  StrainTrans = np.zeros((3,3));
  ViscStress = np.zeros((3,3)); 
  ViscStressTrans = np.zeros((3,3));
  TransT = np.zeros((3,3));
  WorkingMatrix = np.zeros((3,3));
  #Populate the transformation matrix from cartesian normal to B unit
  # aligned cartesian - Li 2018
  Trans[0,0]=-by_p; Trans[0,1]=-bx_p*bz_pp; Trans[0,2]=bx_pp;

  Trans[1,0]= bx_p; Trans[1,1]=-by_p*bz_pp; Trans[1,2] = by_pp;

  Trans[2,0]= 0;    Trans[2,1]= bx_p*bx_pp\
                                +by_p*by_pp; Trans[2,2] =bz_pp;
  #Populate the transpose of the transformation matrix
  TransT[0,0]=Trans[0,0]; TransT[0,1]=Trans[1,0]; TransT[0,2]=Trans[2,0];
  TransT[1,0]=Trans[0,1]; TransT[1,1]=Trans[1,1]; TransT[1,2]=Trans[2,1];
  TransT[2,0]=Trans[0,2]; TransT[2,1]=Trans[1,2]; TransT[2,2]=Trans[2,2];
  """
  if verbosity > 1:
    print("Check transformation matrix is correct");
  """
  if not (TransT != np.transpose(Trans)).any(): pass 
  else:
    print("transformation matrix failed"); pdb.set_trace();
  #Populate strain rate tensor in B unit aligned cartesian frame
  Strain[0,0] = 2*du_dx - 2./3.*divu;
  Strain[0,1] = du_dy + dv_dx;
  Strain[0,2] = dw_dx + du_dz;
  Strain[1,0] = Strain[0,1];
  Strain[1,1] = 2*dv_dy - 2./3.*divu;
  Strain[1,2] = dv_dz + dw_dy;
  Strain[2,0] = Strain[0,2];
  Strain[2,1] = Strain[1,2];
  Strain[2,2] = 2*dw_dz - 2./3.*divu;
  """
  if verbosity > 1:
    print(f"Strain = {Strain}")
  """
  for i_disp in range(3): # set to zero
      for j_disp in range(3):
          WorkingMatrix[i_disp,j_disp] = 0.;
          StrainTrans[i_disp,j_disp] = 0.;  
          #if verbosity > 1:
          #  print("Livescue on ");
          Strain[i_disp,j_disp] = - Strain[i_disp,j_disp] ; 

  #if verbosity > 1:
  #  print(f"Check negative of strain rate calculated\n{Strain}");
  #  print("Test each of the B field conditions to try break the code");

  if (B < EffectiveZero) or (isoOveride): # Do we have a special case of B=0?
    #if isoOveride:  print("IsoOveride enforced visc stress tensor")
    ViscStress[0,0] = -eta0*Strain[0,0];
    
    ViscStress[0,1] = -eta0*Strain[0,1];
    ViscStress[1,0] = ViscStress[0,1];

    ViscStress[0,2] = -eta0*Strain[0,2];
    ViscStress[2,0] = ViscStress[0,2];

    ViscStress[1,1] = -eta0*Strain[1,1];
    ViscStress[1,2] = -eta0*Strain[1,2];
    ViscStress[2,1] = ViscStress[1,2];

    ViscStress[2,2] = -eta0*Strain[2,2];

  elif ( (abs(Bx) < EffectiveZero) and (abs(By) < EffectiveZero) and (abs(Bz) > EffectiveZero) ) :# case of all z aligned B field - no tranformation
    ViscStress[0,0]=-1/2*eta0*(Strain[0,0] + Strain[1,1]) \
                    - 1/2*eta1*(Strain[0,0] - Strain[1,1])\
                    -eta3*(Strain[0,1]);
    
    ViscStress[0,1]=-eta1*Strain[0,1] + 1/2*eta3*(Strain[0,0] - Strain[1,1]);

    ViscStress[1,0]= ViscStress[0,1];

    ViscStress[0,2]=-eta2*Strain[0,2] - eta4*Strain[1,2];

    ViscStress[2,0]= ViscStress[0,2];

    ViscStress[1,1]=-1/2*eta0*(Strain[0,0] + Strain[1,1])\
                    -1/2*eta1*(Strain[1,1] - Strain[0,0]) + eta3*Strain[0,1];

    ViscStress[1,2]=-eta2*Strain[1,2] + eta4*Strain[0,2];

    ViscStress[2,1]=ViscStress[1,2];

    ViscStress[2,2]=-eta0*Strain[2,2];

  else:# generic case of non zero and non cartesian z orientated field 

    # Multiplying Q' (Transpose) by W (StressStrain)
    """
    for i_disp in range(3):
        for j_disp in range(3):
            for k_disp in range(3):
                WorkingMatrix[i_disp,j_disp] += TransT[i_disp,k_disp]*Strain[k_disp,j_disp];
    #TODO replace with just matmul when checks complete:            
    if (np.matmul( TransT, Strain) != WorkingMatrix).any(): sys.exit("matrix multiplcation error:TransT*Strain")
    """
    WorkingMatrix = np.matmul( TransT, Strain)

    # Multiplying Q'W by Q
    """
    for i_disp in range(3): 
        for j_disp in range(3):
            for k_disp in range(3):
              StrainTrans[i_disp,j_disp] += WorkingMatrix[i_disp,k_disp] * Trans[k_disp,j_disp];
   
    if (np.matmul(WorkingMatrix, Trans) != StrainTrans).any(): sys.exit("matrix multiplcation error: QtWQ")
    """
    StrainTrans = np.matmul(WorkingMatrix, Trans)
  
    #Populate visc stress tensor in cartesian normal frame
    ViscStressTrans[0,0]=-1/2*eta0*\
            (StrainTrans[0,0] + StrainTrans[1,1])\
            -1/2*eta1*\
            (StrainTrans[0,0] - StrainTrans[1,1])\
            -eta3*(StrainTrans[0,1]);
    
    ViscStressTrans[0,1]=-eta1*StrainTrans[0,1]\
            +1/2*eta3*\
            (StrainTrans[0,0] - StrainTrans[1,1]);
  
    ViscStressTrans[1,0]= ViscStressTrans[0,1];
  
    ViscStressTrans[0,2]=-eta2*StrainTrans[0,2]\
            - eta4*StrainTrans[1,2];
  
    ViscStressTrans[2,0]= ViscStressTrans[2,0];
  
    ViscStressTrans[1,1]=-1/2*eta0*\
            (StrainTrans[0,0] + StrainTrans[1,1])\
            -1/2*eta1*\
            (StrainTrans[1,1] - StrainTrans[0,0])\
            +eta3*StrainTrans[0,1];
  
    ViscStressTrans[1,2]=-eta2*StrainTrans[1,2] +\
            eta4*StrainTrans[0,2];
  
    ViscStressTrans[2,1]=ViscStressTrans[1,2];
  
    ViscStressTrans[2,2]=-eta0*StrainTrans[2,2];
  
    for  i_disp in range(3):
        for  j_disp in range(3):
            WorkingMatrix[i_disp,j_disp] = 0.;
            ViscStress[i_disp,j_disp] = 0.;
    
    # Multiplying Q  Trans  by PI' (ViscStressTrans)
    """
    for i_disp in range(3):
      for j_disp in range(3):
        for k_disp in range(3):
          WorkingMatrix[i_disp,j_disp] += Trans[i_disp,k_disp]*ViscStressTrans[k_disp,j_disp];
     
    if (np.matmul(Trans, ViscStressTrans) != WorkingMatrix).any(): pdb.set_trace(); sys.exit("matrix multiplcationis fucked")
    """
    WorkingMatrix = np.matmul(Trans, ViscStressTrans)

    # Multiplying Q*PI' by Q^T
    """
    for i_disp in range(3):
      for j_disp in range(3):
        for k_disp in range(3):
          ViscStress[i_disp,j_disp] += WorkingMatrix[i_disp,k_disp] * TransT[k_disp,j_disp]; 
    if (np.matmul(WorkingMatrix, TransT) != ViscStress).any(): pdb.set_trace();sys.exit("matrix multiplcationis fucked")
    """
    ViscStress = np.matmul(WorkingMatrix, TransT)

    #Storing
    #NOTE STORAGE ACCORDING TO THE TAUXX, TAUYY, TAUZZ, TAUXY OR TAUYX,
    # TAUYZ OR TAUZY, TAUXZ OR TAUZX
  ViscTens[0] = ViscStress[0,0];
  ViscTens[1] = ViscStress[1,1];
  ViscTens[2] = ViscStress[2,2];
  ViscTens[3] = ViscStress[0,1];
  ViscTens[4] = ViscStress[1,2];
  ViscTens[5] = ViscStress[0,2];

  q_flux = np.zeros((3));

  if isoOveride: # do the normal version
    #print("isoOveride active q")
    if (name == eleName):
      q_flux[0] = beta1*duVec[0] -kappa1*dT_dVec[0]
      q_flux[1] = beta1*duVec[1] -kappa1*dT_dVec[1]
      q_flux[2] = beta1*duVec[2] -kappa1*dT_dVec[2]

    elif(name == ionName):
      q_flux[0] = -kappa1*dT_dVec[0]
      q_flux[1] = -kappa1*dT_dVec[1]
      q_flux[2] = -kappa1*dT_dVec[2]
  else: # do the anisotropic 
    if (name == eleName):
      #qu_e_temp =  ;
      #qt_e_temp = ;
      #q_flux[0] = qt_e_temp + qu_e_temp ;
      q_flux[0] = beta1*u_para[0] + beta2*u_perp[0] + beta3*u_chev[0] + \
                  -kappa1*TG_para[0] - kappa2*TG_perp[0] - kappa3*TG_chev[0]
      #TODO fix up the backward differences for the higher dimensions.
  
      q_flux[1] = beta1*u_para[1] + beta2*u_perp[1] + beta3*u_chev[1] \
                 -kappa1*TG_para[1] -kappa2*TG_perp[1] - kappa3*TG_chev[1];
  
      q_flux[2] = beta1*u_para[2] + beta2*u_perp[2] + beta3*u_chev[2] \
                 -kappa1*TG_para[2] -kappa2*TG_perp[2] - kappa3*TG_chev[2];
    elif(name == ionName):
      q_flux[0] = -kappa1*TG_para[0] - kappa2*TG_perp[0] + kappa3*TG_chev[0];
  
      q_flux[1] = -kappa1*TG_para[1] - kappa2*TG_perp[1] + kappa3*TG_chev[1];
  
      q_flux[2] = -kappa1*TG_para[2] - kappa2*TG_perp[2]+ kappa3*TG_chev[2];
  
  if (verbosity > 4):
    print("qVector\t" , q_flux[0], "\nqVector\t", q_flux[1], "\nqVector\t", q_flux[2]);
    print( "ViscTens\t" , ViscTens[0] ,"\nviscTens\t" , ViscTens[1], "\nviscTens\t" , 
                          ViscTens[2] , "\nviscTens\t" ,ViscTens[3], "\nviscTens\t" , 
                          ViscTens[4] , "\nviscTens\t" , ViscTens[5]);

  return (ViscTens, q_flux)

def braginskiiSourceWrapper(din):
  xcells = din['xcells']; ycells = din['ycells']
  ionName = din["ionName"]; eleName = din["eleName"]
  srcDst = {}

  #TODO trying to vectorise the function
  srcDst[ionName], srcDst[eleName] = \
    braginskiiSource(din["Q"], din["QD"], din["ionName"], din["eleName"], din["fieldName"], 
                     din["Debye"], din["Larmor"], din["lightspeed"], din["xref"], 
                     din["n0ref"], din["mref"], din["rhoref"], din["Tref"], din["uref"], 
                     din["dx"], din["isoOveride"], din["verbosity"])
  return srcDst

def braginskiiSource(Q, QD, ionName, electronName, fieldName, Debye, Larmor, lightspeed, 
                     x_ref, n0_ref, m_ref, rho_ref, T_ref, u_ref, dx, 
                     isoOveride=False, verbosity = 1):

  Density=0; Xvel = 1; Yvel=2; Zvel=3; Prs=4; Temp=5; Alpha=6;
  xD = 0; yD = 1; zD = 2; xB = 3; yB = 4; zB = 5; muIdx = 6; epIdx = 7;

  EffectiveZero = 1e-14;

  ep = Q[fieldName][epIdx]
  Ex = Q[fieldName][xD]/ep
  Ey = Q[fieldName][yD]/ep
  Ez = Q[fieldName][zD]/ep
  Bx = Q[fieldName][xB]
  By = Q[fieldName][yB]
  Bz = Q[fieldName][zB]
  B_xyz = np.array([Bx, By, Bz]);

  B = Bx*Bx + By*By + Bz*Bz;

  """
  if isoOveride:
    Bx[:] = 0.; By[:] = 0. ; Bz[:] = 0. ;
    B_xyz[:] = 0;
    B[:] = 0. 
    pdb.set_trace()
  """

  # NOTE that below here everything was designed to be a cell based operation, have to check the operations that van be vectorised and those which cannot

  #print(f"dimensions of B:\t{B.shape}")

  #TODO Vectorise the magnetic field sorting
  if ((B < 0.0).any()):
    print("MFP_braginskii.cpp ln 203 - Negative magnetic field error");
    sys.ext("ln 204 MFP_braginskii.cpp");
  else:
    B_pp = np.where(B < EffectiveZero, 0., 1./np.sqrt(B) )
    B_p  = np.where(B < EffectiveZero, 0., np.where(\
      np.logical_and( 
        np.logical_and((np.absolute(Bx)<EffectiveZero), (np.absolute(By)<EffectiveZero)), 
        (np.absolute(Bz)>EffectiveZero)),
      0., 1/np.sqrt(By*By + Bx*Bx)))


    """
    if (B < EffectiveZero):
      #print("\nzero magnetic field\n");
      B_pp = 0.;
      B_p  = 0.;
    elif ((xB < EffectiveZero) and (yB < EffectiveZero) and (zB > EffectiveZero)):
      print("\nzero x and y magnetic field\n");
      B_pp = 1/sqrt(B); # B prime prime 
      B_p  = 0.;
    else :  
      B_pp = 1/np.sqrt(B); # B prime prime 
      B_p  = 1/np.sqrt(Bx*Bx + By*By); # B prime 
    """
  bx_pp = Bx*B_pp; bx_p = Bx*B_p;
  by_pp = By*B_pp; by_p = By*B_p;
  bz_pp = Bz*B_pp; 

  #TODO vectorise
  B_unit = np.array([bx_pp, by_pp, bz_pp]); 

  # need to run this function to get the coefficients. 
  """These are the viscosity and conductivity coefficients
  T_i, eta0, eta1, eta2, eta3, eta4, kappa1, kappa2, kappa3 = \
  get_ion_coeffs(stateProps[ionName], stateProps[fieldName], stateProps[eleName],
                 stateU[offsetMap[ionName]:offsetMap[ionName]+stateProps[ionName].nVars], 
                 stateU[offsetMap[eleName]:offsetMap[eleName]+stateProps[eleName].nVars], B_xyz)


  T_e, eta0, eta1, eta2, eta3, eta4, kappa1, kappa2, kappa3, beta1, beta2, beta3  =\
  get_ion_coeffs(stateProps[eleName], stateProps[fieldName], stateProps[ionName],
                 stateU[offsetMap[eleName]:offsetMap[eleName]+stateProps[eleName].nVars], 
                 stateU[offsetMap[ionName]:offsetMap[ionName]+stateProps[ionName].nVars], B_xyz)
  """
  #TODO this will have to be a differently shaped structure now 
  u_para = {}  # Velocity parallel to B_unit
  u_perp = {}  # Velocity perpendicular to B_unit
  u_chev = {}  # unit vector perp to u and B_unit
  TG_para = {} #Temp grad parallel to B_unit 
  TG_perp = {} #Temp grad perpendicular to B_unit
  TG_chev = {} #unit vector perp to gradT and B_unit

  #Real B_p=0.,B_pp=0.,bx_pp=0.,by_pp=0.,bz_pp=0.,bx_p=0.,by_p=0., xB, yB, zB;
  # electron primitives 
  primEle = Q[electronName]
  rho_a =   primEle[Density];
  u_a =     primEle[Xvel];
  v_a =     primEle[Yvel];
  w_a =     primEle[Zvel];
  p_a =     primEle[Prs];
  T_a =     primEle[Temp];
  alpha_a = primEle[Alpha];

  m_a = primEle['mass']
  q_a = primEle['charge']
  n_a = rho_a/m_a;
  q_a2 = q_a*q_a;
  gam_a= primEle['gamma']

  
  #pull out the gradient values from the slopes stuff
  #dT_dVec = np.zeros((3)); #TODO change structure 
  
  """ Not used in source terms 
  drho_dx = primGrad[eleName][Density];
  du_dx   = primGrad[eleName][Xvel];
  dv_dx   = primGrad[eleName][Yvel];
  dw_dx   = primGrad[eleName][Zvel];
  dp_dx   = primGrad[eleName][Prs];
  """
  try:  
    dT_dVec = {}
    dT_dVec[0] = QD['Tedx'] # is this the handle used? dont think so 
    dT_dVec[1] = QD['Tedy']
    dT_dVec[2] = QD['Tedy']*0.
    #print("temp gradients used in sources")
  except:
    dT_dVec = {}
    dT_dVec[0] = np.gradient(primEle[Temp], dx, axis = 1);
    dT_dVec[1] = np.gradient(primEle[Temp], dx, axis = 0);
    dT_dVec[2] = dT_dVec[1]*0. 

  """
  drho_dy=0; du_dy=0; dv_dy=0; dw_dy=0; dp_dy=0; dT_dy=0; # the three dimensional stuff is needed bc i havent derived the algebraic expression for each property *(these are needed for the coordinate transformation)
  
  drho_dz=0; du_dz=0; dv_dz=0; dw_dz=0; dp_dz=0; dT_dz=0;
  """

  primIon = Q[ionName]
  rho_b =   primIon[Density];
  u_b =     primIon[Xvel];
  v_b =     primIon[Yvel];
  w_b =     primIon[Zvel];
  p_b =     primIon[Prs];
  T_b =     primIon[Temp];
  alpha_b = primIon[Alpha];

  m_b = primIon['mass']
  q_b = primIon['charge']
  n_b = rho_b/m_b;
  q_b2 = q_b*q_b;
  gam_b= primIon['gamma']

  duVec = {0:u_a - u_b, 1:v_a - v_b, 2:w_a - w_b}; #dU2 = du*du + dv*dv + dw*dw;

  #Braginskii directionality stuff. Real dot_B_unit_TG, dot_B_unit_U ;#temp variables
  dot_B_unit_U = bx_pp*duVec[0]+ by_pp*duVec[1]+ bz_pp*duVec[2];
  dot_B_unit_TG = bx_pp*dT_dVec[0] + by_pp*dT_dVec[1] + bz_pp*dT_dVec[2];

  #TODO vectorise
  for i_disp in range(3):#i_disp - i disposable 
    u_para[i_disp] = B_unit[i_disp]*dot_B_unit_U ;
    TG_para[i_disp]= B_unit[i_disp]*dot_B_unit_TG ;
    u_perp[i_disp] = duVec[i_disp] - u_para[i_disp];
    TG_perp[i_disp]= dT_dVec[i_disp] - TG_para[i_disp];  
    if (i_disp==0):
      u_chev[i_disp] = B_unit[1]*duVec[2]-B_unit[2]*duVec[1];
      TG_chev[i_disp]= B_unit[1]*dT_dVec[2]-B_unit[2]*dT_dVec[1];
    elif (i_disp==1):
      u_chev[i_disp] = -(B_unit[0]*duVec[2]-B_unit[2]*duVec[0]);
      TG_chev[i_disp]= -(B_unit[0]*dT_dVec[2]-B_unit[2]*dT_dVec[0]);
    else:
      u_chev[i_disp] = B_unit[0]*duVec[1]-B_unit[1]*duVec[0];
      TG_chev[i_disp]= B_unit[0]*dT_dVec[1]-B_unit[1]*dT_dVec[0];
  if (verbosity > 4):
    print("\nBunit\t", B_unit[0], "\nBunit\t", B_unit[1], "\nBunit\t", B_unit[2])
    print("u_rel\t", duVec[0], "\nu_rel\t", duVec[1], "\nu_rel\t", duVec[2])
    print("dot_B_unit_U\t", dot_B_unit_U)
    print("dot_B_unit_TG\t", dot_B_unit_TG)
    print("u_para[0]\t", u_para[0], "\nu_para\t", u_para[1], "\nu_para\t", u_para[2])
    print("u_perp[0]\t", u_perp[0], "\nu_perp\t", u_perp[1], "\nu_perp\t", u_perp[2])
    print("u_chev[0]\t", u_chev[0], "\nu_chev[1]\t", u_chev[1], "\nu_chev[2]\t", u_chev[2])
    print("TG_para[0]\t", TG_para[0], "\nTG_para[1]\t", TG_para[1], "\nTG_para[2]\t", TG_para[2])
    print("TG_perp[0]\t", TG_perp[0], "\nTG_perp[1]\t", TG_perp[1], "\nTG_perp[2]\t", TG_perp[2])
    print("TG_chev[0]\t", TG_chev[0], "\nTG_chev[1]\t", TG_chev[1],"\nTG_chev[2]\t", TG_chev[2])

  #---------------Braginskii Momentum source 
  #Real alpha_0, alpha_1, alpha_2, beta_0, beta_1, beta_2, t_c_a;
  p_lambda = get_coulomb_logarithm();

  #note a referes to electrons and b refers to ions 
  alpha_0, alpha_1, alpha_2, beta_0, beta_1, beta_2, t_c_a = \
    get_alpha_beta_coefficients(m_a, T_a, q_a, q_b, n_a, p_lambda, Debye, Larmor, x_ref, 
    n0_ref, m_ref, rho_ref, T_ref, u_ref, Bx, By, Bz, verbosity);

  if (verbosity>4):
    print("alpha_0\t", alpha_0, "\nalpha_1\t", alpha_1, "\nalpha_2\t", alpha_2   
             , "\nbeta_0\t", beta_0, "\nbeta_1\t", beta_1, "\nbeta_2\t", beta_2 
             , "\nt_c_a\t", t_c_a)

  #print(f"primIon\n{primIon}\nprimEle\t{primEle}")
  #print("Testing Daryl's collision frequency")

  #c1 = 0.02116454531141366/(n0*(Debye**4)); # sqrt(2)/(12*pi^(3/2))*...
  #c2 = 0.4135669939329334; # (2/(9*pi))**(1/3.0)
  #c3 = 2.127692162140974*n0; # (4/3)*n0*sqrt(8/pi)

  #m_ei = m_a*m_b/(m_a+m_b);
  #coeff_1 = c1*((q_a**2*q_b**2*p_lambda)/(m_ei*m_b));
  #nu_ie_rambo = n_a*coeff_1*(c2* np.sum(duVec**2) + T_a/m_a + T_b/m_b)**(-1.5);
  #coeff_1 = c1*((q_a**2*q_b**2*p_lambda)/(m_ei*m_a));
  #nu_ei_rambo = n_a*coeff_1*(c2* np.sum(duVec**2) + T_a/m_a + T_b/m_b)**(-1.5);
  #print(f"nu_brag:\t{1/t_c_a}\tnu_ie_rambo\t{nu_ie_rambo}\tnu_ei_rambo\t{nu_ei_rambo}") 

  #print("Using ion-electron collision rate")
  #t_c_a = 1/nu_ie_rambo;
  R_u = {} #np.zeros((3)); 
  R_T = {} #np.zeros((3));

  if isoOveride:
    #print("isoOveride active source terms")
    R_u[0] = -alpha_0*duVec[0]
    R_u[1] = -alpha_0*duVec[1]
    R_u[2] = -alpha_0*duVec[2]
    #Thermal force
    R_T[0] = -beta_0*dT_dVec[0]
    R_T[1] = -beta_0*dT_dVec[1]
    R_T[2] = -beta_0*dT_dVec[2]
   
  else:
    #frictional force
    R_u[0] = -alpha_0*u_para[0] - alpha_1*u_perp[0] + alpha_2*u_chev[0];
    R_u[1] = -alpha_0*u_para[1] - alpha_1*u_perp[1] + alpha_2*u_chev[1];
    R_u[2] = -alpha_0*u_para[2] - alpha_1*u_perp[2] + alpha_2*u_chev[2];
    #Thermal force
    R_T[0] = -beta_0*TG_para[0] - beta_1*TG_perp[0] - beta_2*TG_chev[0];
    R_T[1] = -beta_0*TG_para[1] - beta_1*TG_perp[1] - beta_2*TG_chev[1];
    R_T[2] = -beta_0*TG_para[2] - beta_1*TG_perp[2] - beta_2*TG_chev[2];

  #Thermal equilibration
  Q_delta = 3*m_a/m_b*n_a/t_c_a*(T_a-T_b);
  #Q_fric  = (R_u[0]+R_T[0])*duVec[0] + (R_u[1]+R_T[1])*duVec[1] + (R_u[2]+R_T[2])*duVec[2] 
  Q_fric  = (R_u[0]+R_T[0])*u_b + (R_u[1]+R_T[1])*v_b + (R_u[2]+R_T[2])*w_b

  Q_i = Q_delta  + Q_fric
  Q_e = -Q_delta - Q_fric ;
  
  #print("Q_e added to ion also")
  #print("Check that the parrallel component is applied for each direction in the event of a zero magnetic field"); #pdb.set_trace()
 

  if verbosity > 4:
    print(f"\nQ_fric:\t{Q_fric}\nQ_delta:\t{Q_delta}")
    print(f"R_u:\t{R_u[0]}\nR_u:\t{R_u[1]}\nR_u:\t{R_u[2]}")
    print(f"R_T:\t{R_T[0]}\nR_T:\t{R_T[1]}\nR_T:\t{R_T[2]}")
    print(f"dT_dVec\t{dT_dVec[0]}\ndT_dVec\t{dT_dVec[1]}\ndT_dVec\t{dT_dVec[2]}")
  srcDst_ele = np.zeros((Ex.shape[0], Ex.shape[1], 5));
  srcDst_ion = np.zeros((Ex.shape[0], Ex.shape[1], 5));
  Xmom, Ymom, Zmom, EdenQfric, EdenQdelta = 0, 1, 2, 3, 4

  srcDst_ele[:, :, Xmom] = R_u[0]+R_T[0];
  srcDst_ele[:, :, Ymom] = R_u[1]+R_T[1];
  srcDst_ele[:, :, Zmom] = R_u[2]+R_T[2];
  srcDst_ele[:, :, EdenQfric]  = -Q_fric ;
  srcDst_ele[:, :, EdenQdelta] = -Q_delta ;
  #note here the b is the ion
  srcDst_ion[:, :, Xmom] = -R_u[0] - R_T[0];
  srcDst_ion[:, :, Ymom] = -R_u[1] - R_T[1];
  srcDst_ion[:, :, Zmom] = -R_u[2] - R_T[2];
  srcDst_ion[:, :, EdenQfric]  = Q_fric
  srcDst_ion[:, :, EdenQdelta] = Q_delta 

  for i in range(4):
    if ( (np.isnan(srcDst_ion[i])).any() ):
      pdb.set_trace()
      sys.exit("Source term nan error")#pdb.set_trace()
    if ( (np.isnan(srcDst_ele[i])).any() ):
      pdb.set_trace()
      sys.exit("Source term nan error")#pdb.set_trace()
  #print(srcDst_ion[:,:,0] + srcDst_ele[:,:,0])

  return srcDst_ion, srcDst_ele

def getIonCoeffs(ionName, eleName, emName, Q, xi, yi, B_xyz, Debye, Larmor, n0_ref, 
                 x_ref, u_ref, verbosity=1):
  """  
  This function will be isntantiated within a sub process which will have access to
  parent process name space that will include the variable index variables e.g. 
  Density etc
  """
  Density=0; Xvel = 1; Yvel=2; Zvel=3; Prs=4; Temp=5; Alpha=6;
  x_D = 0; y_D = 1; z_D = 2; x_B = 3; y_B = 4; z_B = 5; muIdx = 6; epIdx = 7;

  if verbosity > 9:  print("Debug and test getIonCoeffs")
  #Extract and assign parameters from Q_i and Q_e
  #--- electron state and properties required for calcs -------Note move this
  primEle = Q[eleName]; 
  alpha_e = primEle[Alpha][xi, yi];
  charge_e= primEle['charge'][xi, yi]; # electron propertis
  mass_e  = primEle['mass'][xi, yi];
  T_e     = primEle[Temp][xi, yi];
  nd_e    = primEle[Density][xi, yi]/mass_e;
  #--- ionand properties required for calcs
  primIon = Q[ionName]; 
  alpha_i = primIon[Alpha][xi, yi];
  charge_i= primIon['charge'][xi, yi]; # electron propertis
  mass_i  = primIon['mass'][xi, yi];
  T_i     = primIon[Temp][xi, yi];
  nd_i    = primIon[Density][xi, yi]/mass_i;
  #Magnetic field
  Bx = B_xyz[0]; By = B_xyz[1]; Bz = B_xyz[2];

  # See page 215 (document numbering) of Braginskii's original transport paper
  #Real t_collision_ion, lambda_i, p_lambda, omega_ci, omega_p;
  p_lambda = get_coulomb_logarithm();
  # absence of the boltzmann constant due to usage of nondimensional variables.
  # Note that here charge_i = e^4*Z^4
  t_ref = x_ref/u_ref;
  t_collision_ion = (Debye)**(4)*n0_ref\
                    *(12*np.sqrt(mass_i)*(3.14159*T_i)**(3./2.))\
                    /(p_lambda * (charge_i)**(4) * nd_i);

  omega_ci = charge_i * np.sqrt(B_xyz[0]**2 + B_xyz[1]**2 + B_xyz[2]**2)/mass_i/Larmor;
  omega_p  = np.sqrt(nd_i*charge_i*charge_i/mass_i/Debye/Debye) ;
  
  if ((1/t_collision_ion < omega_ci/10/2/3.14159) and 
    (1/t_collision_ion < omega_p/10/2/3.14159)): 
    if False: #(verbosity > 1):
      print("\ngetIonCoeffs ln 572 --- Ion collision frequency limited to" + 
          " minimum of plasma and cyclotron frequency");
      print(f"1/tau_i = {1/t_collision_ion}\tomega_ci = {omega_ci}\tomega_p = {omega_p}"); 
    t_collision_ion = 1/min(omega_ci/2/3.14159, omega_p/2/3.14159) ;

  x_coef = omega_ci*t_collision_ion;

  #TODO fix up coefficients here also with tabled depending atomic number

  delta_kappa = x_coef*x_coef*x_coef*x_coef + 2.700*x_coef*x_coef + 0.677;
  delta_eta   = x_coef*x_coef*x_coef*x_coef + 4.030*x_coef*x_coef + 2.330;
  delta_eta2  = 16*x_coef*x_coef*x_coef*x_coef + 4*4.030*x_coef*x_coef + 2.330;

  eta0 = 0.96*nd_i*T_i*t_collision_ion ;#* n0_ref;
  eta2 = nd_i*T_i*t_collision_ion*(6./5.*x_coef*x_coef+2.23)/delta_eta;
  eta1 = nd_i*T_i*t_collision_ion*(6./5.*(2*x_coef)*(2*x_coef)+2.23)/delta_eta2;
  eta4 = nd_i*T_i*t_collision_ion*x_coef*(x_coef*x_coef + 2.38)/delta_eta;
  eta3 = nd_i*T_i*t_collision_ion*(2*x_coef)*((2*x_coef)*(2*x_coef) + 2.38)/delta_eta2;
  #if (verbosity >= 4 ):
  #    print("\nIon 5 dynamic viscosity coefficients\n", eta0, eta1, eta2, eta3, eta4);

  #Any of the viscosities greater than the parallel viscosity?
  if ((eta0<0) or (eta1<0) or (eta2<0) or (eta3<0) or (eta4<0)): 
      if (eta0<=0):
        print("\neta0 = ", eta0 ,"\n");
      elif (eta1<=0):
        print( "\neta1 = ", eta1 ,"\n");
      elif (eta2<=0):
        print("\neta2 = ", eta2 ,"\n");
      elif (eta3<=0):
        print("\neta3 = ", eta3 ,"\n");
      elif (eta4<=0):
        print("\neta4 = ", eta4 ,"\n");
      
      sys.exit("mfp_viscous.cpp ln: 334 - Braginski Ion coefficients are non-physical");
  
  #From Braginskii OG paper page 250 of paper in journal heading 4
  # Kinetics of a simple plasma (Quantitative Analyis)
  # TODO Add in flexibility for different atomic numbers of the ion species used,
  # see Braginskii
  kappa1 = 3.906*nd_i*T_i*t_collision_ion/mass_i;
  kappa2 = (2.*x_coef*x_coef + 2.645)/delta_kappa*nd_i*T_i*t_collision_ion/mass_i;
  kappa3 = (5./2.*x_coef*x_coef + 4.65)*x_coef*nd_i*T_i*t_collision_ion/mass_i/delta_kappa;

  #if (verbosity >= 4 ):
  #  print("\nIon heat conductivity coefficients\n", kappa1, kappa2, kappa3);

  if ((kappa1<0) or (kappa2<0) or (kappa3 < 0)):
    sys.exit("mfp_viscous.cpp ln: 350 - Braginski Ion coefficients are non-physical");
  if (kappa1 < kappa2):
    kappa2 = kappa1;
    if (verbosity >= 4 ):
      print("\nmfp_viscous.cpp ln: 401 - ion kappa2 exceed kappp1\n");
  if (kappa1 < kappa3):
    kappa3 = kappa1;
    if (verbosity >= 4 ):
      print("\nmfp_viscous.cpp ln: 404 - ion kappa3 exceed kappp1\n");
  return (T_i, eta0, eta1, eta2, eta3, eta4, kappa1, kappa2, kappa3) ; 

def getElectronCoeffs(ionName, eleName, emName, Q, xi, yi, B_xyz, Debye, Larmor, \
                      n0_ref, x_ref, u_ref, verbosity=1): 
    Density=0; Xvel = 1; Yvel=2; Zvel=3; Prs=4; Temp=5; Alpha=6;
    x_D = 0; y_D = 1; z_D = 2; x_B = 3; y_B = 4; z_B = 5; muIdx = 6; epIdx = 7;

    #Real mass_i,mass_e,charge_i,charge_e,T_i,nd_i,nd_e,alpha_e,alpha_i;
    #State &istate = GD::get_state(idx);
    if verbosity > 9:  print("\nDebug and test getElectronCoeffs")
    #Extract and assign parameters from Q_i and Q_e
    #--- electron state and properties required for calcs -------Note move this
    #--- ionand properties required for calcs
    primIon = Q[ionName]; 
    alpha_i = primIon[Alpha][xi, yi];
    charge_i= primIon['charge'][xi, yi]; # electron propertis
    mass_i  = primIon['mass'][xi, yi];
    T_i     = primIon[Temp][xi, yi];
    nd_i    = primIon[Density][xi, yi]/mass_i;
    #--- ion state and properties required for calcs
    primEle = Q[eleName]; 
    alpha_e = primEle[Alpha][xi, yi];
    charge_e= primEle['charge'][xi, yi]; # electron propertis
    mass_e  = primEle['mass'][xi, yi];
    T_e     = primEle[Temp][xi, yi];
    nd_e    = primEle[Density][xi, yi]/mass_e;
    #Magnetic field
    Bx=B_xyz[0]; By=B_xyz[1]; Bz=B_xyz[2];

    # See page 215 (document numbering) of Braginskii's original transport paper
    #Real t_collision_ele, lambda_e, p_lambda, omega_ce, omega_p;
    p_lambda = get_coulomb_logarithm();

    t_collision_ele = (Debye)**(4)*n0_ref\
                      *(6*np.sqrt(2*mass_e)*(3.14159*T_e)**(3./2.)) /\
                      (p_lambda*(charge_e)**(4)*(charge_i/-charge_e)*nd_e); 

    #if (verbosity >= 4): 
    #  print("\nln 535 viscous t_c_e = ", t_collision_ele);

    omega_ce = -charge_e * np.sqrt(B_xyz[0]**2 + B_xyz[1]**2 + B_xyz[2]**2)/mass_e/Larmor;
    omega_p  = np.sqrt(nd_e*charge_e*charge_e/mass_e/Debye/Debye) ;

    #print("\nProps:", charge_e, B_xyz[0], B_xyz[1], B_xyz[2], mass_e, Larmor);

    if ((1/t_collision_ele < omega_ce/10/2/3.14159) and \
        (1/t_collision_ele < omega_p/10/2/3.14159)):
      if (verbosity > 9):
        #print("\ngetElectronCoeffs ln 676 --- Electron collision frequency limited " +
        #    "to minimum of plasma and cyclotron frequency");
        print(f"1/tau_e = {1/t_collision_ele}\tomega_ce = {omega_ce}\tomega_p = {omega_p}"); 

      t_collision_ele = 1/min(omega_ce/2/3.14159, omega_p/2/3.14159) ;

      if (verbosity > 9):
        print(f"New tau_e = {t_collision_ele}"); 

    x_coef = omega_ce*t_collision_ele;
    # TODO fix up these table 2 page 251 BT
    delta_0=3.7703; delta_1=14.79;
    delta_kappa= x_coef*x_coef*x_coef*x_coef+delta_1*x_coef*x_coef + delta_0;
    delta_eta  = x_coef*x_coef*x_coef*x_coef+13.8*x_coef*x_coef + 11.6;
    delta_eta2 = 16*x_coef*x_coef*x_coef*x_coef+4*13.8*x_coef*x_coef + 11.6;

    if (verbosity > 4): print(f"nd_e  = {nd_e}\tT_e = {T_e}\tt_collision_ele={t_collision_ele}")

    eta0 = 0.733*nd_e *T_e * t_collision_ele;
    eta2 = nd_e *T_e*t_collision_ele*(2.05*x_coef*x_coef+8.5)/delta_eta;
    eta1 = nd_e *T_e*t_collision_ele*(2.05*(2*x_coef)*(2*x_coef)+8.5)/delta_eta2;
    eta4 = -nd_e*T_e*t_collision_ele*x_coef*(x_coef*x_coef+7.91)/delta_eta;
    eta3 = -nd_e*T_e*t_collision_ele*(2*x_coef)*((2*x_coef)*(2*x_coef)+7.91)/delta_eta2;

    if (False) and (x_coef < 1e-8 ): #hall correction abandoned 
      eta2 = eta2/x_coef/x_coef;
      eta1 = eta1/x_coef/x_coef;
      eta4 = eta4/x_coef;
      eta3 = eta3/x_coef;

    #maybe benchmark to see if pow should be used idunnobruda
    #if (verbosity >= 4 ):
    #    print("\nElectron 5 dynamic viscosity coefficients\n", eta0, eta1, eta2, eta3, eta4);

    #From Braginskii OG paper page 250 of paper in journal heading 4
    # Kinetics of a simple plasma (Quantitative Analyis)
    #TODO change coefficient values for different Z values
    # Currently set for a hydrogen  plasma
    BT_gamma_0=11.92/3.7703;BT_gamma_0_p=11.92;BT_gamma_1_p=4.664;BT_gamma_1_pp=5./2.;
    BT_gamma_0_pp=21.67;

    kappa1 = nd_e*T_e*t_collision_ele/mass_e*BT_gamma_0;
    kappa2 = (BT_gamma_1_p*x_coef*x_coef+BT_gamma_0_p)/\
              delta_kappa*nd_e*T_e*t_collision_ele/mass_e;
    kappa3=(BT_gamma_1_pp*x_coef*x_coef+BT_gamma_0_pp)*x_coef*nd_e*T_e*t_collision_ele\
           /mass_e/delta_kappa;

    if ((kappa1<0.) or (kappa2<0.) or (kappa3 < 0.)):
      #amrex::Abort("mfp_viscous.cpp ln: 673 - Braginski Ion coefficients are non-physical");
      print("\nmfp_viscous.cpp ln: 718 - Braginski Ion coefficients are non-physical\n");
      print("\n", kappa1, "\n", kappa2, "\n", kappa3, "\nomega_ce = ", omega_ce);

    if (kappa1 < kappa2):
      if (verbosity >= 4):
        print("mfp_viscous.cpp ln: 688 - electron kappa2 exceed kappp1");
      kappa2 = kappa1;
    if (kappa1 < kappa3):
      if (verbosity >= 4):
        print("mfp_viscous.cpp ln: 694 - electron kappa3 exceed kappp1");
      kappa3 = kappa1;

    #--- beta terms for the thermal component of thermal heat flux of the electrons.
    b_0 = 0.711; b_0_pp = 3.053; b_0_p=2.681; b_1_p=5.101; b_1_pp=3./2.;
    beta1 = nd_e*b_0*T_e;
    beta2 = nd_e*(b_1_p*x_coef*x_coef+b_0_p)/delta_kappa*T_e;
    beta3 = nd_e*x_coef*(b_1_pp*x_coef*x_coef+b_0_pp)/delta_kappa*T_e;
    return (T_e, eta0, eta1, eta2, eta3, eta4, kappa1, kappa2, kappa3, beta1, beta2, beta3);

def get_alpha_beta_coefficients(mass_e, T_e, charge_e, charge_i, nd_e, p_lambda, Debye, Larmor, 
                                x_ref, n0_ref, m_ref, rho_ref, T_ref, u_ref, Bx, By, Bz, verbosity):
  # collision time nondimensional
  pi_num = 3.14159265358979323846;
  t_ref = x_ref/u_ref;
  t_c_e = (Debye)**(4)*n0_ref\
                    *(6*np.sqrt(2*mass_e)*(pi_num*T_e)**(3./2.)) /\
                    (p_lambda*(charge_e)**(4)*(charge_i/-charge_e)*nd_e); 

  omega_ce = -charge_e * np.sqrt( Bx*Bx + By*By + Bz*Bz ) / mass_e / Larmor; 
  omega_p  = np.sqrt(nd_e*charge_e*charge_e/mass_e/Debye/Debye) ;

  t_c_e = np.where( np.logical_and( np.less(1/t_c_e, omega_ce/10/2/pi_num), np.less(1/t_c_e, omega_p/10/2/pi_num)), 1/np.minimum(omega_ce/2/pi_num, omega_p/2/pi_num), t_c_e)
  """
  if ((1/t_c_e < omega_ce/10/2/pi_num) and \
      (1/t_c_e < omega_p/10/2/pi_num)):
    if False: #(verbosity > 1):
      print("\nget_alpha_beta ln 749 --- Electron collision frequency limited " +
          "to minimum of plasma and cyclotron frequency");
      print(f"1/tau_e = {1/t_c_e}\tomega_ce = {omega_ce}\tomega_p = {omega_p}"); 

    t_c_e = 1/min(omega_ce/2/pi_num, omega_p/2/pi_num) ;

    if False: #(verbosity > 1):
      print(f"New 1/tau_e = {1/t_c_e}"); 
  """
  x_coef = omega_ce*t_c_e;
  delta_0 =3.7703; delta_1 = 14.79;
  delta_coef = x_coef*x_coef*x_coef*x_coef + delta_1*x_coef*x_coef + delta_0;# TODO delta0 tables 

  a_0 = 0.5129; a_0_p =1.837; a_1_p =6.416; a_0_pp =0.7796; a_1_pp =1.704;

  alpha_0 = mass_e*nd_e/t_c_e*a_0;
  alpha_1 = mass_e*nd_e/t_c_e*( 1 - (a_1_p*x_coef*x_coef + a_0_p)/delta_coef);
  alpha_2 = mass_e*nd_e/t_c_e*x_coef*(a_1_pp*x_coef*x_coef+a_0_pp)/delta_coef;

  b_0 = 0.711; b_0_pp = 3.053; b_0_p=2.681; b_1_p=5.101; b_1_pp=3./2.;
  beta_0 = nd_e*b_0;
  beta_1 = nd_e*(b_1_p*x_coef*x_coef+b_0_p)/delta_coef;
  beta_2 = nd_e*x_coef*(b_1_pp*x_coef*x_coef+b_0_pp)/delta_coef;
  if  (verbosity >= 4):
    print("1/tau_e\t", 1/t_c_e, "\nomega_ce\t", omega_ce, "\nomega_p\t", omega_p, 
          "\nmass_e\t", mass_e, "\nnd_e\t", nd_e);

  return (alpha_0, alpha_1, alpha_2, beta_0, beta_1, beta_2, t_c_e)

def get_coulomb_logarithm():
    return 10.;


##################################-End Braginskii functions-###################################
def get_interface_mask(inputs):
    name = inputs['name']; tol = inputs['name']
    try:
      cmap = inputs['cmap']
    except:
      cmap =mpl.cm.gray_r
    h5file = inputs['h5file']
    level = inputs['level']
    max_res = inputs['max_res']
    window = inputs['window']

    rc = ReadBoxLib( h5file, level, window)
    cmap_gray_r2rgba = cmap
    try:
      x, tracer_array = rc.get_flat('alpha-%s'%name) 
    except:
      x, tracer_raw = rc.expression("tracer-%s"%name);
      x, rho = rc.expression("rho-%s"%name);
      tracer_array = tracer_raw/rho
     
    t_h = 0.5 + tol
    t_l = 0.5 - tol
    mask = (t_l <  tracer_array) & (tracer_array < t_h) 
    del tracer_raw, tracer_array, x_rho, y_rho, x_tracer, y_tracer, rho; gc.collect()
    return mask

def exactTimeFrame(inputs):
  t_interest = inputs['t_interest']
  fileNames = inputs['fileNames']
  level = 0 #inputs['level'] # only need lowest level to minimuse read time and get tiem 
  window = inputs['window']

  t_n = 0; tol = 0.
  for i in range(len(fileNames)):
    try:
      rc = ReadBoxLib(fileNames[i], 0, [[0., 0.01], [0, 0.01]])
    except:
      print("can't open file:\n", fileNames[i])
      pdb.set_trace()

    if i == 0:
      tol = abs(t_interest - rc.time)
    else:
      tol_new = abs(t_interest - rc.time)
      if tol > tol_new:
        t_n = i; tol = tol_new ;
      if tol < tol_new:
        #t_n = i
        break
    rc.close()
  return t_n, fileNames[t_n]

def find_time_frames(inputs):
    """
    1 - Every second of the contour data frames is used (useful for reducing the number
        contour plots, of those produced, displayed). 
    2 - set time inputs which we want to find the closest frame

    3 - Some default time ranges for the VEL and EM plots which have been used in the past  
    """
    frameType = inputs['frameType']
    fileNames = inputs['fileNames']
    ETF_inputs = {} 
    if frameType == 1: 
      dataIndex = inputs['dataIndex']
      n_time_slices = inputs['n_time_slices']  
      time_slices = inputs['time_slices']
      data_index = [ data_index[2*i] for i in range(int(len(data_index)/2)+1) ]
      n_time_slices = len(data_index)
      time_slices = range(n_time_slices) # [ time_slices[2*i] for i in range(int(len(data_index)/2)) ]
      data_index[0] = int(data_index[1]/2)
      return data_index, n_time_slices, time_slices 

    elif frameType == 2:
      times_wanted = inputs['times_wanted']
      ETF_inputs['fileNames'] = fileNames 
      ETF_inputs['level'] = inputs['level'] 
      ETF_inputs['window'] = inputs['window']
      data_index = []

      for t_interest in times_wanted:
        ETF_inputs['t_interest'] = t_interest 
        data_index.append(exactTimeFrame(ETF_inputs)[0])

      n_time_slices = len(data_index)
      time_slices = range(n_time_slices)      
      return data_index, n_time_slices, time_slices 

    elif frameType == 3: 
      ETF_inputs['fileNames'] = fileNames 
      ETF_inputs['level'] = level 
      ETF_inputs['max_res'] = max_res
      ETF_inputs['window'] = window

      data_index = []

      plot_input_dic = {0:(True, False, 0.005, 0.035), 1:( False, True, 0.005, 0.035), 
                      2:(True, False,0.04, 0.07 ), 3:(False, True, 0.04, 0.07 ), 
                      3:(False, True, 0.03, 0.05), 4:(False, True, 0.095, 0.121), 
                      5:(True, False, 0.05, 0.081), 6:(False, True, 0.2, 0.351), 
                      7:(False, True, 0.376, 0.526), 8:(True, False, 0.376, 0.526),
                      9:(True, False, 0.2, 0.351),# used in confirmation for ds=10 SRMI
                      10:(False, True, 0.05, 1.), 11:(True, False, 0.005, 0.15), 
                      12:(True, False, 0.0025, 1.005), 13:(True, False, 0.2, 0.235), 
                      14:(False, True, 0.6, 0.85),
                       }
      input_index = 12

      start_time = plot_input_dic[input_index][2]
      finish_time = plot_input_dic[input_index][3]
      steps = 6
      step = (finish_time - start_time)/steps
      print('step:', step)
      data_times = [ start_time + step*i for i in range(steps+1) ]  
    
      for t_interest in data_times:
        ETF_inputs['t_interest'] = t_interest 
        data_index.append(exactTimeFrame(ETF_inputs)[0])

      print("EM/VEL contour indexes:", data_index)
      n_time_slices = len(data_index)
      time_slices = range(n_time_slices)      
      return data_index, n_time_slices, time_slices 

def interfacePeak(prop, threshold, iStart, iEnd, jIndex):
  """jIndex == False indicates an average, see usage for tempalte"""
  propList = []
  propSwitch = False; 
  for i in range(iStart, iEnd):
    if not jIndex:
      if prop[i-iStart] > threshold and not propSwitch and len(propList) < 2: 
        propSwitch = True; propList.append(i)
      if prop[i-iStart] < threshold  and propSwitch and len(propList) < 2:
        propList.append(i)
    else:
      if prop[i,j] > threshold and not propSwitch and len(propList) < 2: 
        propSwitch = True; propList.append(i)
      if prop[i,j] < threshold and propSwitch and len(propList) < 2:
        propList.append(i)

  return propList

def interfacePeakSigned(prop, thresholdNP, iStart, iEnd, iN, iP, jIndex):  
  """
  iN - negative peak 
  iP - positive peak 
  jIndex == False indicates an average, see usage for tempalte
  """
  thresholdN = thresholdNP[0]
  thresholdP = thresholdNP[1]
  if thresholdN < 0 or thresholdP < 0: xxxx

  propList = []
  propSwitch = False; 
  for i in range(iStart, iEnd):
    if not jIndex: pVal = prop[i-iStart] # assign thev alue to test 
    else: pVal = prop[i,j]

    if iN < iP : # is the negative peak comes before the positive peak 
      if i < iN + iStart and pVal < 0 and abs(pVal) > thresholdN and not propSwitch:
        propSwitch = True; propList.append(i)
      if pVal < thresholdP and pVal> 0 and propSwitch and len(propList) < 2 and i > iP + iStart:
        propList.append(i)
        break 
    elif iP < iN:
      if i < iP +iStart and pVal > 0 and pVal > thresholdP and not propSwitch:
        propSwitch = True; propList.append(i)
      if pVal < 0 and abs(pVal) < thresholdN and propSwitch and len(propList) < 2 and i > iN + iStart:
        propList.append(i)
        break 
     
    else: xxx #Somethings wrong

  return propList




if __name__ == '__main__':
  print('Hello, \n begin test.')














  print('End test.')
