import numpy as np, pdb, sys
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

from matplotlib import animation
from matplotlib.animation import FuncAnimation

cmd_folder = "../../../vis" #/home/kyriakos/Documents/Code/000_cerberus_dev/githubRelease-cerberus/cerberus/vis"

analyticalHartmann_folder = "/media/H_drive/000_PhD/002_Hartmann_flow/"

if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
    
from get_boxlib import ReadBoxLib, get_files
import PHM_MFP_Solver_Post_functions_v6 as PHM_MFP

if analyticalHartmann_folder not in sys.path:
    sys.path.insert(0, analyticalHartmann_folder)
import hartmannFlow_v2 as HA

#Data extraction 
def nPlot(f):
  return int(f.split('plt')[1])
 
def getPlotData(fileInterval):
  # get a list of all the files in this directory
  files = get_files('.', include=['.plt'], exclude=["temp", "chk"], times=[], tol=1e-4, get_all=True)
  #sort files
  
  files = sorted(files, key=nPlot)
  #files = [files[len(files)-1]]
  timeFilter =False 
  filesUse = []
  print(f"\n=====\nFiltering every {intervalVal} file\n=====\n")

  for i in range(len(files)):
    if i%fileInterval != 0: 
      #print(f"Exclude file at step:\t {nPlot(files[i])}")
      continue
    filesUse.append(files[i])
  files = filesUse.copy() + [files[i]]

  plotData = {}; timeIndex = 0
  if timeFilter: print("\n=====\nFilter times\n=====\n")
  for f in files:
    print(f"Extracting:\t{f}")
    rh5 = ReadBoxLib(f, max_level=-1)
    t = rh5.time; 
    if timeFilter and t < 10.: 
      print(f"Exclude time time: {t}\t at step:\t {nPlot(f)}")
      continue
   
    plotData[timeIndex] = {}
    plotData[timeIndex]['t'] = t;

    plotData[timeIndex]['x'], rhoc_i, plotData[timeIndex]['ux_i'], plotData[timeIndex]['uy_i'], plotData[timeIndex]['uz_i'], Jx_i, Jy_i, Jz_i = PHM_MFP.get_charge_number_density(rh5, 'ions')
 
    x, rhoc_e, plotData[timeIndex]['ux_e'], plotData[timeIndex]['uy_e'], plotData[timeIndex]['uz_e'], Jx_e, Jy_e, Jz_e = PHM_MFP.get_charge_number_density(rh5, 'electrons')
    plotData[timeIndex]['rhoc'] = rhoc_i + rhoc_e
    plotData[timeIndex]['Jx'] = Jx_e + Jx_i;
    plotData[timeIndex]['Jy'] = Jy_e + Jy_i;
    plotData[timeIndex]['Jz'] = Jz_e + Jz_i;
    x, plotData[timeIndex]['Bx'] = rh5.get("x_B-field")
    x, plotData[timeIndex]['By'] = rh5.get("y_B-field")
    x, plotData[timeIndex]['Bz'] = rh5.get("z_B-field")
    x, plotData[timeIndex]['Dx'] = rh5.get("x_D-field")
    x, plotData[timeIndex]['Dy'] = rh5.get("y_D-field")
    x, plotData[timeIndex]['Dz'] = rh5.get("z_D-field")
    x, plotData[timeIndex]['dBxdx'] = rh5.get("x_B-field-dx")
    x, plotData[timeIndex]['dBydx'] = rh5.get("y_B-field-dx")
    x, plotData[timeIndex]['dBzdx'] = rh5.get("z_B-field-dx")
    x, plotData[timeIndex]['P_i'] = rh5.get("p-ions")
    x, plotData[timeIndex]['P_e'] = rh5.get("p-electrons")

    timeIndex += 1
    rh5.close()

  # get reference values from any file 
  rh5 = ReadBoxLib(files[0], max_level=1)
  x_ref = rh5.data["x_ref"]
  rho_ref= rh5.data["rho_ref"]
  m_ref= rh5.data["m_ref"]
  T_ref = rh5.data["T_ref"]
  skin_nd = rh5.data["skin_depth"]
  beta = rh5.data["beta"]
  rh5.close()
  print(x_ref, rho_ref, m_ref, T_ref, skin_nd, beta)

  #print(f"Final tiem By\n{plotData[timeIndex-1]['By']}")
  #print(f"Final tiem By\n{plotData[timeIndex-1]['Bz']}")
  return plotData, timeIndex, x_ref, rho_ref, m_ref, T_ref, skin_nd, beta # len(files)
  
def plotData(inputData, timeIndex):
  #Plotting 
  fluidTitles = [r"$J_x$", r"$J_y$", r"$J_z$", 
                 r"$u^e_x$", r"$u^e_y$", r"$u^e_z$", 
                 r"$u^i_x$", r"$u^i_y$", r"$u^i_z$"]; 
  fluidHandles = ['Jx', 'Jy', 'Jz', 'ux_e', 'uy_e', 'uz_e','ux_i', 'uy_i', 'uz_i']

  fieldTitles = [r"$B_x$", r"$B_y$", r"$B_z$", 
                 r"$\partial B_x/\partial x$", r"$\partial B_y/\partial x$", 
                 r"$\partial B_z \partial x$", r"$\rho_c$"]

  fieldHandles = ['Bx', 'By', 'Bz','dBxdx', 'dBydx', 'dBzdx', 'rhoc']
  padBoi = 40
  # fluid props
  fig_fluid = plt.figure(figsize=(60, 60)); ax_d = [] # derived solutions
  ax_fluid = []
  for (i, prop) in enumerate(fluidHandles):
      i_ax = i+1
      ax_fluid.append(fig_fluid.add_subplot(3,3, i_ax))
      ax_fluid[i_ax-1].plot(inputData[timeIndex]['x'], inputData[timeIndex][prop])
      ax_fluid[i_ax-1].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
      ax_fluid[i_ax-1].tick_params(axis='y', labelsize=10)
      ax_fluid[i_ax-1].set_ylabel(fluidTitles[i], fontsize=5)
      if i == 0: ax_fluid[i_ax-1].set_title(f"t = {inputData[timeIndex]['t']:.3f}")

    #if i < 2: ax_fluid[i_ax-1].set_xticks([])
    #else: ax_d[i_ax-1].set_xlabel('x')
  
  fig_fluid.tight_layout(pad=padBoi)
  #EM rops
  fig_field = plt.figure(figsize=(60, 60)); ax_d = [] # derived solutions
  ax_field = []
  for (i, prop) in enumerate(fieldHandles):
      i_ax = i+1
      ax_field.append(fig_field.add_subplot(3,3, i_ax))
      ax_field[i_ax-1].plot(inputData[timeIndex]['x'], inputData[timeIndex][prop])
      ax_field[i_ax-1].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
      ax_field[i_ax-1].tick_params(axis='y', labelsize=10)
      ax_field[i_ax-1].set_ylabel(fieldTitles[i], fontsize=5)
      if i == 0: ax_field[i_ax-1].set_title(f"t = {inputData[timeIndex]['t']:.3e}")
    #if i < 2: ax_field[i_ax-1].set_xticks([])
    #else: ax_d[i_ax-1].set_xlabel('x')
  fig_field.tight_layout(pad=padBoi)
  
 
  plt.show()

  return 

def plotAnimatedData(dataFrame): 
  global dataToPlot, ax_field, interval;
  fieldTitles = [r"$B_x$", r"$B_y$", r"$B_z$", r"$J_x$", r"$J_y$", r"$J_z$", 
                   r"$\partial B_x/\partial x$", r"$\partial B_y/\partial x$", 
                   r"$\partial B_z \partial x$", [r"$P_i$", r"$P_e$"], 
                  [r"$u^e_y$", r"$u^i_y$"], [r"$u^e_z$", r"$u^i_z$"]]
  fieldHandles = ['Bx', 'By', 'Bz', 'Jx', 'Jy', 'Jz','dBxdx', 'dBydx', 'dBzdx', ['P_i', 'P_e'], 
                  ['uy_e', 'uy_i'], ['uz_e', 'uz_i']]

  for (i, prop) in enumerate(fieldHandles):
      i_ax = i+1
      ax_field[i_ax-1].cla()
      loBound = minB[i] * (1-0.1*minB[i]/np.abs(minB[i]))
      hiBound = maxB[i] * (1+0.1*maxB[i]/np.abs(maxB[i]))
      if globalYlim:
        ax_field[i_ax-1].set_ylim(loBound, hiBound)
        #ax_field[i_ax-1].set_ylim(0.1*minB[i]/np.abs(minB[i]) + minB[i], 1.1*maxB[i])

      if type(prop) == type([]):
        maxP = 0; minP = 0; 
        if 'P_' in prop[0]: minP = 1;

        for (j, pr) in enumerate(prop):
          if '_e' in pr: lStyle = '-'
          else: lStyle = '--'
          ax_field[i_ax-1].plot(dataToPlot[dataFrame]['x'], 
                        dataToPlot[dataFrame][pr], label=fieldTitles[i][j], linestyle=lStyle)  
          ax_field[i_ax-1].legend()
          maxP = max( maxP, max(dataToPlot[dataFrame][pr]))
          minP = min( minP, min(dataToPlot[dataFrame][pr]))
      else:
        ax_field[i_ax-1].plot(dataToPlot[dataFrame]['x'], dataToPlot[dataFrame][prop])  
        ax_field[i_ax-1].set_ylabel(fieldTitles[i], fontsize=10)
        maxP = max(dataToPlot[dataFrame][prop])
        minP = min(dataToPlot[dataFrame][prop])
      ax_field[i_ax-1].set_yticks(np.linspace(minP, maxP, 5)) #, format='.2e')
      ax_field[i_ax-1].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))

  ax_field[1].set_title(f"t = {dataToPlot[dataFrame]['t']:.3e} / {dataToPlot[nFiles-1]['t']:.3e}")
  return 

def plotMillerComparisonAnimated(dataFrame): 
  global dataToPlot, ax1, interval, millerNormalise, plotAnalytical;
  #TODO make clear where this is and that it must be set
  v_plate = 5e-3; n0_nd = 1.; B0_nd = 1 ; T0_nd = 1

  #print(f"Hard coded v_plate:\t{v_plate}\tn0_nd:\t{n0_nd}")

  if not millerNormalise:
    Titles = [[r"$u^e_x$", r"$u^e_y$", r"$u^e_z$", r"$u^i_x$", r"$u^i_y$", r"$u^i_z$"],
              [r"$J_z$"], [r"$B_y$", r"$B_z$"]]
    millerNormaliseFactor = [1., 1., 1.]
  else:
    Titles = [[r"$u^e_x/V$", r"$u^e_y/V$", r"$u^e_z/V$",r"$u^i_x/V$",r"$u^i_y/V$",r"$u^i_z/V$"],
              [r"$J_z/nV$"], [r"$B_ye5$", r"$B_ze5$"]]
    millerNormaliseFactor = [v_plate, v_plate*n0_nd, 1e-5]

  if plotAnalytical:
    HA.defineSol(n0_nd, T0_nd, B0_nd, v_plate)

    xs_an, sFundamental_an, sx_an, Eigenvalues_an, Q_an, normaliseValues_an, skin_nd_an, \
    beta_an = HA.findFundamentalSol()

    ux_an = HA.findVelocitySol(sx_an, xs_an, Q_an)

    xs_an, xs_B_an, Jy_vals_an, Jz_vals_an, By_vals_an, Bz_vals_an, dBydx_vals_an, dBzdx_vals_an = \
      HA.findDerivedSol(xs_an, ux_an)
    props_an  = [[ux_an[0,:]/millerNormaliseFactor[0], ux_an[1,:]/millerNormaliseFactor[0], 
                  ux_an[2,:]/millerNormaliseFactor[0], ux_an[3,:]/millerNormaliseFactor[0]], 
                 [Jz_vals_an/millerNormaliseFactor[1]], 
                 [By_vals_an/millerNormaliseFactor[2], Bz_vals_an/millerNormaliseFactor[2]]]
    Titles_an = [[r"$An u^e_y/V$", r"$An u^e_z/V$", r"$An u^i_y/V$",r"$An u^i_z/V$"],
                 [r"$An J_z/nV$"], [r"$An B_ye5$", r"$An B_ze5$"]]
    lines_an = [['-', '-', '--', '--'], ['-'], ['-', '-']]
    colours_an = 'r' 
    x_an = [[xs_an]*4, [xs_an], [xs_B_an, xs_B_an]] 

  Handles = [['ux_e', 'uy_e', 'uz_e', 'ux_i', 'uy_i', 'uz_i'], ['Jz'], ['By', 'Bz']]
  for (i, props) in enumerate(Handles):
    i_ax = i + 1
    ax1[i_ax-1].cla()
    i_ax = i+1
    if globalYlim:
      ax1[i_ax-1].set_ylim(1.1*minVal[i], 1.1*maxVal[i])
    #ax1[i_ax-1].set_ylabel(Titles[i], fontsize=10)
    maxP = 0; minP = 0;
    for (j, prop) in enumerate(props):
        if '_e' in prop: lStyle = '-'
        else: lStyle = '--'
        ax1[i_ax-1].plot(dataToPlot[dataFrame]['x'], 
                        dataToPlot[dataFrame][prop]/millerNormaliseFactor[i], 
                        label=Titles[i][j], linestyle=lStyle)        
        maxP = max( maxP, max(dataToPlot[dataFrame][prop]/millerNormaliseFactor[i]))
        minP = min( minP, min(dataToPlot[dataFrame][prop]/millerNormaliseFactor[i]))
    if plotAnalytical: # plot analytical solution on the same axis 
      for (j_an, prop_an) in enumerate(props_an[i]):
        print(i, j_an)
        ax1[i_ax-1].plot(x_an[i][j_an], prop_an, label=Titles_an[i][j_an], 
          linestyle=lines_an[i][j_an], color=colours_an)

    ax1[i_ax-1].set_yticks(np.linspace(minP, maxP, 5)) #, format='.2e')
    ax1[i_ax-1].legend()
    ax1[i_ax-1].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    #del_x = max(dataToPlot[dataFrame]['x'])-min(dataToPlot[dataFrame]['x'])
    #del_y = 1.1*maxVal[i] - 1.1*minVal[i]
    #del_x_pow = np.log10(np.int(del_x))
    #del_y_pow = np.log10(np.int(del_y))
    #print((del_y/pow(10,del_y_pow) ) / (del_x/pow(10,del_x_pow)))
    #pdb.set_trace()
    #ax1[i_ax-1].set_aspect((del_y/pow(10,del_y_pow) ) / (del_x/pow(10,del_x_pow)) )
    ax1[1].set_title(f"t = {dataToPlot[dataFrame]['t']:.3e} / {dataToPlot[nFiles-1]['t']:.3e}")
  return 

def plotElectromagneticAnimated(dataFrame): 
  global dataToPlot, ax2, interval;
  
  Titles = [r"$B_x$", r"$B_y$", r"$B_z$", r"$D_x$", r"$D_y$", r"$D_z$", r"$\rho_c$"]
  Handles = ['Bx', 'By', 'Bz', 'Dx', 'Dy', 'Dz', 'rhoc']

  for (i, prop) in enumerate(Handles):
    i_ax = i + 1
    ax2[i_ax-1].cla()
    loBound = minVal[i] * (1-0.1*minVal[i]/np.abs(minVal[i]))
    hiBound = maxVal[i] * (1+0.1*maxVal[i]/np.abs(maxVal[i]))

    if globalYlim:
      ax2[i_ax-1].set_ylim(loBound, hiBound)
    #ax[i_ax-1].set_ylim(1.1*minVal[i], 1.1*maxVal[i])
    ax2[i_ax-1].set_ylabel(Titles[i], fontsize=10)
    ax2[i_ax-1].plot(dataToPlot[dataFrame]['x'], 
                        dataToPlot[dataFrame][prop], label=Titles[i])  
    ax2[i_ax-1].legend()
    ax2[i_ax-1].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
  
  ax2[1].set_title(f"t = {dataToPlot[dataFrame]['t']:.3e} / {dataToPlot[nFiles-1]['t']:.3e}")
  return 


if __name__ == "__main__":
  prefix = '20231030_Cerberus_comparisons_'
  globalYlim = False
  plotStill = True
  intervalVal = 1
  frameInterval = 100 #  miliseconds 
  repeatAnimation = True
  print(f"Retrieving data...")
  global dataToPlot, fig_ani, nFiles;
  dataToPlot, nFiles, x_ref, rho_ref, m_ref, T_ref, skin_nd, beta = getPlotData(intervalVal)
  print(f"\t...data retrieved. \n{nFiles} available.\n\nPlotting...")

  maxuye = []; maxuyi = []; maxuze = []; maxuzi = [];
  for timeIndex in range(nFiles):
    maxuye.append((np.gradient(dataToPlot[timeIndex]['uy_e'], dataToPlot[timeIndex]['x'])).max())
    maxuze.append((np.gradient(dataToPlot[timeIndex]['uz_e'], dataToPlot[timeIndex]['x'])).max())
    maxuyi.append((np.gradient(dataToPlot[timeIndex]['uy_i'], dataToPlot[timeIndex]['x'])).max())
    maxuzi.append((np.gradient(dataToPlot[timeIndex]['uz_i'], dataToPlot[timeIndex]['x'])).max())
  print(max(maxuye), "\n", max(maxuyi), "\n", max(maxuze), "\n", max(maxuzi))

  if True: # static figures 
    maxF = max(dataToPlot.keys())
    plotFiles = [0,1,maxF]
    #for i in range(nFiles):
    for i in plotFiles:
      print(f"next i:\t{i}\nn i's:\t{nFiles}")
      plotData(dataToPlot, i)

  global fig_ani1, fig_ani2, fig_ani2, ax_field, ax1, ax2;
  
  print("False here ln 311")
  if False: # Many em properties
    interval = intervalVal
    fieldTitles = [r"$B_x$", r"$B_y$", r"$B_z$", r"$J_x$", r"$J_y$", r"$J_z$", 
                   r"$\partial B_x/\partial x$", r"$\partial B_y/\partial x$", 
                   r"$\partial B_z \partial x$", [r"$P_i$", r"$P_e$"], 
                  [r"$u^e_y$", r"$u^i_y$"], [r"$u^e_z$", r"$u^i_z$"]]
    fieldHandles = ['Bx', 'By', 'Bz', 'Jx', 'Jy', 'Jz','dBxdx', 'dBydx', 'dBzdx', ['P_i', 'P_e'],
                  ['uy_e', 'uy_i'], ['uz_e', 'uz_i']]

    fig_ani1 = plt.figure(1, figsize=(20, 15)); 
    padBoi = 40
    fig_ani1.tight_layout(pad=padBoi)
    maxB = []; minB = []
    for props in fieldHandles:
      if type(props) != type([]): props = [props]
      maxVal = []; minVal = []

      for pr in props:
        maxVal.append( max([ (dataToPlot[i][pr]).max() for i in range(nFiles)]) )
        minVal.append( min([ (dataToPlot[i][pr]).min() for i in range(nFiles)]) )
      maxB.append(max(maxVal))
      minB.append(min(minVal))
      print(f"Property:\t{props}\tMax:\t{maxB[-1]}\tMin:\t{minB[-1]}")
    #print(f"Max val:\t{maxB}\nMin val:\t{minB}")
    ax_field = []
    for (i, prop) in enumerate(fieldHandles):
      i_ax = i+1
      ax_field.append(fig_ani1.add_subplot(4,3, i_ax))
      ax_field[i_ax-1].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
      ax_field[i_ax-1].tick_params(axis='y', labelsize=10)

    if plotStill:
      plotAnimatedData(max(dataToPlot.keys()))
      fig_ani1.savefig(prefix+'B_J_rho.png')
      plt.show()

    ani1 = FuncAnimation(fig=fig_ani1, func=plotAnimatedData, 
            frames=range(max(dataToPlot.keys())), interval=frameInterval, repeat=repeatAnimation)
    #plotAnimatedData(max(dataToPlot.keys()))
    #plt.show()

  if True: #Miller graph comparison 
    #===== Get analytical solution
    global plotAnalytical;
    plotAnalytical = False

    if plotStill == True: plotAnalytical = True

    global millerNormalise;
    millerNormalise = True; 
    interval = intervalVal
    fig_ani2 = plt.figure(2, figsize=(30, 15)); 
    padBoi = 40
    fig_ani2.tight_layout(pad=padBoi)
    if not millerNormalise:
      Titles = [[r"$u^e_x$", r"$u^e_y$", r"$u^e_z$", r"$u^i_x$", r"$u^i_y$", r"$u^i_z$"],
               [r"$J_z$"], [r"$B_y$", r"$B_z$"]]
    else:
      Titles = [[r"$u^e_x/V$", r"$u^e_y/V$", r"$u^e_z/V$", r"$u^i_x/V$", r"$u^i_y/V$", r"$u^i_z/V$"],
               [r"$J_z/nV$"], [r"$B_ye5$", r"$B_ze5$"]]
    Handles = [['ux_e', 'uy_e', 'uz_e', 'ux_i', 'uy_i', 'uz_i'], ['Jz'], ['By', 'Bz']]

    maxVal = []; minVal = []
    ax1 = []; i_ax = 1
    for props in [['uy_e', 'uz_e', 'uy_i', 'uz_i'], ['Jz'], ['By', 'Bz']]:
      val1 = []; val2 = []
      ax1.append(fig_ani2.add_subplot(3,1, i_ax))
      ax1[i_ax-1].tick_params(axis='y', labelsize=10)
      for prop in props:
        val1 += [(dataToPlot[i][prop]).max() for i in range(nFiles)]
        val2 += [(dataToPlot[i][prop]).min() for i in range(nFiles)]
      maxVal.append(max(val1))
      minVal.append(min(val2))
      #print(f"Properties:\t{props}\tMax:\t{maxVal[-1]}\tMin:\t{minVal[-1]}")
      i_ax += 1

    if plotStill:
      plotMillerComparisonAnimated(max(dataToPlot.keys()))
      fig_ani2.savefig(prefix+'_MillerGraphComparison.png')

    ani2 = FuncAnimation(fig=fig_ani2, func=plotMillerComparisonAnimated, 
            frames=range(max(dataToPlot.keys())), interval=frameInterval,repeat=repeatAnimation)
    #plotMillerComparisonAnimated(0)
    plt.show()
 
  print("False here ln 397")
  if False: 
    interval = intervalVal
    fig_ani3 = plt.figure(3, figsize=(30, 30)); 
    padBoi = 40
    fig_ani3.tight_layout(pad=padBoi)

    fieldTitles = [r"$B_x$", r"$B_y$", r"$B_z$", r"$D_x$", r"$D_y$", r"$D_z$", r"$\rho_c$"]
    fieldHandles = ['Bx', 'By', 'Bz', 'Dx', 'Dy', 'Dz', 'rhoc']

    maxVal = []; minVal = []
    ax2 = []; i_ax = 1
    for prop in fieldHandles:
      val1 = []; val2 = []
      if prop != 'rhoc': ax2.append(fig_ani3.add_subplot(3,3, i_ax))
      else: ax2.append(fig_ani3.add_subplot(3,1, 3))

      ax2[i_ax-1].tick_params(axis='y', labelsize=10)
      maxVal.append( max([(dataToPlot[i][prop]).max() for i in range(nFiles)]) ) 
      minVal.append( min([(dataToPlot[i][prop]).min() for i in range(nFiles)]) ) 
      #print(f"Property:\t{prop}\tMax:\t{maxVal[-1]}\tMin:\t{minVal[-1]}")
      i_ax += 1

    if plotStill:
      plotElectromagneticAnimated(max(dataToPlot.keys()))
      fig_ani3.savefig(prefix+'_EM_props.png')

    ani3 = FuncAnimation(fig=fig_ani3, func=plotElectromagneticAnimated, 
            frames=range(max(dataToPlot.keys())), interval=frameInterval,repeat=repeatAnimation)
    #plotElectromagneticAnimated(0)
    plt.show()

  #print("\n\nabout to plot\n\n")
  #fig_ani1.show()
  #fig_ani2.show()
  #fig_ani3.show()

  #pdb.set_trace()

