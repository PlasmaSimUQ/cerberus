
import sys
cmd_folder = "/home/kyriakos/Documents_ubuntu/000_refactor_cerberus/cerberus/vis"
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
    
from get_boxlib import ReadBoxLib, get_files


derived_functions_code_folder = "/home/kyriakos/Documents_ubuntu/000_refactor_cerberus/cerberus/vis"

if derived_functions_code_folder not in sys.path:
  sys.path.insert(0, derived_functions_code_folder)
import PHM_MFP_Solver_Post_functions_v6 as phmmfp

import numpy as np
import pylab as plt
import pdb   
# =============================================================================
# 
# =============================================================================

def exactTimeFrameCerberus(inputs):
  t_interest = inputs['t_interest']
  fileNames = inputs['fileNames']

  max_level=-1
  t_n = 0; tol = 0.
  for i in range(len(fileNames)):
    rc = ReadBoxLib(fileNames[i], max_level)
    if i == 0:
      tol = abs(t_interest - rc.time)
    else:
      tol_new = abs(t_interest - rc.time)
      if tol > tol_new:
        t_n = i; tol = tol_new ;
      #if tol < tol_new:
        #t_n = i
        #break
    rc.close()  
    #if t_n > len(fileNames):
    #pdb.set_trace()
  return t_n, fileNames[t_n]

# =============================================================================
# 
# =============================================================================

properties = ["rho-", "p-", "T-","x_vel-"]

directories = [
"/home/kyriakos/Documents_ubuntu/000_refactor_cerberus/cerberus/Exec/testing/Braginskii_Riemann_2F", 
    ]

namesList = [ ["ion", "electron"], ["ion", "electron"]]

labelList = ["Brag-", "Ideal-"]
data = {}
time_interest = 0.45

for i in range(len(directories)):
  # get a list of all the files in this directory
  direc = directories[i]
  files = get_files(direc, include=["plt"], exclude=["chk"], get_all=False)

  inputs = {}; inputs['t_interest'] = time_interest;
  inputs['fileNames'] = files
  print("Finding time step...")
  stepCerberus, dummy = exactTimeFrameCerberus(inputs);
  #stepCerberus = 100
  print("...found time step", stepCerberus, ".")
  # get all the names of the different states
  #f = files[0]
  #ds = ReadBoxLib(f)
  names = namesList[i]; #["ions", "electrons"] #sorted(ds.names)
  label = labelList[i];
  n_names = len(names)
  n_times = len(files)
  for name in names:
      data[label + name] = {}
      for prop in properties:
          data[label + name][prop] = []
  data[label+"t"] = []
  print("Reading ", label+"data" )
  print("\t", n_names, " names\t", n_times, " n_times")

  fileRead = files[stepCerberus]
  ds = ReadBoxLib(fileRead, max_level = -1)
  
  data[label+"t"] = ds.time
  print("\ttime = ", ds.time)
  for name in names:
      for prop in properties:
          x, v = ds.get(prop + name)
          data[label+"x"] = x[0]
          data[label+name][prop] = v

fig = plt.figure(figsize=(10,10))
ax = {}

for i in range(len(properties)):
  ax[i] = fig.add_subplot(2,2,i+1)
  ax[i].set_ylabel(properties[i][:-1])
  ax[i].set_xlabel("x")
  ax[i].set_xlim([-2, 2])
for i in range(len(directories)):
    print(i )
    label = labelList[i]; names = namesList[i]; #["ions", "electrons"] #sorted(ds.names)
    print(label)
    print(names)
    for name in names:
        print("\t", name)
        for j in range(len(properties)):
            print("\t\t",j )
            ax[j].plot(data[label+"x"], data[label+name][properties[j]], label=label+name)
            #linestyle = "None", marker='x', markersize=2, label=label+name)

for i in range(len(properties)):
  ax[i].legend()

fig.tight_layout()
fig.savefig("Riemann-2F-FullBraginskii-Ideal-Comparison-t=%f.png"%time_interest, dpi=300)
plt.close(fig)
    
    
print("DONE")
