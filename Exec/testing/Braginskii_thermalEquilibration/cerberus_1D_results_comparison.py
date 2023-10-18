
import sys
cmd_folder = "/home/kyriakos/Documents_ubuntu/000_refactor_cerberus/cerberus/vis"
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
    
from get_boxlib import ReadBoxLib, get_files

import numpy as np
import pylab as plt
import pdb   
# =============================================================================
# 
# =============================================================================
properties = ["rho-", "p-", "T-","x_vel-", "y_vel-", "z_vel-"]

plotProperties = []

directories = [
# braginskii  
"/home/kyriakos/Documents_ubuntu/000_refactor_cerberus/cerberus/Exec/testing/Braginskii_thermalEquilibration/", 
#Ghosh

    ]

namesList = [["ion", "electron"]]#, ["A", "B"], ["ion", "electron"], ["ion", "electron"]]
nameStyle = {"ion":"dashed", "electron":"dotted"}#, "A":"dashed", "B":"dotted"}

#namesList = [["A", "B"], ["A", "B"]]
#nameStyle = {"ion":"dashed", "electron":"dotted", "A":"dashed", "B":"dotted"}

#namesList = [["ion", "electron"], ["ion", "electron"], ["A", "B"], ["A", "B"]]
#nameStyle = {"ion":"dashed", "electron":"dotted", "A":"dashed", "B":"dotted"}

#labelList = ["Brag-", "Brag-QiRegress", "Ghosh-", "Ghosh-noIonRu-"]
labelColour = ["r", "b"]
labelList = ["Brag-", "Ghosh-"]

#labelList = ["Ghosh-", "Ghosh-D"]
#labelColour = ["r", 'm', "b", "k"]

data = {}
for i in range(len(directories)):
  print("for dir ", i)
  # get a list of all the files in this directory
  direc = directories[i]
  files = get_files(direc, include=["plt"], exclude=["chk"], get_all=True)
  
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
  print("Reading ", label+"data")
  print("\t", n_names, " names\t", n_times, " n_times")
  for f in files:
      ds = ReadBoxLib(f)
  
      data[label+"t"].append(ds.time)
      
      for name in names:
          for prop in properties:
              x, v = ds.get(prop + name)
              data[label+name][prop].append(v[0,0])

### total eneergy
gam = 5./3.;
for i in range(len(directories)):
  print(f"energy dir:{i}")
  names = namesList[i]; 
  label = labelList[i];
  data[label] = {}
  aa = np.array(data[label + names[0]]['rho-'])
  data[label]["nrg-total-"] = np.zeros( aa.shape )

  for name in names:
      print(f"\tName:{name}")
      mv2 = np.array(data[label + name]["rho-"])*np.array(data[label + name]["x_vel-"])**2 + \
            np.array(data[label + name]["rho-"])*np.array(data[label + name]["y_vel-"])**2 + \
            np.array(data[label + name]["rho-"])*np.array(data[label + name]["z_vel-"])**2;
  
      data[label + name]["nrg-"] = np.array(data[label + name]["p-"])/(gam - 1.0) + \
                                  mv2/2;
      data[label]["nrg-total-"] += data[label + name]["nrg-"]



fig = plt.figure(figsize=(10,10))
ax = {}

properties[0] = 'nrg-total-'
#properties[0] = 'nrg-'
properties = properties[0:4]
print(properties)
for i in range(len(properties)):
  print(f"prop:{i}")
  ax[i] = fig.add_subplot(2,2,i+1)
  ax[i].set_ylabel(properties[i][:-1])
  ax[i].set_xlabel("time")

for i in range(len(directories)):
    print(f"\ndirec:{i}")
    label = labelList[i]; names = namesList[i]; #["ions", "electrons"] #sorted(ds.names)
    useColour = labelColour[i]
    print(f"\t{label}")
    print(f"\t{names}")
    for name in names:
        useLine = nameStyle[name]
        print( "\t\t", name)
        for j in range(len(properties)):
            print(f"\t\t\tprop:{properties[j]}")
            if properties[j] == 'nrg-total-':
              ax[j].plot(data[label+"t"], data[label][properties[j]], color=useColour, \
                linestyle=useLine, label=label+name) # marker='x', markersize=2,              
            else:
              #print "\t\t",j 
              ax[j].plot(data[label+"t"], data[label+name][properties[j]], color=useColour, \
                linestyle=useLine, label=label+name) # marker='x', markersize=2, 

print("\n\n")
for i in range(len(properties)):
    print(i)
    ax[i].legend()

fig.tight_layout()
fig.savefig("2023_refactor_test_1_1_deleteMe.png", dpi=300)
#fig.savefig("20220812_ThermalEquilibration-Case_2_2_mime_100_EnergyConsCorrected.png", dpi=300)
plt.close(fig)
    
    
print("DONE")
