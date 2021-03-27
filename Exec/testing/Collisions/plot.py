
import sys
cmd_folder = "../../../vis"
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
    
from get_boxlib import ReadBoxLib, get_files

import numpy as np
import pylab as plt
import h5py
    
# =============================================================================
# 
# =============================================================================

# get a list of all the files in this directory
files = get_files('.', include=['plt'], get_all=True)

# get all the names of the different states
f = files[0]
ds = ReadBoxLib(f)
names = ds.names

data = {}
for name in names:
    data[name] = {"T":[], "t":[], "style":{"ls":"-"}}

for f in files:
    ds = ReadBoxLib(f)

    for name in names:

        dat = data[name]

        x, T = ds.get("T-%s"%name)
        dat["T"].append(T[0])

        dat["t"].append(ds.time)

        

# reference solution
h5 = h5py.File("Ghosh2019_fig_5d.hdf5", "r")

t = []
T = []
for key, line in h5.items():
    if key == "N":
        continue
    x = line['X'][()]
    y = line['Y'][()]

    if len(x) == 1:
        t.append(x[0])
        T.append(y[0])

h5.close()

t = np.array(t)
T = np.array(T)

I = np.argsort(t)

t = t[I]
T = T[I]

data[r"Rambo & Procassini, 1995"] = {"T":T, "t":t, "style":{"ls":'None', "marker":'o', "mfc":"None", "ms":5, "mec":"black"}}

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1,1,1)

for name, dat in data.items():
    ax.plot(dat["t"], dat["T"], label=name, **dat["style"])

ax.legend()

fig.tight_layout()
fig.savefig("plot.png", dpi=300)
plt.close(fig)
    
    
print("DONE")
