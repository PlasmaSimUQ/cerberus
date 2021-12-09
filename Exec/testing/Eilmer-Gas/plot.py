
import sys
cmd_folder = "../../../vis"
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
    
from get_boxlib import ReadBoxLib, get_files

import numpy as np
import pylab as plt
from cycler import cycler
    
# =============================================================================
# 
# =============================================================================

linestyle_cycler = cycler('linestyle',['-','--',':','-.'])
plt.rc('axes', prop_cycle=linestyle_cycler)

# get a list of all the files in this directory
files = get_files('.', include=['plt'], exclude=['old'], get_all=True)

data = {}

for f in files:
    ds = ReadBoxLib(f)

    names = ds.names

    for sid, name in enumerate(names):

        comp_names = ds.data["hydro_comp_names"][sid]

        sum = 0.0

        n_comp = len(comp_names)

        x, rho = ds.get("rho-%s"%(name))

        for cid, comp in enumerate(comp_names):

            comp = name+"_"+comp

            if comp not in data:
                data[comp] = {"val":[], "t":[]}

            dat = data[comp]

            if n_comp == 1:

                dat["val"].append(rho[0,0])
                dat["t"].append(ds.time)

            else:

                if (cid < len(comp_names)-1):        
                    x, m = ds.get("tracer_%i-%s"%(cid,name))
                    sum += m

                else:
                    m = rho - sum
                
                dat["val"].append(m[0,0])
                dat["t"].append(ds.time)
        
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1,1,1)

for name, dat in data.items():
    ax.semilogy(dat["t"], dat["val"], label=name)

ax.legend()

fig.tight_layout()
fig.savefig("plot.png", dpi=300)
plt.close(fig)
    
    
print("DONE")
