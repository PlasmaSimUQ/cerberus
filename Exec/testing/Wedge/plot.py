
import sys
cmd_folder = "../../../vis"
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
    
from get_boxlib import ReadBoxLib, get_files

import numpy as np
import pylab as plt
from matplotlib.image import NonUniformImage
import matplotlib.ticker as ticker


#==============================================================================
# 
#==============================================================================

fig = plt.figure(figsize=(10,5))

# get a list of all the files in this directory
files = get_files('.', include=['plt'], times=[0.001, 0.002, 0.0301], tol=1e-2)
N = len(files)

for i, f in enumerate(files):

    data = ReadBoxLib(f)
    t = data.time
        
    data = ReadBoxLib(f, max_level=-1)

    xc, vf = data.get("vfrac-fluid", grid="cell")
    xn, r = data.get("rho-fluid", grid="node")

    yn, xn = np.meshgrid(xn[1], xn[0])
    yc, xc = np.meshgrid(xc[1], xc[0])

    # plot stuff


    ax = fig.add_subplot(1,N,i+1)
    pc = ax.pcolormesh(xn, yn, r)
    ax.contour(xc, yc, vf, [0.5], colors=['k'], linewidths=[0.5])
    plt.colorbar(pc, orientation='horizontal')

    ax.set_aspect(1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(r"$\rho$ @ $t=%g$"%t)

fig.tight_layout()

fig.savefig("plot.png", dpi=300)
# plt.show()
    
    
print("DONE")
