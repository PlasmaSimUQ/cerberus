
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
files = get_files('.', include=['plt'])

N = 3

for i in range(N):

    # get data
    f = files[i*int(len(files)/float(N-1))]

    data = ReadBoxLib(f)
    t = data.time
        
    data = ReadBoxLib(f, max_level=-1)

    xc, vf = data.get("vfrac-field", grid="cell")
    xn, xD = data.get("x_D-field", grid="node")

    yn, xn = np.meshgrid(xn[1], xn[0])
    yc, xc = np.meshgrid(xc[1], xc[0])

    # plot stuff


    ax = fig.add_subplot(1,N,i+1)
    pc = ax.pcolormesh(xn, yn, xD)
    ax.contour(xc, yc, vf, [0.5], colors=['k'], linewidths=[0.5])
    plt.colorbar(pc, orientation='horizontal')

    ax.set_aspect(1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(r"$D_x$ @ $t=%g$"%t)

fig.tight_layout()

fig.savefig("plot.png", dpi=300)
# plt.show()
    
    
print("DONE")
