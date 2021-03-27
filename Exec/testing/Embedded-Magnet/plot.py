
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

stuff = {}

# get a list of all the files in this directory
files = get_files('.', include=['plt'])

# get data

f = files[-1]

data = ReadBoxLib(f)
t = data.time
    
data = ReadBoxLib(f, max_level=-1)

x, rho = data.get("rho-ion", grid="node")
xc, Bx = data.get("x_B-field", grid="cell")
xc, By = data.get("y_B-field", grid="cell")
xc, vfi = data.get("vfrac-ion", grid="cell")
xc, vff = data.get("vfrac-field", grid="cell")

y, x = np.meshgrid(x[1], x[0])
yc, xc = np.meshgrid(xc[1], xc[0])

# plot stuff

fig = plt.figure()
ax = fig.add_subplot(111)
pc = ax.pcolormesh(x, y, rho)
ax.streamplot(xc.T, yc.T, Bx.T, By.T, linewidth=0.5, color='black', arrowsize=0.5)

for vf in [vfi, vff]:
    cs = ax.contour(xc, yc, vf, levels=[0.5], colors=['k'])
    for line in cs.allsegs[0]:
        plt.fill(line[:,0], line[:,1])
        
plt.colorbar(pc, label=r"$\rho_e$", orientation='horizontal')

ax.set_aspect(1)
fig.savefig("plot.png", dpi=300)
# plt.show()
    
    
print("DONE")
