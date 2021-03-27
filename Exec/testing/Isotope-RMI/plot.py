
import sys
cmd_folder = "../../../vis"
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
    
from get_boxlib import ReadBoxLib, get_files

import numpy as np
import pylab as plt
from matplotlib.image import NonUniformImage


#==============================================================================
# 
#==============================================================================

files = get_files('.', include=['plt'], get_all=True)

f = files[-1]

window = [[-1.25,0], [1.25,1]]

data = ReadBoxLib(f, limits = window)



xc, re = data.get("rho-electron")
xc, me = data.get("mass-electron")
ne = re/me

xn, ri = data.get("rho-ion", grid="node")
xn, mi = data.get("mass-ion", grid="node")
ni = ri/mi

xc, yc = xc
xn, yn = xn

fig = plt.figure(figsize=(6,4))

axes = []

ax = fig.add_subplot(211); axes.append(ax)
im = NonUniformImage(ax, interpolation='bilinear', extent=(xn[0], yn[0], xn[-1], yn[-1]))
im.set_data(xc, yc, ne.T)
ax.images.append(im)
plt.colorbar(im, label=r"$n_e$")

ax = fig.add_subplot(212); axes.append(ax)
im = NonUniformImage(ax, interpolation='bilinear', extent=(xn[0], yn[0], xn[-1], yn[-1]))
im.set_data(xc, yc, ni.T)
ax.images.append(im)
plt.colorbar(im, label=r"$n_i$")

for ax in axes:
    ax.set_xlim(window[0][0], window[1][0])
    ax.set_ylim(window[0][1], window[1][1])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect(1)

fig.tight_layout()

fig.savefig("IRMI.png", dpi=300)
plt.close(fig)
    
    
print("DONE")