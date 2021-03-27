
import sys
cmd_folder = "../../../vis"
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
    
from get_boxlib import ReadBoxLib, get_files

import numpy as np
import pylab as plt
from matplotlib.image import NonUniformImage
from multiprocessing import Pool


#==============================================================================
# 
#==============================================================================

plt_file = "Carbuncle" #str(sys.argv[1])

# get a list of all the files in this directory
files = get_files('.', include=[plt_file], exclude=["temp", ".png"], times=[], tol=1e-4, get_all=True)
    
rh5 = ReadBoxLib(files[-1], max_level=-1)

component = rh5.names[-1]

x, z = rh5.get("rho-%s"%component)

x, y = x

##
# get a window where stuff is actually happening

di, dj = np.gradient(z)
d = np.sqrt(di**2 + dj**2)
d /= d.max()
di = d.sum(axis=1)
dj = d.sum(axis=0)

lo = 5e-2
lo_i = np.argwhere(di>lo)[0]
hi_i = di.size - np.argwhere(di[::-1]>lo)[0] - 1 

lo_j = np.argwhere(dj>lo)[0]
hi_j = dj.size - np.argwhere(dj[::-1]>lo)[0] - 1

window = [[x[lo_i], x[hi_i]], [y[lo_j], y[hi_j]]]

# =============================================================================
# 
# =============================================================================


fig = plt.figure(figsize=(8,4))
ax = fig.add_subplot(111)


# plot the contour
im = NonUniformImage(ax, interpolation='bilinear', extent=np.ravel(window),
                     cmap="viridis")
im.set_data(x, y, z.T)
ax.images.append(im)
plt.colorbar(im, label=r"$\rho_\mathrm{%s}$"%component, orientation="horizontal")

ax.set_xlim(window[0][0], window[0][1])
ax.set_ylim(window[1][0], window[1][1])
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")
ax.set_aspect(1)

fig.tight_layout()

fig.savefig(plt_file+".png", dpi=300)
plt.close(fig)
    
    
print("DONE")
