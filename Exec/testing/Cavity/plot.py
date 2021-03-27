
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

component = "fluid"

# get a list of all the files in this directory
files = get_files('.', include=['plt'], times=[x for x in range(-10,0,1)])# get_all=True)

# get tracer particle data

def get_tracer_data(f):
    ds = ReadBoxLib(f)
    print(f," @ ",ds.time)
    
    idata, rdata = ds.get_particles(component)

    return idata, rdata

nproc = 1

if nproc > 1:
    pool = Pool(processes=nproc)
    DATA = pool.map(get_tracer_data, files)
    pool.close()
    pool.join()
else:
    DATA = []
    for f in files:
        DATA.append(get_tracer_data(f))

# due to the invalidation of particles we have a varying number of particles
# throughout the simulation. Therefore, we need to match up our particles over time
# in order to plot the streak lines

# find the unique ids
idx = set()

for pairs in DATA[0][0]:
    idx.add(tuple(pairs))

traces = {}
for i in idx:
    x = []
    y = []
    

    for d in DATA:
        ii = np.where((d[0][:,0] == i[0]) & (d[0][:,1] == i[1]))
        if len(ii[0]) > 0:
            ii = ii[0][0]
        else:
            continue
        x.append(d[1][ii,0])
        y.append(d[1][ii,1])

    traces[i] = {"x":x, "y":y}


    
# get fluid variables
    
ds = ReadBoxLib(files[-1], max_level=-1)

_, rho = ds.get("rho-%s"%component)
xc, u = ds.get("x_vel-%s"%component)
xn, v = ds.get("y_vel-%s"%component, grid="node")


vel = np.sqrt(u**2 + v**2)

x, y = xc


# =============================================================================
# 
# =============================================================================


fig = plt.figure(figsize=(4,4))
ax = fig.add_subplot(111)


# plot the contour
im = NonUniformImage(ax, interpolation='bilinear', extent=(xn[0][0], xn[0][0], xn[0][-1], xn[1][-1]),
                     cmap="viridis")
im.set_data(x, y, vel.T)
ax.images.append(im)
plt.colorbar(im, label=r"$\left| u_\mathrm{%s}\right|$"%component, orientation="horizontal")

# plot particle traces
for k, t in traces.items():
    ax.plot(t["x"], t["y"],'w-',lw=0.25, alpha=0.5)
    if len(t["x"]) > 1:
        ax.plot([t["x"][-1]],[t["y"][-1]],'wo', ms=0.25, alpha=0.5)

ax.set_xlim(xn[0][0], xn[0][-1])
ax.set_ylim(xn[1][0], xn[1][-1])
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")
ax.set_aspect(1)

fig.tight_layout()

fig.savefig("plot.png", dpi=300)
plt.close(fig)
    
    
print("DONE")
